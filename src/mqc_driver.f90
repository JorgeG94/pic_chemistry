!! Main calculation driver module for metalquicha
module mqc_driver
   !! Routes a configured system to the fragmented (many-body expansion) or
   !! unfragmented path, and holds the run types that take neither: MAKEFP, EFP
   !! and SAPT.
   use pic_types, only: int32, int64, dp, default_int
   use pic_mpi_lib, only: comm_t, abort_comm, bcast, allgather
   use mqc_resources, only: resources_t
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use mqc_method_types, only: needs_serial_execution
   use mqc_mbe_fragment_distribution_scheme, only: unfragmented_calculation, distributed_unfragmented_hessian
   use mqc_many_body_expansion, only: many_body_expansion_t, mbe_context_t, gmbe_context_t, &
                                      fmo_context_t
   use mqc_method_config, only: method_config_t
   ! GMBE functions are now called via type-bound procedures in gmbe_context_t
   use mqc_validate, only: validate_system, validate_terms
   use mqc_fraglist, only: fraglist_t
   use mqc_frag_utils, only: get_nfrags, create_monomer_list, generate_fragment_list, generate_intersections, &
                             gmbe_enumerate_pie_terms, binomial, combine, apply_distance_screening, &
                             sort_fragments_by_size, generate_mbe_term_list
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, &
                                    check_system_geometry, &
                                    build_fragment_from_indices, build_fragment_from_atom_list
   use mqc_config_adapter, only: driver_config_t, config_to_driver, config_to_system_geometry, &
                                 check_counterpoise_support
   use mqc_method_types, only: method_type_to_string
   use mqc_calc_types, only: calc_type_to_string, CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, &
                             CALC_TYPE_OPTIMIZE, CALC_TYPE_CONFORMERS, &
                             CALC_TYPE_HESSIAN, CALC_TYPE_MAKEFP
   use mqc_config_types, only: bond_t, mqc_config_t
   use mqc_scf_types, only: scf_numerics_t
   use mqc_mbe, only: compute_gmbe
   use mqc_result_types, only: calculation_result_t
   use mqc_scf_common, only: lindep_tally_t, lindep_collect_begin, lindep_collect_end, &
                             report_linear_dependence_tally
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_bridge, only: czt_set_eri_path
   use mqc_fingerprint, only: calculation_fingerprint
   use mqc_checkpoint, only: checkpoint_t
   use mqc_io_helpers, only: set_molecule_suffix, get_output_json_filename
   use mqc_json, only: merge_multi_molecule_json
   use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_NONE
   use mqc_json_writer, only: write_json_output
   implicit none
   private

   public :: run_calculation  !! Main entry point for all calculations
   public :: run_multi_molecule_calculations  !! Multi-molecule calculation dispatcher

contains

   subroutine run_calculation(resources, config, sys_geom, bonds, result_out, all_ranks_write_json, &
                              supplied_terms, n_supplied_terms, write_output)
      !! Dispatch a configured system to the calculation it asks for
      !!
      !! `config%nlevel` chooses fragmented or unfragmented; MAKEFP, EFP and
      !! SAPT take neither path. `result_out` and `write_output` are
      !! independent: a caller can take the energy back, write the files, or
      !! both.
      !!
      !! `supplied_terms` hands over a term list instead of generating one.
      !! Distance screening and sorting are skipped with it: the list is taken
      !! as given, in the order given. Only rank 0 needs it -- workers build
      !! each fragment from `sys_geom`, which every rank must already have.
      !!
      !! **The list must be closed under subsets.** An n-body term's delta is
      !! its energy less every proper subset's delta, so keeping a trimer whose
      !! dimers were screened away fails the lookup rather than approximating
      !! anything. `fraglist_t%close_subsets` exists for this.
      use mqc_method_types, only: METHOD_TYPE_EFP2, METHOD_TYPE_SAPT0, &
                                  METHOD_TYPE_SAPT2
      type(resources_t), intent(in) :: resources
      !! Resources container (MPI comms, etc.)
      type(driver_config_t), intent(in) :: config
      !! Driver configuration
      type(system_geometry_t), intent(in) :: sys_geom
      !! System geometry and fragment info
      type(bond_t), intent(in), optional :: bonds(:)
      !! Bond connectivity information
      type(calculation_result_t), intent(out), optional :: result_out
      !! Optional result output
      logical, intent(in), optional :: all_ranks_write_json
      !! If true, all ranks write JSON (for multi-molecule)
      logical, intent(in), optional :: write_output
         !! Write the JSON summary and fragment breakdown. Default true.
         !! `skip_json_output` in the config still overrides it.
      integer, intent(in), optional :: supplied_terms(:, :)
         !! (n_terms, max_level), 1-based monomer indices, zero-padded. Rank 0 only.
      integer(int64), intent(in), optional :: n_supplied_terms

      ! Local variables
      integer :: max_level   !! Maximum fragment level (nlevel from config)
      integer :: i  !! Loop counter
      type(json_output_data_t) :: json_data  !! Cached output data for centralized JSON writing
      logical :: should_write_json  !! Whether this rank should write JSON
      logical :: wants_output       !! Whether the caller asked for files at all
      type(error_t) :: geometry_error  !! Two atoms in the same place, if any
      type(lindep_tally_t) :: lindep   !! One report for every fragment SCF

      ! TODO(mqc): refactor

      ! Set max_level from config
      max_level = config%nlevel

      ! Log method-specific settings (rank 0 only)
      if (resources%mpi_comms%world_comm%rank() == 0) then
         call config%method_config%log_settings()
      end if

      ! The four-centre integral path, chosen once for every quartet of the
      ! run and on every rank, since the workers build their own Fock
      ! matrices. Refused here, before any integral, if the build lacks it.
      block
         type(error_t) :: path_error
         call czt_set_eri_path(config%method_config%scf%eri_path, path_error)
         if (path_error%has_error()) then
            if (resources%mpi_comms%world_comm%rank() == 0) then
               call logger%error(path_error%get_message())
            end if
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if
      end block

      if (resources%mpi_comms%world_comm%rank() == 0 .and. max_level > 0) then
         call logger%info("============================================")
         call logger%info("Loaded geometry:")
         call logger%info("  Total monomers: "//to_char(sys_geom%n_monomers))
         ! Zero is the "variable-sized fragments" sentinel, not a count -- say so
         ! rather than printing "0 atoms per monomer", which reads as nonsense.
         if (sys_geom%atoms_per_monomer > 0) then
            call logger%info("  Atoms per monomer: "//to_char(sys_geom%atoms_per_monomer))
         else
            call logger%info("  Atoms per monomer: variable")
         end if
         call logger%info("  Fragment level: "//to_char(max_level))
         call logger%info("  Total atoms: "//to_char(sys_geom%total_atoms))
         call logger%info("============================================")
      end if

      ! Before anything is dispatched, and on every path. Two atoms in the same
      ! place make the electron count wrong, and the per-fragment check that
      ! would eventually notice does so a long way into a large expansion, or
      ! never at all if screening drops the fragment holding both copies.
      call check_system_geometry(sys_geom, geometry_error)
      if (geometry_error%has_error()) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%error(geometry_error%get_message())
         end if
         ! `abort_comm` rather than `error stop`: under MPI the latter reaches
         ! MPI_ABORT before anything written above it is flushed, so the run
         ! dies with a rank number and no reason.
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      ! An optimization is a loop over calculations and is driven from above
      ! this routine, so reaching here with one means a caller took a path that
      ! does not know about it -- the multi-molecule loop, a session, the C API.
      ! Refused by name: left to fall through, `calc_type` reaches the method
      ! layer as an unhandled value and the run dies in a backtrace.
      if (config%calc_type == CALC_TYPE_CONFORMERS) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%error('driver "conformers" is not available on this path.')
            call logger%error("  Sampling drives run_calculation rather than being "// &
                              "driven by it, so it")
            call logger%error("  works for a single molecule from an input deck, and "// &
                              "not from a")
            call logger%error("  multi-molecule deck, a session or the C API.")
         end if
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      if (config%calc_type == CALC_TYPE_OPTIMIZE) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%error('driver "Optimize" is not available on this path.')
            call logger%error("  It is driven above run_calculation, so it works for a "// &
                              "single molecule")
            call logger%error("  from an input deck, and not yet from a multi-molecule "// &
                              "deck, a session")
            call logger%error("  or the C API.")
         end if
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      ! Resolved before the branches below rather than beside the JSON write at
      ! the end, because those branches write their own summary and return
      ! without ever reaching it.
      wants_output = .true.
      if (present(write_output)) wants_output = write_output

      ! MAKEFP writes a file and returns no energy, so it takes neither the
      ! fragmented nor the unfragmented path: a fragment potential is built from
      ! the whole system by definition, and there is no result_t to fill.
      if (config%calc_type == CALC_TYPE_MAKEFP) then
         call run_makefp(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                         result_out)
         return
      end if

      ! Quantum nuclei take neither path either: a proton's orbital is solved
      ! with every electron of the system, so there is nothing to fragment.
      if (config%method_config%neo%active) then
         call run_neo(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                      wants_output, result_out)
         return
      end if

      ! EFP takes neither path either: the fragments already carry their
      ! wavefunctions, so there is no SCF to fragment and nothing for a
      ! many-body expansion to expand.
      if (config%method_config%method_type == METHOD_TYPE_EFP2) then
         call run_efp(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                      wants_output, result_out)
         return
      end if

      ! EFMO takes neither path either. It is a fragmented method, but not a
      ! many-body expansion: there is no term list, no subset closure and no
      ! embedding loop, only monomers, a pair split and one induction over
      ! everything. The whole expression is assembled inside the backend.
      if (config%expansion_kind == "efmo") then
         call run_efmo_energy(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                              wants_output, result_out)
         return
      end if

      ! SAPT takes neither path either: it returns the interaction between two
      ! monomers rather than the energy of one system.
      if (config%method_config%method_type == METHOD_TYPE_SAPT0 .or. &
          config%method_config%method_type == METHOD_TYPE_SAPT2) then
         call run_sapt(config, sys_geom, resources%mpi_comms%world_comm%rank(), &
                       wants_output, result_out)
         return
      end if

      ! Warn if overlapping fragments flag is set but nlevel=0
      if (config%allow_overlapping_fragments .and. max_level == 0) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%warning("allow_overlapping_fragments is set to true, but nlevel=0")
            call logger%warning("Running unfragmented calculation - overlapping fragments flag will be ignored")
         end if
      end if

      if (max_level == 0) then
         ! One thread, for the methods that need it rather than for all of them.
         ! tblite run threaded corrupts a result instead of failing, so xTB is
         ! clamped here; the ab initio path keeps the threads it was given.
         !
         ! The clamp is not restored afterwards. `omp_set_num_threads(1)` makes
         ! `omp_get_max_threads()` report 1, so the value has to be saved before
         ! the call to be recoverable; an unfragmented run ends here, and the
         ! one path that does need the count back saves it first. See
         ! mqc_serial_fragment_processor.
         if (needs_serial_execution(config%method_config%method_type)) then
            call omp_set_num_threads(1)
         end if
         if (present(result_out)) then
            call run_unfragmented_calculation(resources%mpi_comms%world_comm, sys_geom, config, &
                                              result_out, json_data=json_data)
         else
            call run_unfragmented_calculation(resources%mpi_comms%world_comm, sys_geom, config, &
                                              json_data=json_data)
         end if
      else
         ! json_data is collected whether or not it will be written, because it
         ! is also where a fragmented result comes from -- the expansion has no
         ! other route back to a calculation_result_t.
         !
         ! Fold every fragment SCF's linear-dependence report into one, rather
         ! than a block per fragment. Unfragmented runs are left alone: there is
         ! one SCF, so the per-SCF report is already the summary. The window is
         ! per rank, so with MPI each rank reports its own fragments.
         call lindep_collect_begin()
         call run_fragmented_calculation(resources, config, sys_geom, bonds, json_data, &
                                         supplied_terms=supplied_terms, &
                                         n_supplied_terms=n_supplied_terms)
         call lindep_collect_end(lindep)
         call report_linear_dependence_tally(lindep, "fragment SCFs")
         if (present(result_out)) call result_from_json(json_data, result_out)
      end if

      ! Stamped whether or not it is written: a caller taking the energy back
      ! without files still needs to know what produced it.
      json_data%fingerprint = calculation_fingerprint(sys_geom, config%method_config, &
                                                      config%calc_type)

      ! Centralized JSON output (rank 0 only by default, or all ranks if all_ranks_write_json is set)
      if (wants_output) then
         ! Check if JSON output should be skipped
         if (config%skip_json_output) then
            if (resources%mpi_comms%world_comm%rank() == 0) then
               call logger%large_info("Skipping JSON output (skip_json_output = true)")
            end if
         else
            ! Determine if this rank should write JSON
            should_write_json = (resources%mpi_comms%world_comm%rank() == 0)
            if (present(all_ranks_write_json)) then
               if (all_ranks_write_json) should_write_json = .true.
            end if

            if (should_write_json) then
               if (json_data%output_mode /= OUTPUT_MODE_NONE) then
                  json_data%fragment_breakdown = config%fragment_breakdown
                  call write_json_output(json_data)
                  call json_data%destroy()
               end if
            end if
         end if
      end if

   end subroutine run_calculation

   subroutine result_from_json(json_data, result_out)
      !! Take a fragmented run's headline numbers back to the caller
      !!
      !! The total, and the gradient and dipole when the run produced them. The
      !! per-fragment detail stays in the CSV.
      type(json_output_data_t), intent(in) :: json_data
      type(calculation_result_t), intent(inout) :: result_out

      ! energy_t computes its total from components; an expansion total has
      ! no correlation breakdown to put in the others, so it lands in scf --
      ! which is what `total()` will then report.
      result_out%energy%scf = json_data%total_energy
      result_out%has_energy = json_data%has_energy
      if (allocated(json_data%gradient)) then
         result_out%gradient = json_data%gradient
         result_out%has_gradient = .true.
      end if
      if (allocated(json_data%dipole)) result_out%dipole = json_data%dipole
   end subroutine result_from_json

   subroutine run_unfragmented_calculation(world_comm, sys_geom, config, result_out, json_data)
      !! Run the whole system as one calculation (`nlevel = 0`)
      !!
      !! Rank 0 alone for a single molecule; every rank for a multi-molecule
      !! run, each with its own. A Hessian across more than one rank goes to
      !! `distributed_unfragmented_hessian`, which spreads the displacements.
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(system_geometry_t), intent(in) :: sys_geom  !! Complete system geometry
      type(driver_config_t), intent(in) :: config  !! Driver configuration (includes method_config, calc_type, etc.)
      type(calculation_result_t), intent(out), optional :: result_out  !! Optional result output
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      ! For Hessian calculations with multiple ranks, use distributed approach
      if (config%calc_type == CALC_TYPE_HESSIAN .and. world_comm%size() > 1) then
         if (world_comm%rank() == 0) then
            call logger%large_info(" ")
            call logger%large_info("Running distributed unfragmented Hessian calculation")
            call logger%large_info("  MPI ranks: "//to_char(world_comm%size()))
            call logger%large_info(" ")
         end if
         call distributed_unfragmented_hessian(world_comm, sys_geom, config, json_data)
         return
      end if

      ! Check if this is multi-molecule mode or single-molecule mode
      ! In multi-molecule mode, each rank processes its own molecule
      ! In single-molecule mode, only rank 0 should work
      if (world_comm%size() == 1 .or. world_comm%rank() == 0) then
         ! Either single-rank calculation, or rank 0 in multi-rank setup
         call logger%info(" ")
         call logger%info("Running unfragmented calculation")
         call logger%info("  Calculation type: "//calc_type_to_string(config%calc_type))
         call logger%info(" ")
         call unfragmented_calculation(sys_geom, config, result_out, json_data)
      else if (sys_geom%total_atoms > 0) then
         ! Multi-molecule mode: non-zero rank with a molecule
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Running unfragmented calculation")
         call unfragmented_calculation(sys_geom, config, result_out, json_data)
      end if

   end subroutine run_unfragmented_calculation

   subroutine run_fragmented_calculation(resources, config, sys_geom, bonds, json_data, &
                                         supplied_terms, n_supplied_terms)
      !! Run a many-body expansion (`nlevel > 0`)
      !!
      !! Builds the term list, groups the ranks into nodes and groups, and hands
      !! the work to an expansion context. `allow_overlapping_fragments` selects
      !! GMBE with its inclusion-exclusion terms instead of plain MBE.

      type(resources_t), intent(in), target :: resources  !! Resources container (MPI comms, etc.)
      type(driver_config_t), intent(in) :: config  !! Driver configuration (includes method_config, calc_type, etc.)
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      type(bond_t), intent(in), optional :: bonds(:)  !! Bond connectivity information
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data
      integer, intent(in), optional :: supplied_terms(:, :)
         !! (n_terms, max_level), 1-based monomer indices, zero-padded. Rank 0 only.
      integer(int64), intent(in), optional :: n_supplied_terms

      ! Local variables extracted from config for readability
      integer :: max_level    !! Maximum fragment level for MBE
      logical :: allow_overlapping_fragments  !! Use GMBE for overlapping fragments
      integer :: max_intersection_level  !! Maximum k-way intersection depth for GMBE

      integer(int64) :: total_fragments  !! Total number of fragments generated (int64 to handle large systems)
      integer(default_int) :: supplied_width  !! Columns the caller actually provided
      type(error_t) :: checkpoint_error
      integer, allocatable :: polymers(:, :)  !! Fragment composition array (fragment, monomer_indices)
      integer :: num_nodes   !! Number of compute nodes
      integer :: i, j        !! Loop counters
      integer, allocatable :: node_leader_ranks(:)  !! Ranks of processes that lead each node
      integer :: global_groups  !! Number of global groups for multi-global coordinator
      integer :: nodes_per_group  !! Nodes per group (optional override)
      integer :: group_size  !! Nodes per group after resolving config
      integer :: group_id, group_start
      integer, allocatable :: group_leader_ranks(:)  !! Group leader rank for each node leader
      integer, allocatable :: group_ids(:)  !! Group id for each node leader
      integer, allocatable :: monomers(:)     !! Temporary monomer list for fragment generation
      integer(int64) :: n_expected_frags  !! Fragment count the combinatorics predict
      integer(int64) :: n_rows      !! Number of rows needed for polymers array (int64 to handle large systems)
      integer :: global_node_rank  !! Global rank if this process leads a node, -1 otherwise
      integer, allocatable :: all_node_leader_ranks(:)  !! Node leader status for all ranks

      ! Polymorphic expansion context for unified MBE/GMBE handling
      class(many_body_expansion_t), allocatable :: expansion

      ! GMBE PIE-based variables
      integer :: n_primaries  !! Number of primary polymers
      integer(int64) :: n_primaries_i64  !! For binomial calculation
      integer, allocatable :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, allocatable :: pie_coefficients(:)  !! PIE coefficient for each term
      integer(int64) :: n_pie_terms  !! Number of unique PIE terms
      type(error_t) :: pie_error  !! Error from PIE enumeration
      type(error_t) :: validation_error  !! Error from semantic validation
      type(fraglist_t) :: supplied_check  !! Supplied terms, for validation only

      ! Extract values from config for readability
      max_level = config%nlevel
      allow_overlapping_fragments = config%allow_overlapping_fragments
      max_intersection_level = config%max_intersection_level
      ! Enumerated on the coordinator alone, far below, but copied into the
      ! expansion on every rank -- so the ranks that never enumerate need a
      ! value to copy rather than whatever the stack happened to hold.
      n_pie_terms = 0_int64

      ! Every input path arrives here, so this is where the system is checked:
      ! the JSON reader, the C interface and a supplied term list alike.
      if (resources%mpi_comms%world_comm%rank() == 0) then
         call validate_system(sys_geom,.not. config%unchecked_input, validation_error, &
                              check_bonds=allocated(sys_geom%bonds))
         if (validation_error%has_error()) then
            call logger%error("invalid system: "//validation_error%get_message())
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if

         ! Ahead of the branch, because each of the three expansions below
         ! ignores counterpoise in its own way and none of them says so.
         call validation_error%clear()
         call check_counterpoise_support(config, validation_error)
         if (validation_error%has_error()) then
            call logger%error(validation_error%get_message())
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if
      end if

      ! Generate fragments -- unless the caller brought their own
      if (present(supplied_terms) .and. present(n_supplied_terms)) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            total_fragments = n_supplied_terms
            allocate (polymers(max(total_fragments, 1_int64), max_level))
            polymers = 0
            ! The caller's array is as wide as their highest term, which need
            ! not be `max_level`: a screen that keeps no trimers hands over a
            ! two-column list for a level-3 expansion, and copying `max_level`
            ! columns from it reads off the end.
            supplied_width = int(size(supplied_terms, 2), default_int)
            if (supplied_width > max_level) then
               ! Wider than the expansion allows. Truncating would silently
               ! turn an n-mer into a smaller one, so refuse instead -- but
               ! only if the extra columns actually hold anything.
               if (any(supplied_terms(1:total_fragments, max_level + 1:supplied_width) /= 0)) then
                  call logger%error("invalid fragment list: a term names more than "// &
                                    to_char(max_level)//" monomers, which is the level "// &
                                    "this expansion was configured for")
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if
               supplied_width = max_level
            end if
            if (total_fragments > 0) then
               polymers(1:total_fragments, 1:supplied_width) = &
                  supplied_terms(1:total_fragments, 1:supplied_width)
            end if
            ! A supplied list comes from outside unscreened and unsorted, so it
            ! is checked before anything is spent on it -- above all for subset
            ! closure, which a reasonable-looking screen breaks silently.
            call supplied_check%replace(polymers, total_fragments, max_level, validation_error)
            if (.not. validation_error%has_error()) then
               call validate_terms(supplied_check, sys_geom,.not. config%unchecked_input, &
                                   validation_error)
            end if
            call supplied_check%destroy()
            if (validation_error%has_error()) then
               call logger%error("invalid fragment list: "//validation_error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            call logger%info("Using a supplied fragment list:")
            call logger%info("  Total fragments: "//to_char(total_fragments))
            call logger%info("  Max level: "//to_char(max_level))
         end if
      else if (resources%mpi_comms%world_comm%rank() == 0) then
         if (allow_overlapping_fragments) then
            ! GMBE mode: PIE-based inclusion-exclusion
            ! GMBE(1): primaries are monomers
            ! GMBE(N): primaries are N-mers (e.g., dimers for N=2)

            ! Generate primaries
            if (max_level == 1) then
               ! GMBE(1): primaries are base monomers
               n_primaries = sys_geom%n_monomers
               allocate (polymers(n_primaries, 1))
               do i = 1, n_primaries
                  polymers(i, 1) = i
               end do
            else
               ! GMBE(N): primaries are all C(M, N) N-tuples
               n_primaries_i64 = binomial(sys_geom%n_monomers, max_level)
               n_primaries = int(n_primaries_i64)
               allocate (monomers(sys_geom%n_monomers))
               allocate (polymers(n_primaries, max_level))
               polymers = 0

               call create_monomer_list(monomers)
               total_fragments = 0_int64
               call combine(monomers, sys_geom%n_monomers, max_level, polymers, total_fragments)
               n_primaries = int(total_fragments)
               deallocate (monomers)

               ! Apply distance-based screening to primaries if cutoffs are provided
               if (max_level > 1) then
                  ! Only screen if primaries are n-mers (not for GMBE(1) where primaries are monomers)
                  total_fragments = int(n_primaries, int64)
                  call apply_distance_screening(polymers, total_fragments, sys_geom, config, max_level)
                  n_primaries = int(total_fragments)
               end if

               ! Sort primaries by size (largest first)
               ! TODO(mqc): with the line below commented out, `total_fragments`
               ! is whatever screening left it. That equals `n_primaries` only
               ! because screening always runs on this branch, so the sort is
               ! correct by accident rather than by construction.
               ! total_fragments = int(n_primaries, int64)
               call sort_fragments_by_size(polymers, total_fragments, max_level)
            end if

         call logger%large_info("Generated "//to_char(n_primaries)//" primary "//to_char(max_level)//"-mers for GMBE("// &
                                   to_char(max_level)//")")

            ! Use DFS to enumerate PIE terms with coefficients
            call gmbe_enumerate_pie_terms(sys_geom, polymers, n_primaries, max_level, max_intersection_level, &
                                          pie_atom_sets, pie_coefficients, n_pie_terms, pie_error)
            if (pie_error%has_error()) then
               call logger%error("GMBE PIE enumeration failed: "//pie_error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

         call logger%large_info("GMBE PIE enumeration complete: "//to_char(n_pie_terms)//" unique subsystems to evaluate")

            ! For now: total_fragments = n_pie_terms (each PIE term is a subsystem to evaluate)
            total_fragments = n_pie_terms
         else
            ! Standard MBE mode. Monomers, then n-mers, then screening, then
            ! the size sort, all of it in `generate_mbe_term_list` -- which a
            ! geometry optimization also calls, so the list it freezes is the
            ! one this would have built.
            call generate_mbe_term_list(sys_geom, config, max_level, polymers, total_fragments)

            call logger%info("Generated fragments:")
            call logger%info("  Total fragments: "//to_char(total_fragments))
            call logger%info("  Max level: "//to_char(max_level))
         end if
      end if

      ! Broadcast total_fragments to all ranks
      call bcast(resources%mpi_comms%world_comm, total_fragments, 1, 0)

      ! Determine node leaders
      global_node_rank = -1
      if (resources%mpi_comms%node_comm%rank() == 0) global_node_rank = resources%mpi_comms%world_comm%rank()

      allocate (all_node_leader_ranks(resources%mpi_comms%world_comm%size()))
      call allgather(resources%mpi_comms%world_comm, global_node_rank, all_node_leader_ranks)

      num_nodes = count(all_node_leader_ranks /= -1)

      if (resources%mpi_comms%world_comm%rank() == 0) then
         call logger%large_info("Running with "//to_char(num_nodes)//" node(s)")
      end if

      allocate (node_leader_ranks(num_nodes))
      i = 0
      do j = 1, resources%mpi_comms%world_comm%size()
         if (all_node_leader_ranks(j) /= -1) then
            i = i + 1
            node_leader_ranks(i) = all_node_leader_ranks(j)
         end if
      end do
      deallocate (all_node_leader_ranks)

      global_groups = config%global_groups
      nodes_per_group = config%nodes_per_group

      if (global_groups > 0 .and. nodes_per_group > 0) then
         if (resources%mpi_comms%world_comm%rank() == 0) then
            call logger%error("Both global_groups and nodes_per_group are set; only one may be specified")
         end if
         call abort_comm(resources%mpi_comms%world_comm, 1)
      end if

      if (nodes_per_group > 0) then
         group_size = nodes_per_group
         global_groups = (num_nodes + group_size - 1)/group_size
      else if (global_groups > 0) then
         if (global_groups > num_nodes) global_groups = num_nodes
         group_size = (num_nodes + global_groups - 1)/global_groups
      else
         global_groups = 1
         group_size = num_nodes
      end if

      allocate (group_leader_ranks(num_nodes))
      allocate (group_ids(num_nodes))
      do i = 1, num_nodes
         group_id = min((i - 1)/group_size + 1, global_groups)
         group_ids(i) = group_id
         group_start = (group_id - 1)*group_size + 1
         if (group_start > num_nodes) group_start = num_nodes
         group_leader_ranks(i) = node_leader_ranks(group_start)
      end do

      if (resources%mpi_comms%world_comm%rank() == 0 .and. num_nodes > 1) then
         call logger%large_info("Multi-global groups: "//to_char(global_groups)//" (nodes_per_group="// &
                                to_char(nodes_per_group)//")")
      end if

      ! Build polymorphic expansion context
      if (config%expansion_kind == "fmo" .or. config%expansion_kind == "ee-mbe") then
         ! FMO or electrostatically embedded MBE. Both are the same machinery,
         ! differing only in what a fragment sees of its neighbours and how the
         ! pieces are added up, so one context serves both.
         allocate (fmo_context_t :: expansion)
         select type (expansion)
         type is (fmo_context_t)
            call expansion%init(config%method_config, config%calc_type)
            allocate (expansion%sys_geom, source=sys_geom)
            if (present(bonds)) then
               if (allocated(expansion%sys_geom%bonds)) deallocate (expansion%sys_geom%bonds)
               allocate (expansion%sys_geom%bonds, source=bonds)
            end if
            call fragment_owner_map(sys_geom, expansion%owner, expansion%n_fragments)
            expansion%basis = config%method_config%basis_set
            expansion%bond_breaking = config%bond_breaking
            expansion%cap_scale = config%cap_scale
            if (config%expansion_kind == "fmo") then
               expansion%esp = "exact"
               expansion%expansion = "fmo"
            else
               expansion%esp = "ptc"
               expansion%expansion = "mbe"
            end if
            ! `embedding` overrides what the expansion implies, which is how a
            ! deck reaches the third pairing the backend supports: esp "none"
            ! with an mbe expansion, a plain many-body expansion through this
            ! module.
            if (trim(config%embedding) == "none") then
               expansion%esp = "none"
            end if

            ! TODO(mqc): refactor this ugly ass code, in general the expansion assignemtn
            expansion%far_field = config%fmo_far_field
            ! The deck's fragmentation level means the same thing here as it
            ! does for MBE: how many fragments at a time.
            expansion%level = max_level
            expansion%resppc = config%fmo_resppc
            expansion%max_outer = config%fmo_max_outer
            expansion%outer_tol = config%fmo_tolerance
            expansion%scf_max_iter = config%fmo_scf_max_iter
            expansion%scf_energy_tol = config%fmo_scf_energy_tol
            expansion%scf_density_tol = config%fmo_scf_density_tol
            ! From `keywords.scf`, the same source the unfragmented path reads.
            ! The three above stay on `keywords.fragmentation`, being
            ! per-fragment by intent.
            expansion%scf_drive%level_shift = config%method_config%scf%level_shift
            expansion%scf_drive%linear_dependence = config%method_config%scf%linear_dependence
            expansion%scf_drive%use_diis = config%method_config%scf%use_diis
            expansion%scf_drive%diis_size = config%method_config%scf%diis_size
            expansion%scf_drive%incremental_fock = config%method_config%scf%incremental_fock
            expansion%scf_drive%accelerator = config%method_config%scf%accelerator
            expansion%scf_drive%convergence_metric = config%method_config%scf%convergence_metric
            expansion%scf_drive%guess = config%method_config%scf%guess
            expansion%scf_drive%allow_crap_scf = config%method_config%scf%allow_crap_scf
            expansion%scf_drive%grad_tol = config%method_config%scf%gradient_convergence
            expansion%resources => resources
            expansion%node_leader_ranks = node_leader_ranks
            expansion%num_nodes = num_nodes
            expansion%global_groups = global_groups
            expansion%nodes_per_group = nodes_per_group
            expansion%group_leader_ranks = group_leader_ranks
            expansion%group_ids = group_ids
         end select
      else if (allow_overlapping_fragments) then
         ! GMBE: allocate gmbe_context_t
         allocate (gmbe_context_t :: expansion)
         select type (expansion)
         type is (gmbe_context_t)
            call expansion%init(config%method_config, config%calc_type)
            allocate (expansion%sys_geom, source=sys_geom)
            if (present(bonds)) then
               ! source=sys_geom above already copies its bonds when the caller
               ! embeds them there (the C API does), so replace rather than
               ! allocate -- allocating an allocated component aborts.
               if (allocated(expansion%sys_geom%bonds)) deallocate (expansion%sys_geom%bonds)
               allocate (expansion%sys_geom%bonds, source=bonds)
            end if
            ! TODO(mqc): `n_pie_terms` is set inside the rank-0 block above, so
            ! every other rank copies an undefined value here. Harmless only
            ! because the coordinator is the sole reader; `total_fragments`
            ! carries the same number and is the one that gets broadcast.
            expansion%n_pie_terms = n_pie_terms
            if (resources%mpi_comms%world_comm%rank() == 0) then
               allocate (expansion%pie_atom_sets, source=pie_atom_sets)
               allocate (expansion%pie_coefficients, source=pie_coefficients)
            end if
            expansion%resources => resources
            expansion%node_leader_ranks = node_leader_ranks
            expansion%num_nodes = num_nodes
            expansion%global_groups = global_groups
            expansion%nodes_per_group = nodes_per_group
            expansion%group_leader_ranks = group_leader_ranks
            expansion%group_ids = group_ids
         end select
      else
         ! Standard MBE: allocate mbe_context_t
         allocate (mbe_context_t :: expansion)
         select type (expansion)
         type is (mbe_context_t)
            call expansion%init(config%method_config, config%calc_type)
            allocate (expansion%sys_geom, source=sys_geom)
            if (present(bonds)) then
               ! source=sys_geom above already copies its bonds when the caller
               ! embeds them there (the C API does), so replace rather than
               ! allocate -- allocating an allocated component aborts.
               if (allocated(expansion%sys_geom%bonds)) deallocate (expansion%sys_geom%bonds)
               allocate (expansion%sys_geom%bonds, source=bonds)
            end if
            expansion%total_fragments = total_fragments
            if (resources%mpi_comms%world_comm%rank() == 0) then
               allocate (expansion%polymers, source=polymers)
            end if
            expansion%max_level = max_level
            expansion%resources => resources
            expansion%node_leader_ranks = node_leader_ranks
            expansion%num_nodes = num_nodes
            expansion%global_groups = global_groups
            expansion%nodes_per_group = nodes_per_group
            expansion%group_leader_ranks = group_leader_ranks
            expansion%group_ids = group_ids
         end select
      end if

      ! Opened on rank 0 only: it is the rank that collects every result, and N
      ! ranks appending to one file would interleave.
      !
      ! GMBE keys its terms by atom set, not by monomer tuple, so the file this
      ! writes would not describe them. Announced rather than opened and left
      ! unused: a checkpoint that silently records nothing is believed.
      if (len_trim(config%checkpoint_file) > 0 .and. allow_overlapping_fragments .and. &
          resources%mpi_comms%world_comm%rank() == 0) then
         call logger%warning("Checkpointing does not support GMBE; this run will write "// &
                             "nothing and resume nothing. Its PIE terms are atom sets "// &
                             "rather than monomer tuples, which the file cannot express.")
      end if

      if (resources%mpi_comms%world_comm%rank() == 0 .and. &
          .not. allow_overlapping_fragments .and. &
          len_trim(config%checkpoint_file) > 0) then
         call expansion%checkpoint%open(trim(config%checkpoint_file), &
                                        calculation_fingerprint(sys_geom, config%method_config, &
                                                                config%calc_type), &
                                        max_level + 1, &
                                        config%calc_type == CALC_TYPE_ENERGY, &
                                        checkpoint_error)
         if (checkpoint_error%has_error()) then
            ! A checkpoint from another calculation is not a warning: its
            ! energies would be spliced into this one and the total would come
            ! out converged and meaningless.
            call logger%error(checkpoint_error%get_message())
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end if
      end if

      ! Execute calculation using polymorphic dispatch
      if (resources%mpi_comms%world_comm%size() == 1) then
         call logger%large_info("Running in serial mode (single MPI rank)")
         call expansion%run_serial(json_data)
      else
         call expansion%run_distributed(json_data)
      end if
      call expansion%checkpoint%close()

      ! Clean up expansion context
      select type (expansion)
      type is (mbe_context_t)
         call expansion%destroy()
      type is (gmbe_context_t)
         call expansion%destroy()
      type is (fmo_context_t)
         call expansion%destroy()
      end select
      deallocate (expansion)

      ! Cleanup
      if (resources%mpi_comms%world_comm%rank() == 0) then
         if (allocated(polymers)) deallocate (polymers)
         if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
         if (allocated(group_leader_ranks)) deallocate (group_leader_ranks)
         if (allocated(group_ids)) deallocate (group_ids)
         if (allocated(pie_atom_sets)) deallocate (pie_atom_sets)
         if (allocated(pie_coefficients)) deallocate (pie_coefficients)
      end if

   end subroutine run_fragmented_calculation

   subroutine fragment_owner_map(sys_geom, owner, n_fragments)
      !! Which fragment each atom belongs to, from the declared monomers
      !!
      !! `fragment_atoms` is 0-indexed and padded to the widest fragment, so
      !! only the first `fragment_sizes(f)` entries of a column mean anything.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, allocatable, intent(out) :: owner(:)
      integer, intent(out) :: n_fragments

      integer :: f, k

      n_fragments = sys_geom%n_monomers
      allocate (owner(sys_geom%total_atoms), source=0)
      do f = 1, n_fragments
         do k = 1, sys_geom%fragment_sizes(f)
            owner(sys_geom%fragment_atoms(k, f) + 1) = f
         end do
      end do
   end subroutine fragment_owner_map

   subroutine run_multi_molecule_calculations(resources, mqc_config)
      !! Run calculations for multiple molecules with MPI parallelization
      !! Each molecule is independent, so assign one molecule per rank
      use mqc_config_types, only: mqc_config_t
      use mqc_config_adapter, only: config_to_system_geometry
      use mqc_error, only: error_t
      use mqc_io_helpers, only: set_molecule_suffix, get_output_json_filename
      use mqc_json, only: merge_multi_molecule_json

      type(resources_t), intent(in) :: resources
      type(mqc_config_t), intent(in) :: mqc_config

      type(driver_config_t) :: config
      type(system_geometry_t) :: sys_geom
      type(resources_t) :: mol_resources
      type(error_t) :: error
      integer :: imol, my_rank, num_ranks, color
      integer :: molecules_processed
      character(len=:), allocatable :: mol_name
      logical :: has_fragmented_molecules
      character(len=256), allocatable :: individual_json_files(:)

      my_rank = resources%mpi_comms%world_comm%rank()
      num_ranks = resources%mpi_comms%world_comm%size()

      ! Allocate array to track individual JSON files for merging
      allocate (individual_json_files(mqc_config%nmol))

      ! Check if any molecules have fragments (nlevel > 0)
      has_fragmented_molecules = .false.
      do imol = 1, mqc_config%nmol
         if (mqc_config%molecules(imol)%nfrag > 0) then
            has_fragmented_molecules = .true.
            exit
         end if
      end do

      if (my_rank == 0) then
         call logger%large_info(" ")
         call logger%large_info("============================================")
         call logger%large_info("Multi-molecule mode: "//to_char(mqc_config%nmol)//" molecules")
         call logger%large_info("MPI ranks: "//to_char(num_ranks))
         if (has_fragmented_molecules) then
            call logger%large_info("Mode: Sequential execution (fragmented molecules detected)")
            call logger%large_info("  Each molecule will use all "//to_char(num_ranks)//" rank(s) for its calculation")
         else if (num_ranks == 1) then
            call logger%large_info("Mode: Sequential execution (single rank)")
         else if (num_ranks > mqc_config%nmol) then
            call logger%large_info("Mode: Parallel execution (one molecule per rank)")
            call logger%large_info("Note: More ranks than molecules - ranks "//to_char(mqc_config%nmol)// &
                                   " to "//to_char(num_ranks - 1)//" will be idle")
         else
            call logger%large_info("Mode: Parallel execution (one molecule per rank)")
         end if
         call logger%large_info("============================================")
         call logger%large_info(" ")
      end if

      ! Determine execution mode:
      ! 1. Sequential: Single rank OR fragmented molecules (each molecule needs all ranks)
      ! 2. Parallel: Multiple ranks AND unfragmented molecules (distribute molecules across ranks)
      molecules_processed = 0

      if (num_ranks == 1 .or. has_fragmented_molecules) then
         ! Sequential mode: process all molecules one after another
         ! Each molecule uses all available ranks for its calculation
         do imol = 1, mqc_config%nmol
            ! Determine molecule name for logging
            if (allocated(mqc_config%molecules(imol)%name)) then
               mol_name = mqc_config%molecules(imol)%name
            else
               mol_name = "molecule_"//to_char(imol)
            end if

            if (my_rank == 0) then
               call logger%info(" ")
               call logger%info("--------------------------------------------")
               call logger%info("Processing molecule "//to_char(imol)//"/"//to_char(mqc_config%nmol)//": "//mol_name)
               call logger%info("--------------------------------------------")
            end if

            ! Convert to driver configuration for this molecule
            call config_to_driver(mqc_config, config, molecule_index=imol, &
                                  node_rank=resources%mpi_comms%node_comm%rank())

            ! Convert geometry for this molecule
            call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
            if (error%has_error()) then
               call error%add_context("mqc_driver:run_multi_molecule_calculation")
               if (my_rank == 0) then
                  call logger%error("Error converting geometry for "//mol_name//": "//error%get_full_trace())
               end if
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            ! Set output filename suffix for this molecule
            call set_molecule_suffix("_"//trim(mol_name))

            ! Run calculation for this molecule
            call run_calculation(resources, config, sys_geom, mqc_config%molecules(imol)%bonds)

            ! Track the JSON filename for later merging
            individual_json_files(imol) = get_output_json_filename()

            ! Clean up for this molecule
            call sys_geom%destroy()

            if (my_rank == 0) then
               call logger%info("Completed molecule "//to_char(imol)//"/"//to_char(mqc_config%nmol)//": "//mol_name)
            end if
            molecules_processed = molecules_processed + 1
         end do
      else
         ! Multiple ranks: distribute molecules across ranks in round-robin fashion
         molecules_processed = 0
         do imol = 1, mqc_config%nmol
            ! This rank processes molecules where (imol - 1) mod num_ranks == my_rank
            if (mod(imol - 1, num_ranks) == my_rank) then
               ! Determine molecule name for logging
               if (allocated(mqc_config%molecules(imol)%name)) then
                  mol_name = mqc_config%molecules(imol)%name
               else
                  mol_name = "molecule_"//to_char(imol)
               end if

               call logger%info(" ")
               call logger%info("--------------------------------------------")
               call logger%info("Rank "//to_char(my_rank)//": Processing molecule "//to_char(imol)// &
                                "/"//to_char(mqc_config%nmol)//": "//mol_name)
               call logger%info("--------------------------------------------")

               ! Convert to driver configuration for this molecule
               call config_to_driver(mqc_config, config, molecule_index=imol, &
                                     node_rank=resources%mpi_comms%node_comm%rank())

               ! Convert geometry for this molecule
               call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
               if (error%has_error()) then
                  call error%add_context("mqc_driver:run_multi_molecule_calculation")
  call logger%error("Rank "//to_char(my_rank)//": Error converting geometry for "//mol_name//": "//error%get_full_trace())
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               ! Set output filename suffix for this molecule
               call set_molecule_suffix("_"//trim(mol_name))

               ! Run calculation for this molecule (all ranks write JSON in parallel mode)
               call run_calculation(resources, config, sys_geom, mqc_config%molecules(imol)%bonds, &
                                    all_ranks_write_json=.true.)

               ! Track the JSON filename for later merging
               individual_json_files(imol) = get_output_json_filename()

               ! Clean up for this molecule
               call sys_geom%destroy()

               call logger%info("Rank "//to_char(my_rank)//": Completed molecule "//to_char(imol)// &
                                "/"//to_char(mqc_config%nmol)//": "//mol_name)

               molecules_processed = molecules_processed + 1
            end if
         end do

         if (molecules_processed == 0) then
            ! Idle rank - no molecules assigned
            call logger%verbose("Rank "//to_char(my_rank)//": No molecules assigned (idle)")
         end if
      end if

      ! Synchronize all ranks
      call resources%mpi_comms%world_comm%barrier()

      ! In parallel execution, rank 0 needs to reconstruct all JSON filenames for merging
      ! since each rank only populated its own entry
      if (my_rank == 0 .and. num_ranks > 1 .and. .not. has_fragmented_molecules) then
         ! Rank 0 constructs filenames for all molecules
         do imol = 1, mqc_config%nmol
            ! Get molecule name
            if (allocated(mqc_config%molecules(imol)%name)) then
               mol_name = mqc_config%molecules(imol)%name
            else
               mol_name = "molecule_"//to_char(imol)
            end if
            ! Construct JSON filename pattern: output_<basename>_<molname>.json
            ! This mirrors what get_output_json_filename() returns after set_molecule_suffix()
            call set_molecule_suffix("_"//trim(mol_name))
            individual_json_files(imol) = get_output_json_filename()
         end do
      end if

      ! Merge individual JSON files into one combined file (rank 0 only)
      if (my_rank == 0) then
         call merge_multi_molecule_json(individual_json_files, mqc_config%nmol)
      end if

      if (my_rank == 0) then
         call logger%large_info(" ")
         call logger%large_info("============================================")
         call logger%large_info("All "//to_char(mqc_config%nmol)//" molecules completed")
         if (has_fragmented_molecules) then
            call logger%large_info("Execution: Sequential (each molecule used all ranks)")
         else if (num_ranks == 1) then
            call logger%large_info("Execution: Sequential (single rank)")
         else if (num_ranks > mqc_config%nmol) then
     call logger%large_info("Execution: Parallel (active ranks: "//to_char(mqc_config%nmol)//"/"//to_char(num_ranks)//")")
         else
            call logger%large_info("Execution: Parallel (all ranks active)")
         end if
         call logger%large_info("============================================")
      end if

   end subroutine run_multi_molecule_calculations

   subroutine run_sapt(config, sys_geom, rank, write_output, result_out)
      !! SAPT0 or SAPT2 between the deck's two fragments
      !!
      !! The monomers are the deck's own `fragments`, so no new keyword is
      !! needed. Rank zero only.
      !!
      !! **Exactly two fragments.** SAPT partitions the Hamiltonian as
      !! `H_A + H_B + V`; there is no slot for a third monomer, so a cluster is
      !! one SAPT calculation per pair. See `validation/check_sapt.f90`, which
      !! walks the pairs of a six-water prism.
      use mqc_czt_bridge, only: run_czt_sapt0, run_czt_sapt2
      use mqc_method_types, only: METHOD_TYPE_SAPT2
      use mqc_program_limits, only: N_SAPT_TERMS, N_SAPT2_TERMS
      use mqc_elements, only: element_number_to_symbol
      use mqc_json_output_types, only: OUTPUT_MODE_UNFRAGMENTED
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out), optional :: result_out

      type(error_t) :: err
      type(json_output_data_t) :: json_data
      real(dp), allocatable :: terms(:)
      logical :: is_sapt2
      integer, allocatable :: z_a(:), z_b(:)
      real(dp), allocatable :: xyz_a(:, :), xyz_b(:, :)
      character(len=8), allocatable :: sym_a(:), sym_b(:)
      integer :: na, nb, i, a
      integer :: charge_a, charge_b, mult_a, mult_b

      if (rank /= 0) return

      if (sys_geom%n_monomers /= 2) then
         call refuse(result_out, "SAPT: the system must have exactly two "// &
                     "fragments; this one has "//to_char(sys_geom%n_monomers)// &
                     ". SAPT is a two-body theory, so a cluster is one "// &
                     "calculation per pair.")
         return
      end if

      ! The monomers' own charges and multiplicities, which decide how many
      ! electrons each SCF is run with. The arrays are optional in
      ! `system_geometry_t`; a deck that omits them means neutral singlets.
      charge_a = 0
      charge_b = 0
      mult_a = 1
      mult_b = 1
      if (allocated(sys_geom%fragment_charges)) then
         charge_a = sys_geom%fragment_charges(1)
         charge_b = sys_geom%fragment_charges(2)
      end if
      if (allocated(sys_geom%fragment_multiplicities)) then
         mult_a = sys_geom%fragment_multiplicities(1)
         mult_b = sys_geom%fragment_multiplicities(2)
      end if

      ! Both monomer references are RHF. Refused here rather than at the
      ! electron-count parity check further in, which catches only the
      ! odd-electron half: a triplet has an even count and would otherwise be
      ! solved as the singlet.
      if (mult_a /= 1 .or. mult_b /= 1) then
         call refuse(result_out, "SAPT: both monomers must be closed shell; "// &
                     "this deck asks for multiplicities "//to_char(mult_a)// &
                     " and "//to_char(mult_b)//". The monomer references here "// &
                     "are RHF, so an open-shell monomer has no wavefunction to "// &
                     "expand about.")
         return
      end if

      na = sys_geom%fragment_sizes(1)
      nb = sys_geom%fragment_sizes(2)
      allocate (z_a(na), xyz_a(3, na), sym_a(na))
      allocate (z_b(nb), xyz_b(3, nb), sym_b(nb))
      do i = 1, na
         a = sys_geom%fragment_atoms(i, 1) + 1        ! stored 0-based
         z_a(i) = sys_geom%element_numbers(a)
         xyz_a(:, i) = sys_geom%coordinates(:, a)
         sym_a(i) = element_number_to_symbol(z_a(i))
      end do
      do i = 1, nb
         a = sys_geom%fragment_atoms(i, 2) + 1
         z_b(i) = sys_geom%element_numbers(a)
         xyz_b(:, i) = sys_geom%coordinates(:, a)
         sym_b(i) = element_number_to_symbol(z_b(i))
      end do

      is_sapt2 = config%method_config%method_type == METHOD_TYPE_SAPT2
      if (is_sapt2) then
         allocate (terms(N_SAPT2_TERMS))
         call run_czt_sapt2(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                            config%method_config%basis_set, &
                            charge_a, charge_b, terms, err)
      else
         allocate (terms(N_SAPT_TERMS))
         call run_czt_sapt0(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                            config%method_config%basis_set, &
                            charge_a, charge_b, terms, err)
      end if
      if (err%has_error()) then
         call refuse(result_out, "SAPT: "//err%get_message())
         return
      end if

      ! The response terms are what enters the total; the uncoupled ones and the
      ! S^2 exchange are printed beside them because that is what another code's
      ! output is usually quoted in.
      call logger%info("============================================")
      if (is_sapt2) then
         call logger%info("  SAPT2 interaction energy, Hartree")
      else
         call logger%info("  SAPT0 interaction energy, Hartree")
      end if
      call logger%info("    electrostatics        "//to_char(terms(1)))
      call logger%info("    exchange              "//to_char(terms(3)))
      call logger%info("      (S^2 approximation) "//to_char(terms(2)))
      call logger%info("    induction             "//to_char(terms(5)))
      call logger%info("      (uncoupled)         "//to_char(terms(4)))
      call logger%info("    exchange-induction    "//to_char(terms(7)))
      call logger%info("      (uncoupled)         "//to_char(terms(6)))
      call logger%info("    dispersion            "//to_char(terms(8)))
      call logger%info("    exchange-dispersion   "//to_char(terms(9)))
      call logger%info("    delta HF              "//to_char(terms(10)))
      call logger%info("    --")
      call logger%info("    supermolecular HF     "//to_char(terms(11)))
      if (is_sapt2) then
         call logger%info("    total SAPT0           "//to_char(terms(12)))
         call logger%info("    --")
         call logger%info("    elst12                "//to_char(terms(13)))
         call logger%info("    exch11                "//to_char(terms(14)))
         call logger%info("    exch12                "//to_char(terms(15)))
         call logger%info("    ind22                 "//to_char(terms(16)))
         call logger%info("    exch-ind22 (scaled)   "//to_char(terms(17)))
         call logger%info("    --")
         call logger%info("    total                 "//to_char(terms(18)))
      else
         call logger%info("    total                 "//to_char(terms(12)))
      end if
      call logger%info("============================================")

      if (present(result_out)) then
         ! The `scf` slot, being the reference `energy_t%total()` sums. This is
         ! an interaction energy and not an SCF energy.
         result_out%energy%scf = terms(size(terms))
         result_out%has_energy = .true.
      end if

      if (write_output .and. .not. config%skip_json_output) then
         json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
         json_data%total_energy = terms(size(terms))
         json_data%has_energy = .true.
         json_data%sapt_terms = terms
         json_data%has_sapt = .true.
         json_data%fragment_breakdown = config%fragment_breakdown
         call write_json_output(json_data)
         call json_data%destroy()
      end if

      deallocate (z_a, xyz_a, sym_a, z_b, xyz_b, sym_b)
   end subroutine run_sapt

   subroutine run_efp(config, sys_geom, rank, write_output, result_out)
      !! The interaction energy of a set of effective fragments
      !!
      !! Each fragment of the deck names a potential; the backend loads them,
      !! places each on the atoms the deck gave it, turns it into that
      !! orientation and evaluates the five EFP2 terms. There is no SCF here:
      !! the wavefunctions were solved when the potentials were made.
      !!
      !! Rank zero only. The work is in `run_czt_efp`, which lives behind
      !! `MQC_ENABLE_CZT`; the stub declines with the same signature and
      !! names the build option.
      use mqc_czt_bridge, only: run_czt_efp
      use mqc_program_limits, only: N_EFP_TERMS
      use mqc_json_output_types, only: OUTPUT_MODE_UNFRAGMENTED
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out), optional :: result_out

      type(error_t) :: err
      type(json_output_data_t) :: json_data
      real(dp) :: terms(N_EFP_TERMS)
      integer :: n

      if (rank /= 0) return

      n = sys_geom%n_monomers
      if (.not. allocated(config%fragment_potentials)) then
         call refuse(result_out, "EFP: no fragment carries a potential -- name "// &
                     "one per fragment in 'fragment_potentials' in the deck, or "// &
                     "through mqc_system_set_fragment_potentials from the C and "// &
                     "Python interfaces")
         return
      end if
      if (size(config%fragment_potentials) /= n) then
         call refuse(result_out, "EFP: the system has "//to_char(n)//" fragments "// &
                     "but "//to_char(size(config%fragment_potentials))//" potentials")
         return
      end if

      call run_czt_efp(config%fragment_potentials, sys_geom%fragment_sizes, &
                       sys_geom%fragment_atoms, sys_geom%coordinates, terms, err)
      if (err%has_error()) then
         call refuse(result_out, "EFP: "//err%get_message())
         return
      end if

      call logger%info("============================================")
      call logger%info("  EFP2 interaction energy, Hartree")
      call logger%info("    electrostatics      "//to_char(terms(1)))
      call logger%info("    polarization        "//to_char(terms(2)))
      call logger%info("    exchange repulsion  "//to_char(terms(3)))
      call logger%info("    dispersion          "//to_char(terms(4)))
      call logger%info("    charge transfer     "//to_char(terms(5)))
      call logger%info("    total               "//to_char(terms(6)))
      call logger%info("============================================")

      if (present(result_out)) then
         ! In the `scf` slot because that is the one `energy_t%total()` sums as
         ! the reference. It is not an SCF energy: no wavefunction is solved
         ! here at all.
         result_out%energy%scf = terms(6)
         result_out%has_energy = .true.
      end if

      ! `UNFRAGMENTED` because that is the shape of what is written -- one
      ! energy for one system -- not a claim that the system has no fragments.
      if (write_output .and. .not. config%skip_json_output) then
         json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
         json_data%total_energy = terms(6)
         json_data%has_energy = .true.
         json_data%fragment_breakdown = config%fragment_breakdown
         call write_json_output(json_data)
         call json_data%destroy()
      end if
   end subroutine run_efp

   subroutine refuse(result_out, message)
      !! Log a refusal and, when there is a caller, put it on the result
      !!
      !! Marking the result matters on the interaction-energy paths: `mqc_run`
      !! reads the energy off a result nothing marked as failed, and a zero
      !! interaction energy is a physically plausible number.
      type(calculation_result_t), intent(inout), optional :: result_out
      character(len=*), intent(in) :: message

      call logger%error(message)
      if (present(result_out)) then
         call result_out%error%set(ERROR_VALIDATION, message)
         result_out%has_error = .true.
      end if
   end subroutine refuse

   subroutine run_efmo_energy(config, sys_geom, rank, write_output, result_out)
      !! An effective fragment molecular orbital energy for a partitioned system
      !!
      !! Rank zero only, as MAKEFP and NEO are: the work is one MAKEFP per
      !! fragment and one SCF per near dimer, all threaded inside the integral
      !! backend, and distributing it over ranks is Phase 4.
      !!
      !! Neither the fragmented nor the unfragmented path applies. EFMO is a
      !! fragmented method with no n-mer list: every fragment's potential is
      !! built, the pairs are split at `R_cut`, the near ones get an in-vacuo
      !! dimer and the far ones four effective-fragment terms, and one induction
      !! runs over all of them. `run_czt_efmo` assembles it and hands back the
      !! breakdown.
      use mqc_czt_bridge, only: run_czt_efmo
      use mqc_method_types, only: METHOD_TYPE_HF
      use mqc_elements, only: element_number_to_symbol
      use mqc_program_limits, only: N_EFMO_TERMS
      use mqc_json_output_types, only: OUTPUT_MODE_UNFRAGMENTED
      use pic_logger, only: verbose_level
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out), optional :: result_out

      type(error_t) :: err
      type(json_output_data_t) :: json_data
      character(len=8), allocatable :: symbols(:)
      integer, allocatable :: owner(:), charges(:)
      real(dp) :: terms(N_EFMO_TERMS)
      real(dp) :: energy
      type(scf_numerics_t) :: efmo_scf
      integer :: i, n_frag, n_qm, n_efp

      integer, parameter :: EFMO_SCF_MAX_ITER = 200
      real(dp), parameter :: EFMO_SCF_ENERGY_TOL = 1.0e-10_dp
      real(dp), parameter :: EFMO_SCF_DENSITY_TOL = 1.0e-8_dp
      real(dp), parameter :: EFMO_SCF_GRAD_TOL = 1.0e-8_dp
         !! What every SCF in an EFMO run is converged to, monomer and dimer
         !! alike, and `make_efp_potential`'s own defaults.
         !!
         !! Not read from `keywords.scf`, and deliberately tighter than a
         !! whole-system run's: the near-dimer correction is
         !! `E_IJ^0 - E_I^0 - E_J^0`, four orders smaller than any of the three,
         !! so a monomer and its dimer converged to 1e-6 leave a correction with
         !! no significant figures. There is no deck key for these yet, since a
         !! looser EFMO is not a cheaper EFMO -- the cost is MAKEFP.

      if (rank /= 0) return

      ! Restricted Hartree-Fock only in this phase. A correlated `E_I^0` runs on
      ! the same orbitals afterwards and is Phase 4; DFT would need a MAKEFP
      ! that is not restricted to a Hartree-Fock reference. Refused by name
      ! rather than silently run as HF.
      if (config%method_config%method_type /= METHOD_TYPE_HF) then
         call refuse(result_out, "EFMO runs restricted Hartree-Fock fragments for "// &
                     "now, and model.method is '"// &
                     trim(method_type_to_string(config%method_config%method_type))// &
                     "'. Correlated fragments and Kohn-Sham fragments are not "// &
                     "implemented yet.")
         return
      end if
      if (config%method_config%scf%unrestricted) then
         call refuse(result_out, "EFMO is closed-shell for now: every fragment and "// &
                     "every dimer is solved with restricted Hartree-Fock, so "// &
                     "keywords.scf.unrestricted cannot be honoured.")
         return
      end if

      call fragment_owner_map(sys_geom, owner, n_frag)
      if (n_frag < 2) then
         call refuse(result_out, "EFMO needs at least two fragments: with one there "// &
                     "is no pair to split and the energy is an ordinary SCF.")
         return
      end if
      if (any(owner == 0)) then
         call refuse(result_out, "EFMO: some atom belongs to no fragment. Every atom "// &
                     "has to be in exactly one, since each fragment's potential is "// &
                     "built from its own atoms alone.")
         return
      end if

      allocate (charges(n_frag), source=0)
      if (allocated(sys_geom%fragment_charges)) charges = sys_geom%fragment_charges

      allocate (symbols(sys_geom%total_atoms))
      do i = 1, sys_geom%total_atoms
         symbols(i) = element_number_to_symbol(sys_geom%element_numbers(i))
      end do

      ! `keywords.scf` drives every SCF here, monomer and dimer alike, the same
      ! way it drives a MakeFP run. The tolerances below are deliberately not
      ! taken from it: a fragment potential and the dimer whose difference
      ! against it is an interaction energy are converged tighter than a
      ! whole-system run would be, and `make_efp_potential` already holds those
      ! defaults.
      efmo_scf%level_shift = config%method_config%scf%level_shift
      efmo_scf%linear_dependence = config%method_config%scf%linear_dependence
      efmo_scf%use_diis = config%method_config%scf%use_diis
      efmo_scf%diis_size = config%method_config%scf%diis_size
      efmo_scf%incremental_fock = config%method_config%scf%incremental_fock
      efmo_scf%accelerator = config%method_config%scf%accelerator
      efmo_scf%convergence_metric = config%method_config%scf%convergence_metric
      efmo_scf%allow_crap_scf = config%method_config%scf%allow_crap_scf

      call logger%info("Running EFMO over "//to_char(n_frag)//" fragments")

      call run_czt_efmo(sys_geom%element_numbers, symbols, sys_geom%coordinates, owner, &
                        charges, config%method_config%basis_set, &
                        config%method_config%efmo%rcut, &
                        config%method_config%efmo%charge_transfer, efmo_scf, &
                        EFMO_SCF_MAX_ITER, EFMO_SCF_ENERGY_TOL, EFMO_SCF_DENSITY_TOL, &
                        EFMO_SCF_GRAD_TOL, trim(config%method_config%scf%guess), &
                        energy, terms, n_qm, n_efp, err, &
                        verbose=(logger%log_level >= verbose_level), &
                        vdwscl=config%method_config%efp%vdw_scale, &
                        quadrupole_blocks=config%method_config%efp%quadrupole_blocks, &
                        dynamic_tol=config%method_config%efp%dynamic_tolerance, &
                        dynamic_maxiter=config%method_config%efp%dynamic_maxiter, &
                        response=config%method_config%efp%response, &
                        allow_crap_response=config%method_config%efp%allow_crap_response, &
                        response_batch=config%method_config%efp%response_batch)
      if (err%has_error()) then
         call refuse(result_out, "EFMO: "//err%get_message())
         return
      end if

      if (present(result_out)) then
         ! In the `scf` slot, which is what `energy_t%total()` reports. It is
         ! not one SCF energy -- it is a sum over fragments, dimers and pair
         ! terms -- but an expansion total has no correlation breakdown to put
         ! anywhere else, which is what the MBE path does with its own.
         result_out%energy%scf = energy
         result_out%has_energy = .true.
      end if

      ! `UNFRAGMENTED` because that is the shape of what is written -- one
      ! energy for one system, with a named breakdown beside it -- not a claim
      ! that the system has no fragments.
      if (write_output .and. .not. config%skip_json_output) then
         json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
         json_data%total_energy = energy
         json_data%has_energy = .true.
         json_data%efmo_terms = terms
         json_data%efmo_qm_dimers = n_qm
         json_data%efmo_efp_dimers = n_efp
         json_data%has_efmo = .true.
         json_data%fragment_breakdown = config%fragment_breakdown
         call write_json_output(json_data)
         call json_data%destroy()
      end if
   end subroutine run_efmo_energy

   subroutine run_makefp(config, sys_geom, rank, result_out)
      !! Build an effective fragment potential for the whole system and write it
      !!
      !! Rank zero only: the work is one SCF and a set of response solves, all
      !! threaded inside the integral backend, and the write is a single file.
      use mqc_calc_types, only: CALC_TYPE_MAKEFP
      use mqc_elements, only: element_number_to_symbol
      use mqc_czt_bridge, only: run_czt_makefp
      use mqc_io_helpers, only: get_basename
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      type(calculation_result_t), intent(out), optional :: result_out
         !! Carries the failure, if there is one. There is no energy to put on
         !! it -- MAKEFP produces a file -- but a caller that gets neither an
         !! energy nor an error has no way to tell a written potential from a
         !! backend that declined.

      type(error_t) :: err
      character(len=8), allocatable :: symbols(:)
      character(len=:), allocatable :: path, name
      integer :: i
      real(dp), allocatable :: named_energy_tol, named_density_tol, named_grad_tol
      integer, allocatable :: named_max_iter
      type(scf_numerics_t) :: makefp_scf  !! SCF settings the MAKEFP run uses

      if (rank /= 0) return

      allocate (symbols(sys_geom%total_atoms))
      do i = 1, sys_geom%total_atoms
         symbols(i) = element_number_to_symbol(sys_geom%element_numbers(i))
      end do

      ! Named off the deck, the way the JSON output is: a run on water.json
      ! leaves water.efp beside it.
      name = trim(get_basename())
      path = trim(name)//".efp"

      call logger%info("Building an effective fragment potential")

      if (config%method_config%scf%energy_convergence_set) then
         named_energy_tol = config%method_config%scf%energy_convergence
      end if
      if (config%method_config%scf%density_convergence_set) then
         named_density_tol = config%method_config%scf%density_convergence
      end if

      if (config%method_config%scf%gradient_convergence > 0.0_dp) then
         named_grad_tol = config%method_config%scf%gradient_convergence
      end if

      if (config%method_config%scf%max_iter_set) then
         named_max_iter = config%method_config%scf%max_iter
      end if
      ! Everything else about how the SCF runs goes down whole. These are the
      ! settings a MakeFP deck could set and then watch do nothing.
      makefp_scf%level_shift = config%method_config%scf%level_shift
      makefp_scf%linear_dependence = config%method_config%scf%linear_dependence
      makefp_scf%use_diis = config%method_config%scf%use_diis
      makefp_scf%diis_size = config%method_config%scf%diis_size
      makefp_scf%incremental_fock = config%method_config%scf%incremental_fock
      makefp_scf%accelerator = config%method_config%scf%accelerator
      if (config%method_config%scf%density_fitting) then
         ! TODO(mqc): make good
         call run_czt_makefp(sys_geom%element_numbers, symbols, sys_geom%coordinates, &
                             config%method_config%basis_set, name, path, err, &
                             charge=sys_geom%charge, verbose=.true., &
                             aux_basis=trim(config%method_config%scf%aux_basis_set), &
                             guess=trim(config%method_config%scf%guess), &
                             energy_tol=named_energy_tol, &
                             density_tol=named_density_tol, &
                             grad_tol=named_grad_tol, &
                             scf_in=makefp_scf, max_iter_in=named_max_iter, &
                             vdwscl=config%method_config%efp%vdw_scale, &
                             quadrupole_blocks=config%method_config%efp%quadrupole_blocks, &
                             dynamic_tol=config%method_config%efp%dynamic_tolerance, &
                             dynamic_maxiter=config%method_config%efp%dynamic_maxiter, &
                             response=config%method_config%efp%response, &
                             allow_crap_response=config%method_config%efp%allow_crap_response, &
                             response_batch=config%method_config%efp%response_batch)
      else
         call run_czt_makefp(sys_geom%element_numbers, symbols, sys_geom%coordinates, &
                             config%method_config%basis_set, name, path, err, &
                             charge=sys_geom%charge, verbose=.true., &
                             guess=trim(config%method_config%scf%guess), &
                             energy_tol=named_energy_tol, &
                             density_tol=named_density_tol, &
                             grad_tol=named_grad_tol, &
                             scf_in=makefp_scf, max_iter_in=named_max_iter, &
                             vdwscl=config%method_config%efp%vdw_scale, &
                             quadrupole_blocks=config%method_config%efp%quadrupole_blocks, &
                             dynamic_tol=config%method_config%efp%dynamic_tolerance, &
                             dynamic_maxiter=config%method_config%efp%dynamic_maxiter, &
                             response=config%method_config%efp%response, &
                             allow_crap_response=config%method_config%efp%allow_crap_response, &
                             response_batch=config%method_config%efp%response_batch)
      end if
      if (err%has_error()) then
         call refuse(result_out, "MAKEFP failed: "//err%get_message())
         return
      end if
      call logger%info("Wrote "//trim(path))
   end subroutine run_makefp

   subroutine run_neo(config, sys_geom, rank, write_output, result_out)
      !! A nuclear-electronic orbital energy for the whole system
      !!
      !! Rank zero only, like MAKEFP: the work is a few coupled SCFs, all
      !! threaded inside the integral backend. Neither the fragmented nor the
      !! unfragmented path applies -- a quantised proton is solved with every
      !! electron of the system, and there is one energy at the end.
      use mqc_czt_bridge, only: run_czt_neo
      use mqc_method_types, only: METHOD_TYPE_HF, METHOD_TYPE_DFT
      use mqc_elements, only: element_number_to_symbol
      use mqc_json_output_types, only: OUTPUT_MODE_UNFRAGMENTED
      use mqc_program_limits, only: MAX_LINE_LENGTH
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: rank
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out), optional :: result_out

      type(error_t) :: err
      type(json_output_data_t) :: json_data
      character(len=8), allocatable :: symbols(:)
      logical, allocatable :: quantum(:)
      real(dp), allocatable :: named_energy_tol, named_density_tol
      integer, allocatable :: named_max_iter
      type(scf_numerics_t) :: neo_scf  !! How every electronic SCF in the run is driven
      real(dp) :: energy
      integer :: i, k
      character(len=MAX_LINE_LENGTH) :: line
      character(len=:), allocatable :: functional

      if (rank /= 0) return
      select case (config%method_config%method_type)
      case (METHOD_TYPE_HF)
         functional = ""
         if (len_trim(config%method_config%neo%epc) > 0) then
            call refuse(result_out, "keywords.neo.epc is an electron-proton correlation "// &
                        "functional and needs model.method dft")
            return
         end if
      case (METHOD_TYPE_DFT)
         functional = trim(config%method_config%dft%functional)
      case default
         call refuse(result_out, "keywords.neo: quantum nuclei are implemented for "// &
                     "model.method hf and dft only so far")
         return
      end select

      allocate (symbols(sys_geom%total_atoms), quantum(sys_geom%total_atoms))
      do i = 1, sys_geom%total_atoms
         symbols(i) = element_number_to_symbol(sys_geom%element_numbers(i))
      end do
      quantum = .false.
      if (allocated(config%method_config%neo%quantum_indices)) then
         do k = 1, size(config%method_config%neo%quantum_indices)
            i = config%method_config%neo%quantum_indices(k)
            if (i < 1 .or. i > sys_geom%total_atoms) then
               write (line, "(A,I0,A,I0,A)") "keywords.neo.quantum_nuclei: atom index ", &
                  i - 1, " is outside this system's ", sys_geom%total_atoms, &
                  " atoms (indices are 0-based)"
               call refuse(result_out, trim(line))
               return
            end if
            quantum(i) = .true.
         end do
      end if
      if (allocated(config%method_config%neo%quantum_symbols)) then
         do k = 1, size(config%method_config%neo%quantum_symbols)
            do i = 1, sys_geom%total_atoms
               if (same_symbol(symbols(i), config%method_config%neo%quantum_symbols(k))) then
                  quantum(i) = .true.
               end if
            end do
         end do
      end if
      if (count(quantum) == 0) then
         call refuse(result_out, "keywords.neo.quantum_nuclei names no atom of this system")
         return
      end if

      write (line, "(A,I0,A,A)") "Nuclear-electronic orbital "// &
         trim(merge("DFT         ", "Hartree-Fock", len(functional) > 0))//": ", count(quantum), &
         " quantum nucleus/nuclei in the ", trim(config%method_config%neo%nuclear_basis)
      call logger%info(trim(line)//" basis")
      if (len(functional) > 0) then
         line = "  functional "//functional
         if (len_trim(config%method_config%neo%epc) > 0) then
            line = trim(line)//", electron-proton correlation epc"//trim(config%method_config%neo%epc)
         else
            line = trim(line)//", no electron-proton correlation functional"
         end if
         call logger%info(trim(line))
      end if

      if (config%method_config%scf%energy_convergence_set) then
         named_energy_tol = config%method_config%scf%energy_convergence
      end if
      if (config%method_config%scf%density_convergence_set) then
         named_density_tol = config%method_config%scf%density_convergence
      end if
      if (config%method_config%scf%max_iter_set) then
         named_max_iter = config%method_config%scf%max_iter
      end if
      ! Everything else about how the SCF runs goes down whole, as MAKEFP does
      ! it: a deck could set these and otherwise watch them do nothing.
      neo_scf%level_shift = config%method_config%scf%level_shift
      neo_scf%linear_dependence = config%method_config%scf%linear_dependence
      neo_scf%use_diis = config%method_config%scf%use_diis
      neo_scf%diis_size = config%method_config%scf%diis_size
      neo_scf%incremental_fock = config%method_config%scf%incremental_fock
      neo_scf%accelerator = config%method_config%scf%accelerator
      call run_czt_neo(sys_geom%element_numbers, symbols, sys_geom%coordinates, &
                       config%method_config%basis_set, &
                       trim(config%method_config%neo%nuclear_basis), quantum, &
                       sys_geom%charge, energy, err, verbose=.true., &
                       energy_tol=named_energy_tol, density_tol=named_density_tol, &
                       max_iter=named_max_iter, functional=functional, &
                       grid_level=config%method_config%dft%grid_level, &
                       epc=trim(config%method_config%neo%epc), scf_in=neo_scf)
      if (err%has_error()) then
         call refuse(result_out, "NEO failed: "//err%get_message())
         return
      end if

      write (line, "(A,F20.10)") "NEO"//trim(merge("-DFT", "-HF ", len(functional) > 0))// &
         " total energy: ", energy
      call logger%info(trim(line))
      if (present(result_out)) then
         result_out%energy%scf = energy
         result_out%has_energy = .true.
      end if
      if (write_output .and. .not. config%skip_json_output) then
         json_data%output_mode = OUTPUT_MODE_UNFRAGMENTED
         json_data%total_energy = energy
         json_data%has_energy = .true.
         json_data%fragment_breakdown = config%fragment_breakdown
         call write_json_output(json_data)
         call json_data%destroy()
      end if
   end subroutine run_neo

   pure function same_symbol(a, b) result(same)
      !! Element symbols compared without regard to case or padding
      character(len=*), intent(in) :: a, b
      logical :: same
      character(len=len(a)) :: la
      character(len=len(b)) :: lb
      integer :: i
      la = a
      lb = b
      do i = 1, len(la)
         if (la(i:i) >= "A" .and. la(i:i) <= "Z") la(i:i) = achar(iachar(la(i:i)) + 32)
      end do
      do i = 1, len(lb)
         if (lb(i:i) >= "A" .and. lb(i:i) <= "Z") lb(i:i) = achar(iachar(lb(i:i)) + 32)
      end do
      same = trim(adjustl(la)) == trim(adjustl(lb))
   end function same_symbol

end module mqc_driver
