!! Adapter module to convert mqc_config_t to internal driver structures
module mqc_config_adapter
   !! Converts a parsed input deck into the structures the driver and the
   !! methods read: `driver_config_t` and `system_geometry_t`.
   use pic_types, only: dp, int32
   use mqc_config_types, only: mqc_config_t
   use mqc_physical_fragment, only: system_geometry_t, to_bohr
   use mqc_elements, only: element_symbol_to_number
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_string_utils, only: int_to_text
   use mqc_calculation_keywords, only: hessian_keywords_t, aimd_keywords_t, scf_keywords_t
   use mqc_optimizer_types, only: optimizer_settings_t, &
                                  coordinates_from_string, algorithm_from_string, &
                                  hessian_update_from_string, target_from_string, &
                                  OPT_COORDS_UNKNOWN, OPT_ALGO_UNKNOWN, &
                                  OPT_HESSIAN_UPDATE_ENGINE, OPT_TARGET_MINIMUM, &
                                  OPT_TARGET_SADDLE, OPT_COORDS_CARTESIAN, &
                                  algorithm_to_string, algorithm_finds_saddle, &
                                  neb_ends_from_string, NEB_ENDS_FROZEN, &
                                  saddle_method_from_string, SADDLE_METHOD_PRFO, &
                                  MIN_NEB_IMAGES
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_CCSD_T, METHOD_TYPE_GFN1, &
                               METHOD_TYPE_GFN2, METHOD_TYPE_EFP2
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: driver_config_t  !! Minimal config for driver
   public :: config_to_driver, config_to_system_geometry
   public :: get_logger_level  !! Convert log level string to integer
   public :: check_fragment_overlap  !! Check for overlapping fragments (for testing)
   public :: check_counterpoise_support  !! Refuse a counterpoise this expansion cannot honour

   type :: driver_config_t
      !! Runtime configuration for the driver (internal use only)
      ! Core calculation settings
      integer(int32) :: calc_type   !! Calculation type constant

      ! Method configuration (includes XTB solvation, DFT settings, etc.)
      type(method_config_t) :: method_config  !! Complete method configuration
      character(len=:), allocatable :: checkpoint_file  !! Append/resume path, empty for neither

      ! Fragmentation settings
      integer :: nlevel = 0         !! Fragmentation level (0 = unfragmented)
      logical :: allow_overlapping_fragments = .false.  !! Enable GMBE for overlapping fragments
      character(len=16) :: expansion_kind = "mbe"  !! "mbe", "fmo" or "ee-mbe"
      character(len=16) :: embedding = ""
         !! What field the fragments sit in: "none" turns the embedding off and
         !! leaves a plain many-body expansion. Empty means the deck said nothing
         !! and the expansion picks its own default.
      character(len=16) :: bond_breaking = "none"
         !! How a cut covalent bond is represented; "none" refuses one
      real(dp) :: cap_scale = 1.0_dp
         !! Where a cap sits along the bond it closes
      character(len=16) :: counterpoise = "none"   !! "none" or "vmfc"
      character(len=16) :: fmo_far_field = "mulliken"  !! mulliken, chelpg or ignore
      real(dp) :: fmo_resppc = 2.0_dp    !! Point-charge cutoff; negative disables it
      integer :: fmo_max_outer = 50
      real(dp) :: fmo_tolerance = 1.0e-7_dp
      integer :: fmo_scf_max_iter = 100         !! Inner per-fragment SCF iteration cap
      real(dp) :: fmo_scf_energy_tol = 1.0e-9_dp   !! Inner SCF energy convergence
      real(dp) :: fmo_scf_density_tol = 1.0e-7_dp  !! Inner SCF density convergence
      integer :: max_intersection_level = 999  !! Maximum k-way intersection depth for GMBE (default: no limit)
      real(dp), allocatable :: fragment_cutoffs(:)  !! Distance cutoffs for n-mer screening (Angstrom)
      integer :: global_groups = 0
      integer :: nodes_per_group = 0

      ! Calculation-specific keywords (structured)
      type(hessian_keywords_t) :: hessian  !! Hessian calculation keywords
      type(aimd_keywords_t) :: aimd        !! AIMD calculation keywords
      type(optimizer_settings_t) :: optimization  !! Geometry optimization keywords
      ! TODO(mqc): `scf` is written by `config_to_driver` and read nowhere in
      ! the tree; every backend takes its SCF settings from
      ! `method_config%scf`. Two structures named scf, one of them wired up.
      type(scf_keywords_t) :: scf          !! SCF calculation keywords

      character(len=256), allocatable :: fragment_potentials(:)
         !! The effective fragment potential describing each fragment, in
         !! fragment order and empty where a fragment carries none. Fixed-length
         !! because this is a flat array of paths and a deferred-length one
         !! cannot be.

      ! Output control
      logical :: skip_json_output = .false.  !! Skip JSON output for large calculations
      logical :: unchecked_input = .false.   !! Validation warns instead of refusing
      character(len=16) :: fragment_breakdown = "csv"  !! Per-fragment table: csv, json or none
   end type driver_config_t

contains

   subroutine config_to_driver(mqc_config, driver_config, molecule_index, node_rank, &
                               n_fragments, error)
      !! Convert mqc_config_t to minimal driver_config_t
      !! Extracts only the fields needed by the driver
      !! If molecule_index is provided, uses that molecule's fragment count
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(out) :: driver_config
      integer, intent(in), optional :: molecule_index  !! Which molecule to use (for multi-molecule mode)
      integer, intent(in), optional :: node_rank  !! Node-local MPI rank, for GPU binding
      type(error_t), intent(inout), optional :: error
         !! Refuses a setting spelled in a way this program does not know.
         !! Optional, so the callers which cannot fail are unchanged; an
         !! unrecognised spelling then reaches `run_geometry_optimization`,
         !! which refuses it there.
      integer, intent(in), optional :: n_fragments
         !! How many fragments the system actually has, when the config does
         !! not say. A settings-only document carries no molecules, so its own
         !! `nfrag` is 0, and the rule below would read that as "unfragmented"
         !! and compute the whole system as one piece -- which converges and
         !! looks right. Pass the count from wherever the geometry came from.

      integer :: nfrag_to_use

      ! Build method configuration
      driver_config%method_config%method_type = mqc_config%method
      ! **True, and not derived from the log level.** This flag answers *whose*
      ! calculation this is, not how much of it to print: an unfragmented run
      ! IS the calculation the user asked for, so its method may speak. How
      ! much it says is the level's business -- the iteration table emits at
      ! `info`, the energy components at `large_info`.
      !
      ! It was pinned `.false.` here, deferring to the fragment scheduler, which
      ! meant an unfragmented run reached no scheduler and every method's output
      ! was unreachable at every level. `E(CASSCF)` and `E(CASCI)` were
      ! invisible for that reason, and so was the SCF iteration table.
      !
      ! Inner SCFs that this run owns -- an atomic guess, a basis-ladder rung, a
      ! Fukui ion -- pass `.false.` themselves at their own call sites, which is
      ! the only thing a flag can express and a level cannot.
      driver_config%method_config%verbose = .true.

      ! Node-local rank, so several ranks on one node land on distinct GPUs.
      if (present(node_rank)) driver_config%method_config%device_rank = node_rank

      ! Basis sets. Ignored by the semi-empirical methods, required by HF/DFT.
      ! The auxiliary set is not optional for the cuEST backend: J and K are
      ! always density-fitted there, so leaving it unset would fail at run time
      ! rather than silently fall back.
      if (allocated(mqc_config%basis)) then
         driver_config%method_config%basis_set = mqc_config%basis
      end if
      driver_config%method_config%scf%cartesian = mqc_config%cartesian
      if (allocated(mqc_config%ecp)) then
         driver_config%method_config%ecp_set = mqc_config%ecp
      end if
      if (allocated(mqc_config%aux_basis)) then
         driver_config%method_config%scf%aux_basis_set = mqc_config%aux_basis
         driver_config%method_config%scf%aux_basis_named = .true.
      end if
      if (allocated(mqc_config%functional)) then
         driver_config%method_config%dft%functional = mqc_config%functional
      end if

      ! Configure XTB solvation settings
      call driver_config%method_config%xtb%configure( &
         use_cds=mqc_config%use_cds, &
         use_shift=mqc_config%use_shift, &
         dielectric=mqc_config%dielectric, &
         cpcm_nang=mqc_config%cpcm_nang, &
         cpcm_rscale=mqc_config%cpcm_rscale, &
         solvent=mqc_config%solvent, &
         solvation_model=mqc_config%solvation_model)

      ! Copy calc_type
      driver_config%calc_type = mqc_config%calc_type

      ! Determine fragment count
      if (present(n_fragments)) then
         nfrag_to_use = n_fragments
      else if (present(molecule_index)) then
         ! Multi-molecule mode: use specific molecule's fragment count
         if (molecule_index < 1 .or. molecule_index > mqc_config%nmol) then
            nfrag_to_use = 0
         else
            nfrag_to_use = mqc_config%molecules(molecule_index)%nfrag
         end if
      else
         ! Single molecule mode (backward compatible)
         nfrag_to_use = mqc_config%nfrag
      end if

      ! Carry each fragment's potential through, in fragment order. Blank where a
      ! fragment has none, which is a quantum fragment rather than a missing entry.
      call copy_fragment_potentials(mqc_config, driver_config, molecule_index)

      ! Set fragmentation level
      ! For unfragmented calculations (nfrag=0), nlevel must be 0
      if (nfrag_to_use == 0) then
         driver_config%nlevel = 0
      else
         driver_config%nlevel = mqc_config%frag_level
      end if

      ! **`keywords.fragmentation.method` decides the expansion.** It used to be
      ! required, validated for presence and then never read -- the choice came
      ! from `expansion` and `allow_overlapping_fragments` instead, so every
      ! deck in the tree said "MBE" including the ones running FMO, and
      ! `method: "wibble"` was accepted just as happily.
      !
      ! The two older keys still work and now have to AGREE. Refusing a
      ! contradiction rather than picking a winner: a deck that says both is a
      ! deck someone half-migrated, and whichever lost would lose silently.
      ! `error` is optional here, so the resolver takes it optionally too and
      ! the refusal is only raised for a caller that asked to hear about it.
      ! Passing it unconditionally referenced an absent dummy and segfaulted a
      ! test that had nothing to do with fragmentation.
      call resolve_fragmentation_method(mqc_config, driver_config, error)
      if (present(error)) then
         if (error%has_error()) return
      end if
      if (allocated(mqc_config%counterpoise)) then
         driver_config%counterpoise = mqc_config%counterpoise
      end if
      if (allocated(mqc_config%embedding)) then
         driver_config%embedding = mqc_config%embedding
      end if
      if (allocated(mqc_config%bond_breaking)) then
         driver_config%bond_breaking = mqc_config%bond_breaking
      end if
      driver_config%cap_scale = mqc_config%cap_scale
      ! `expansion_kind` is NOT copied from the deck here any more --
      ! `resolve_fragmentation_method` above derives it from `method` and has
      ! already refused a deck whose `expansion` disagrees. Copying it again
      ! would put the deprecated key back in charge.
      if (allocated(mqc_config%fmo_far_field)) then
         driver_config%fmo_far_field = mqc_config%fmo_far_field
      end if
      driver_config%fmo_resppc = mqc_config%fmo_resppc
      driver_config%fmo_max_outer = mqc_config%fmo_max_outer
      driver_config%fmo_tolerance = mqc_config%fmo_tolerance
      driver_config%fmo_scf_max_iter = mqc_config%fmo_scf_max_iter
      driver_config%fmo_scf_energy_tol = mqc_config%fmo_scf_energy_tol
      driver_config%fmo_scf_density_tol = mqc_config%fmo_scf_density_tol

      ! Set GMBE maximum intersection level
      driver_config%max_intersection_level = mqc_config%max_intersection_level

      ! Copy fragment distance cutoffs if present
      if (allocated(mqc_config%fragment_cutoffs)) then
         allocate (driver_config%fragment_cutoffs(size(mqc_config%fragment_cutoffs)))
         driver_config%fragment_cutoffs = mqc_config%fragment_cutoffs
      end if

      driver_config%global_groups = mqc_config%global_groups
      driver_config%nodes_per_group = mqc_config%nodes_per_group

      ! Set calculation-specific keywords
      driver_config%hessian%displacement = mqc_config%hessian_displacement
      ! ...and to the method, which is what actually runs the displacements.
      driver_config%method_config%hessian_displacement = mqc_config%hessian_displacement
      driver_config%method_config%hessian_response_tol = mqc_config%hessian_response_tol
      driver_config%method_config%hessian_response_max_iter = &
         mqc_config%hessian_response_max_iter
      driver_config%method_config%hessian_response_batch = mqc_config%hessian_response_batch
      driver_config%hessian%temperature = mqc_config%hessian_temperature
      driver_config%hessian%pressure = mqc_config%hessian_pressure
      driver_config%aimd%dt = mqc_config%aimd_dt
      driver_config%aimd%nsteps = mqc_config%aimd_nsteps
      driver_config%aimd%initial_temperature = mqc_config%aimd_initial_temperature
      driver_config%aimd%output_frequency = mqc_config%aimd_output_frequency
      driver_config%optimization%max_steps = mqc_config%opt_max_steps
      driver_config%optimization%gradient_tolerance = mqc_config%opt_gradient_tolerance
      driver_config%optimization%energy_tolerance = mqc_config%opt_energy_tolerance
      driver_config%optimization%max_step = mqc_config%opt_max_step
      driver_config%optimization%lbfgs_memory = mqc_config%opt_lbfgs_memory
      driver_config%optimization%print_level = mqc_config%opt_print_level
      driver_config%optimization%trajectory = mqc_config%opt_trajectory
      driver_config%optimization%freeze_terms = mqc_config%opt_freeze_terms
      driver_config%optimization%hess_end = mqc_config%opt_hess_end
      driver_config%optimization%timestep = mqc_config%opt_timestep
      driver_config%optimization%friction = mqc_config%opt_friction
      driver_config%optimization%friction_factor = mqc_config%opt_friction_factor
      driver_config%optimization%friction_rising = mqc_config%opt_friction_rising
      if (allocated(mqc_config%opt_frozen_atoms)) then
         driver_config%optimization%frozen_atoms = mqc_config%opt_frozen_atoms
      end if
      call set_optimization_constraints(mqc_config, driver_config)
      call set_optimization_vocabulary(mqc_config, driver_config, error)

      driver_config%scf%max_iterations = mqc_config%scf_maxiter
      driver_config%scf%convergence_threshold = mqc_config%scf_tolerance
      ! And again into the method's own SCF settings, which is the copy every
      ! backend actually reads -- `configure_hf` and the xTB calculator both
      ! take max_iter and the tolerances from `method_config%scf`.
      driver_config%method_config%scf%max_iter = mqc_config%scf_maxiter
      driver_config%method_config%scf%max_iter_set = mqc_config%scf_maxiter_set
      driver_config%method_config%scf%energy_convergence = mqc_config%scf_tolerance
      driver_config%method_config%scf%density_convergence = mqc_config%scf_density_tolerance
      driver_config%method_config%scf%gradient_convergence = mqc_config%scf_gradient_tolerance
      driver_config%method_config%scf%level_shift = mqc_config%scf_level_shift
      driver_config%method_config%scf%use_diis = mqc_config%scf_use_diis
      driver_config%method_config%scf%diis_size = mqc_config%scf_diis_size
      driver_config%method_config%scf%incremental_fock = mqc_config%scf_incremental_fock
      if (allocated(mqc_config%scf_accelerator)) then
         driver_config%method_config%scf%accelerator = mqc_config%scf_accelerator
      end if
      if (allocated(mqc_config%scf_convergence_metric)) then
         driver_config%method_config%scf%convergence_metric = &
            mqc_config%scf_convergence_metric
      end if
      driver_config%method_config%scf%linear_dependence = mqc_config%scf_linear_dependence
      ! Carry across whether the deck actually named them, not just what they
      ! came out as: MAKEFP, whose own default is stricter, cannot otherwise
      ! tell "the user wants 1e-6" from "nobody said anything".
      driver_config%method_config%scf%energy_convergence_set = mqc_config%scf_tolerance_set
      driver_config%method_config%scf%density_convergence_set = &
         mqc_config%scf_density_tolerance_set
      ! The MAKEFP group. No presence flags here: these have one default each,
      ! named once in `mqc_calculation_defaults` and read by both the deck and
      ! the solver, so a silent deck carries the solver's own number back to it.
      driver_config%method_config%efp%dynamic_tolerance = mqc_config%efp_dynamic_tolerance
      driver_config%method_config%efp%dynamic_maxiter = mqc_config%efp_dynamic_maxiter
      driver_config%method_config%efp%allow_crap_response = mqc_config%efp_allow_crap_response
      driver_config%method_config%efp%response_batch = mqc_config%efp_response_batch
      driver_config%method_config%efp%response = mqc_config%efp_response
      driver_config%method_config%efp%vdw_scale = mqc_config%efp_vdw_scale
      driver_config%method_config%efp%quadrupole_blocks = mqc_config%efp_quadrupole_blocks
      ! EFMO. `rcut` comes from `keywords.fragmentation` and the switch from
      ! `keywords.efmo`; both land in one object here so the backend reads one
      ! place.
      driver_config%method_config%efmo%rcut = mqc_config%efmo_rcut
      driver_config%method_config%efmo%charge_transfer = mqc_config%efmo_charge_transfer
      ! Quantum nuclei. The two lists are exclusive and one of them is absent.
      driver_config%method_config%neo%active = mqc_config%neo_active
      driver_config%method_config%neo%nuclear_basis = mqc_config%neo_nuclear_basis
      driver_config%method_config%neo%epc = mqc_config%neo_epc
      if (allocated(mqc_config%neo_quantum_indices)) then
         driver_config%method_config%neo%quantum_indices = mqc_config%neo_quantum_indices
      end if
      if (allocated(mqc_config%neo_quantum_symbols)) then
         driver_config%method_config%neo%quantum_symbols = mqc_config%neo_quantum_symbols
      end if
      if (allocated(mqc_config%checkpoint_file)) then
         driver_config%checkpoint_file = mqc_config%checkpoint_file
      end if
      driver_config%method_config%scf%allow_crap_scf = mqc_config%allow_crap_scf
      driver_config%method_config%scf%unrestricted = mqc_config%unrestricted
      driver_config%method_config%scf%density_fitting = mqc_config%scf_density_fitting
      driver_config%method_config%corr%freeze_core = mqc_config%corr_freeze_core
      driver_config%method_config%corr%n_frozen_core = mqc_config%corr_n_frozen_core
      driver_config%method_config%corr%use_df = mqc_config%corr_density_fitting
      driver_config%method_config%corr%use_scs = mqc_config%corr_scs
      driver_config%method_config%corr%scs_ss = mqc_config%corr_scs_ss
      driver_config%method_config%corr%scs_os = mqc_config%corr_scs_os
      driver_config%method_config%dft%grid_level = mqc_config%dft_grid_level
      driver_config%method_config%dft%nlc_grid_level = &
         mqc_config%dft_nlc_grid_level
      driver_config%method_config%dft%screening_tolerance = &
         mqc_config%dft_screening_tolerance
      driver_config%method_config%dft%block_size = mqc_config%dft_block_size
      driver_config%method_config%backend = mqc_config%backend
      if (allocated(mqc_config%charges_scheme)) then
         driver_config%method_config%properties%charges_scheme = &
            mqc_config%charges_scheme
      end if
      if (allocated(mqc_config%bonding_analysis)) then
         driver_config%method_config%properties%bonding_analysis = &
            mqc_config%bonding_analysis
      end if
      driver_config%method_config%properties%bonding_threshold = &
         mqc_config%bonding_threshold
      driver_config%method_config%properties%bonding_energy = &
         mqc_config%bonding_energy
      if (allocated(mqc_config%fukui_population)) then
         driver_config%method_config%properties%fukui_population = &
            mqc_config%fukui_population
      end if
      if (allocated(mqc_config%fukui_guess)) then
         driver_config%method_config%properties%fukui_guess = mqc_config%fukui_guess
      end if
      ! One assignment, and no test: the reader already resolved this against
      ! `keywords.scf`, so every field holds the value the ions should use.
      driver_config%method_config%properties%fukui_scf = mqc_config%fukui_scf
      driver_config%method_config%properties%bonding_no_sharing = &
         mqc_config%bonding_no_sharing
      driver_config%method_config%properties%bonding_restrict_localization = &
         mqc_config%bonding_restrict_localization
      if (allocated(mqc_config%bonding_no_sharing_ci)) then
         driver_config%method_config%properties%bonding_no_sharing_ci = &
            mqc_config%bonding_no_sharing_ci
      end if

      driver_config%method_config%pcm%enabled = mqc_config%pcm_enabled
      driver_config%method_config%pcm%method = mqc_config%pcm_method
      driver_config%method_config%pcm%dielectric = mqc_config%pcm_dielectric
      driver_config%method_config%pcm%angular_points = mqc_config%pcm_angular_points
      driver_config%method_config%pcm%radii_scale = mqc_config%pcm_radii_scale
      driver_config%method_config%pcm%zeta = mqc_config%pcm_zeta
      driver_config%method_config%pcm%tolerance = mqc_config%pcm_tolerance
      driver_config%method_config%pcm%max_iter = mqc_config%pcm_max_iter
      ! Only overridden when a deck asked; -1 leaves the level in charge, and
      ! the grid builder refuses one without the other rather than half-applying.
      if (mqc_config%dft_radial_points > 0 .and. mqc_config%dft_angular_points > 0) then
         driver_config%method_config%dft%radial_points = mqc_config%dft_radial_points
         driver_config%method_config%dft%angular_points = mqc_config%dft_angular_points
         ! Negative level means "the counts are in charge", which is what the
         ! banner and the grid builder both read to decide which was asked for.
         driver_config%method_config%dft%grid_level = -1
      end if
      driver_config%method_config%cc%max_iter = mqc_config%cc_maxiter
      driver_config%method_config%cc%amplitude_convergence = mqc_config%cc_tolerance
      driver_config%method_config%cc%use_diis = mqc_config%cc_diis
      driver_config%method_config%cc%diis_size = mqc_config%cc_diis_size
      driver_config%method_config%cc%spin_adapted = mqc_config%cc_spin_adapted
      ! The method name settles the triples unless a deck said otherwise:
      ! "ccsd(t)" and "ccsd" are separate method types, so the distinction
      ! survives the parse, unlike RI.
      driver_config%method_config%cc%include_triples = &
         (mqc_config%method == METHOD_TYPE_CCSD_T)
      if (mqc_config%cc_triples_set) then
         driver_config%method_config%cc%include_triples = mqc_config%cc_triples
      end if
      if (allocated(mqc_config%mcscf_avas_orbitals)) then
         driver_config%method_config%mcscf%avas_orbitals = mqc_config%mcscf_avas_orbitals
      end if
      driver_config%method_config%mcscf%avas_threshold = mqc_config%mcscf_avas_threshold
      driver_config%method_config%mcscf%full_valence = mqc_config%mcscf_full_valence
      if (allocated(mqc_config%mcscf_ormas_subspaces)) then
         driver_config%method_config%mcscf%ormas_subspaces = &
            mqc_config%mcscf_ormas_subspaces
         driver_config%method_config%mcscf%ormas_min_electrons = &
            mqc_config%mcscf_ormas_min_electrons
         driver_config%method_config%mcscf%ormas_max_electrons = &
            mqc_config%mcscf_ormas_max_electrons
      end if
      driver_config%method_config%mcscf%n_active_electrons = &
         mqc_config%mcscf_n_active_electrons
      driver_config%method_config%mcscf%n_active_orbitals = &
         mqc_config%mcscf_n_active_orbitals
      driver_config%method_config%mcscf%n_inactive_orbitals = &
         mqc_config%mcscf_n_inactive_orbitals
      ! Already settled against the method spelling by the reader, which is the
      ! only place that spelling exists -- unlike the triples above, whose two
      ! method types survive the parse and so are resolved here.
      driver_config%method_config%mcscf%optimize_orbitals = &
         mqc_config%mcscf_optimize_orbitals
      driver_config%method_config%mcscf%max_macro_iter = mqc_config%mcscf_max_macro_iter
      driver_config%method_config%mcscf%orbital_convergence = &
         mqc_config%mcscf_orbital_convergence
      ! `guess_type` is the current spelling and `scf_guess` the superseded one.
      ! The schema has already refused a deck that sets both, so whichever is
      ! allocated is the one the user meant.
      if (allocated(mqc_config%scf_guess)) then
         driver_config%method_config%scf%guess = mqc_config%scf_guess
      end if
      if (allocated(mqc_config%guess_type)) then
         driver_config%method_config%scf%guess = mqc_config%guess_type
      end if
      if (allocated(mqc_config%guess_steps)) then
         driver_config%method_config%scf%guess_steps = mqc_config%guess_steps
      end if
      if (allocated(mqc_config%scf_eri_path)) then
         driver_config%method_config%scf%eri_path = mqc_config%scf_eri_path
      end if

      ! Output control
      driver_config%skip_json_output = mqc_config%skip_json_output
      driver_config%unchecked_input = mqc_config%unchecked_input
      if (allocated(mqc_config%fragment_breakdown)) then
         driver_config%fragment_breakdown = mqc_config%fragment_breakdown
      end if

   end subroutine config_to_driver

   subroutine set_optimization_vocabulary(mqc_config, driver_config, error)
      !! Turn the two spelled optimization settings into constants
      !!
      !! The vocabulary lives with the type that holds the constants rather
      !! than in the reader. An unrecognised spelling is refused here by name,
      !! this being the last place that still has the string the user typed.
      use pic_io, only: to_char
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(inout) :: driver_config
      type(error_t), intent(inout), optional :: error

      integer :: want

      if (allocated(mqc_config%opt_coordinates)) then
         driver_config%optimization%coordinates = &
            coordinates_from_string(mqc_config%opt_coordinates)
         if (driver_config%optimization%coordinates == OPT_COORDS_UNKNOWN) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "Unknown keywords.optimization.coordinates: '"// &
                              trim(mqc_config%opt_coordinates)// &
                              "'. Use cartesian, hdlc or dlc.")
            end if
            return
         end if
      end if

      if (allocated(mqc_config%opt_algorithm)) then
         driver_config%optimization%algorithm = &
            algorithm_from_string(mqc_config%opt_algorithm)
         if (driver_config%optimization%algorithm == OPT_ALGO_UNKNOWN) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "Unknown keywords.optimization.algorithm: '"// &
                              trim(mqc_config%opt_algorithm)// &
                              "'. Use lbfgs, cg, cg-auto, sd, prfo, nr or damped.")
            end if
            return
         end if
      end if

      ! "auto" is a spelling of "leave DL-FIND's own choice alone", so it maps
      ! to the sentinel rather than to a scheme -- which is why the refusal
      ! below tests for the parse failure and not for the sentinel.
      if (allocated(mqc_config%opt_hessian_update)) then
         driver_config%optimization%hessian_update = &
            hessian_update_from_string(mqc_config%opt_hessian_update)
         if (driver_config%optimization%hessian_update < OPT_HESSIAN_UPDATE_ENGINE) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "Unknown keywords.optimization.hessian_update: '"// &
                              trim(mqc_config%opt_hessian_update)// &
                              "'. Use none, powell, bofill or auto.")
            end if
            return
         end if
      end if

      ! Parsed into a local and only then stored: `error` is optional, and a
      ! caller that omits it would otherwise carry the parse failure into the
      ! driver, where a sentinel below OPT_TARGET_MINIMUM reads as "not a
      ! saddle" and runs a minimisation.
      if (allocated(mqc_config%opt_target)) then
         want = target_from_string(mqc_config%opt_target)
         if (want < OPT_TARGET_MINIMUM) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "Unknown keywords.optimization.target: '"// &
                              trim(mqc_config%opt_target)// &
                              "'. Use minimum or saddle.")
            end if
            return
         end if
         driver_config%optimization%target = want
      end if

      if (allocated(mqc_config%opt_endpoint)) then
         driver_config%optimization%endpoint = mqc_config%opt_endpoint
      end if
      driver_config%optimization%n_images = mqc_config%opt_n_images
      driver_config%optimization%neb_spring = mqc_config%opt_neb_spring

      if (allocated(mqc_config%opt_neb_ends)) then
         driver_config%optimization%neb_ends = neb_ends_from_string(mqc_config%opt_neb_ends)
         if (driver_config%optimization%neb_ends < NEB_ENDS_FROZEN) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "Unknown keywords.optimization.neb_endpoints: '"// &
                              trim(mqc_config%opt_neb_ends)// &
                              "'. Use frozen, perpendicular or free.")
            end if
            return
         end if
      end if

      ! The path keywords describe a band, and a band needs two structures. A
      ! deck that sets images or a spring without an endpoint would otherwise
      ! get the ordinary single-structure optimization, which is different work
      ! from the one asked for.
      if (.not. allocated(driver_config%optimization%endpoint)) then
         if (mqc_config%opt_n_images > 0 .or. allocated(mqc_config%opt_neb_ends) .or. &
             mqc_config%opt_neb_spring >= 0.0_dp) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "keywords.optimization sets a path option (images, "// &
                              "neb_spring or neb_endpoints) but no 'endpoint'. A "// &
                              "chain-of-states run needs a second structure to run "// &
                              "between; give the product geometry as "// &
                              "keywords.optimization.endpoint.")
            end if
            return
         end if
      else
         ! Two images are the endpoints and nothing in between, so there is no
         ! band to relax. DL-FIND refuses it too, but from inside the engine and
         ! with a message about its own variable.
         if (mqc_config%opt_n_images > 0 .and. mqc_config%opt_n_images < MIN_NEB_IMAGES) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "keywords.optimization.images is "// &
                              trim(int_to_text(mqc_config%opt_n_images))// &
                              ". A path needs at least "// &
                              trim(int_to_text(MIN_NEB_IMAGES))//": two of them are the "// &
                              "endpoints, so fewer leaves nothing between them to relax.")
            end if
            return
         end if
      end if

      driver_config%optimization%dimer_separation = mqc_config%opt_dimer_separation
      driver_config%optimization%dimer_max_rotations = mqc_config%opt_dimer_max_rotations
      driver_config%optimization%dimer_rotation_tolerance = mqc_config%opt_dimer_rot_tol

      if (allocated(mqc_config%opt_saddle_method)) then
         driver_config%optimization%saddle_method = &
            saddle_method_from_string(mqc_config%opt_saddle_method)
         if (driver_config%optimization%saddle_method < SADDLE_METHOD_PRFO) then
            if (present(error)) then
               call error%set(ERROR_VALIDATION, &
                              "Unknown keywords.optimization.saddle_method: '"// &
                              trim(mqc_config%opt_saddle_method)// &
                              "'. Use prfo or dimer.")
            end if
            return
         end if
      end if

      driver_config%optimization%zero_modes = mqc_config%opt_zero_modes
      driver_config%optimization%soft_mode_threshold = mqc_config%opt_soft_modes
      driver_config%optimization%connect = mqc_config%opt_connect
      driver_config%optimization%connect_distort = mqc_config%opt_connect_distort

      ! Following a saddle downhill needs a saddle to follow. Asked for
      ! alongside a minimisation, running it anyway would report success for
      ! half of what the deck asked.
      if (mqc_config%opt_connect .and. &
          driver_config%optimization%target /= OPT_TARGET_SADDLE) then
         if (present(error)) then
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization.connect asks for the minima on either "// &
                           "side of a transition state, but target is 'minimum'. Set "// &
                           "target to 'saddle'.")
         end if
         return
      end if

      ! A connect run does not stop at the saddle -- it carries on downhill
      ! twice -- so the geometry it finishes with is the second minimum. A
      ! Hessian there describes that minimum and reports, correctly, that it
      ! has no imaginary frequencies, which would read as the saddle having
      ! been checked and failed.
      if (mqc_config%opt_connect .and. mqc_config%opt_hess_end) then
         if (present(error)) then
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization sets both connect and hess_end. A "// &
                           "connect run ends at one of the minima rather than at the "// &
                           "saddle, so the final Hessian would describe the wrong "// &
                           "structure. Run the saddle search with hess_end first, then "// &
                           "connect from the geometry it found.")
         end if
         return
      end if

      ! A saddle is not something every algorithm can look for. Steepest
      ! descent and the quasi-Newton minimisers go downhill by construction and
      ! report a minimum however the target is spelled, so asking them for a
      ! saddle is a deck that cannot be satisfied. Refused rather than silently
      ! switched to one that can.
      if (driver_config%optimization%target == OPT_TARGET_SADDLE) then
         ! The dimer finds a saddle by construction -- it inverts the force
         ! along the softest direction -- so the algorithm beside it is only
         ! translating and rotating the pair, and a minimiser is the right
         ! thing there. The P-RFO requirement applies to the P-RFO path only.
         if (driver_config%optimization%saddle_method == SADDLE_METHOD_PRFO .and. &
             .not. algorithm_finds_saddle(driver_config%optimization%algorithm)) then
            if (present(error)) then
               ! Named or defaulted, said differently. A deck with no
               ! `algorithm` key at all gets the default, and a message that
               ! quotes it back as though the user had chosen it sends them
               ! looking for a line their file does not contain.
               if (allocated(mqc_config%opt_algorithm)) then
                  call error%set(ERROR_VALIDATION, &
                                 "keywords.optimization.target is 'saddle' but the algorithm is '"// &
                                 algorithm_to_string(driver_config%optimization%algorithm)// &
                                 "', which descends to a minimum. Use prfo, which "// &
                                 "maximises along one mode and minimises along the rest, or nr.")
               else
                  call error%set(ERROR_VALIDATION, &
                                 "keywords.optimization.target is 'saddle' but no algorithm was "// &
                                 "named, so the default '"// &
                                 algorithm_to_string(driver_config%optimization%algorithm)// &
                                 "' applies, and it descends to a minimum. Set "// &
                                 "keywords.optimization.algorithm to prfo, which maximises along "// &
                                 "one mode and minimises along the rest, or nr.")
               end if
            end if
            return
         end if
         ! Cartesian P-RFO maximises along the lowest eigenvalue, and in
         ! Cartesian coordinates the lowest six are translations and rotations.
         ! It follows one of those instead of the reaction mode and wanders.
         ! Warned rather than refused: from a guess close enough that the
         ! reaction mode is already the lowest eigenvalue this can be
         ! satisfied, where a downhill algorithm cannot. Cartesian is also the
         ! default, so a deck that never mentioned coordinates lands here.
         if (driver_config%optimization%coordinates == OPT_COORDS_CARTESIAN .and. &
             driver_config%optimization%saddle_method == SADDLE_METHOD_PRFO .and. &
             driver_config%optimization%zero_modes < 0) then
            if (allocated(mqc_config%opt_coordinates)) then
               call logger%warning("  keywords.optimization: a saddle search in Cartesian "// &
                                   "coordinates follows a rotation, not the reaction mode.")
            else
               call logger%warning("  keywords.optimization: no coordinates were named, so this "// &
                                   "saddle search runs in Cartesians, where it follows a "// &
                                   "rotation rather than the reaction mode.")
            end if
            call logger%warning("     Expect it to wander for all "// &
                                to_char(driver_config%optimization%max_steps)// &
                                " steps without converging.")
            call logger%warning("     Use coordinates 'dlc' (or 'hdlc' for a cluster), or "// &
                                "set zero_modes to 6 (5 if linear) so the")
            call logger%warning("     translations and rotations are recognised as zero "// &
                                "rather than followed.")
         end if
      end if

   end subroutine set_optimization_vocabulary

   subroutine set_optimization_constraints(mqc_config, driver_config)
      !! Gather the constrained coordinates into the optimizer's own type
      !!
      !! The reader has already checked that each type is one this program
      !! knows and that its atom count matches, so this only rearranges the two
      !! parallel arrays into the one array of records the engine iterates.
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(inout) :: driver_config

      integer :: i, n

      if (.not. allocated(mqc_config%opt_constraint_kinds)) return
      n = size(mqc_config%opt_constraint_kinds)
      if (n < 1) return

      allocate (driver_config%optimization%constraints(n))
      do i = 1, n
         driver_config%optimization%constraints(i)%kind = mqc_config%opt_constraint_kinds(i)
         driver_config%optimization%constraints(i)%atoms = mqc_config%opt_constraint_atoms(:, i)
      end do
   end subroutine set_optimization_constraints

   subroutine check_counterpoise_support(driver_config, error)
      !! Refuse a counterpoise request the chosen expansion cannot honour
      !!
      !! Counterpoise is carried by ghosted rows in the MBE term list and by
      !! the `is_ghost` mask on the fragments those rows build. Three things
      !! read neither, and each fails differently:
      !!
      !!   * **GMBE** and **FMO/EE-MBE** build their own term lists and never
      !!     call `generate_mbe_term_list`, so the setting is not seen and the
      !!     answer is a valid uncorrected one -- the number the deck asked
      !!     this program *not* to produce, with nothing about it saying so.
      !!   * the **xTB and EFP** paths take `element_numbers` as the atom list
      !!     and never consult `is_ghost`, so a ghosted monomer is computed as
      !!     though its ghost centres were nuclei, against an electron count
      !!     derived without them. That is a different molecule with the wrong
      !!     charge, not an approximation.
      !!
      !! An unrecognised spelling is refused too: once it reaches the term list
      !! it cannot be told from `none`.
      type(driver_config_t), intent(in) :: driver_config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: scheme

      scheme = trim(driver_config%counterpoise)
      if (scheme == "none" .or. len(scheme) == 0) return

      if (scheme /= "vmfc") then
         call error%set(ERROR_VALIDATION, &
                        "Unknown keywords.fragmentation.counterpoise: '"//scheme// &
                        "'. Use vmfc, or none.")
         return
      end if

      if (driver_config%allow_overlapping_fragments) then
         call error%set(ERROR_VALIDATION, &
                        "counterpoise is not available for GMBE. GMBE builds its "// &
                        "terms by inclusion-exclusion over overlapping primaries "// &
                        "rather than from the subset list counterpoise ghosts, so "// &
                        "the request would be ignored and an uncorrected energy "// &
                        "returned. Use MBE, or drop counterpoise.")
         return
      end if

      if (trim(driver_config%expansion_kind) /= "mbe") then
         call error%set(ERROR_VALIDATION, &
                        "counterpoise is not available for expansion '"// &
                        trim(driver_config%expansion_kind)//"'. An embedded "// &
                        "expansion builds its own term list and would ignore the "// &
                        "request. Use the plain expansion, or drop counterpoise.")
         return
      end if

      select case (driver_config%method_config%method_type)
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_EFP2)
         call error%set(ERROR_VALIDATION, &
                        "counterpoise needs a method with a basis set to ghost. "// &
                        "The semi-empirical and EFP paths take every centre as a "// &
                        "real atom, so a ghosted term would be computed as a "// &
                        "different molecule rather than as a basis-set correction. "// &
                        "Use an ab initio method, or drop counterpoise.")
      case default
         ! An ab initio method: the cenzontle path reads `is_ghost` and honours it.
      end select

   end subroutine check_counterpoise_support

   subroutine copy_fragment_potentials(mqc_config, driver_config, molecule_index)
      !! The per-fragment potential paths, in fragment order
      !!
      !! Left unallocated when no fragment carries one, so the driver can tell
      !! "not an EFP calculation" from "an EFP calculation whose potentials are
      !! all blank" without a second flag.
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(inout) :: driver_config
      integer, intent(in), optional :: molecule_index

      integer :: n, k
      logical :: any_potential

      if (allocated(driver_config%fragment_potentials)) then
         deallocate (driver_config%fragment_potentials)
      end if

      if (present(molecule_index)) then
         if (molecule_index < 1 .or. molecule_index > mqc_config%nmol) return
         n = mqc_config%molecules(molecule_index)%nfrag
         if (n == 0) return
         allocate (driver_config%fragment_potentials(n))
         driver_config%fragment_potentials = ""
         any_potential = .false.
         do k = 1, n
            if (.not. allocated(mqc_config%molecules(molecule_index)%fragments(k)%potential)) cycle
            driver_config%fragment_potentials(k) = &
               mqc_config%molecules(molecule_index)%fragments(k)%potential
            any_potential = .true.
         end do
      else
         n = mqc_config%nfrag
         if (n == 0) return
         allocate (driver_config%fragment_potentials(n))
         driver_config%fragment_potentials = ""
         any_potential = .false.
         do k = 1, n
            if (.not. allocated(mqc_config%fragments(k)%potential)) cycle
            driver_config%fragment_potentials(k) = mqc_config%fragments(k)%potential
            any_potential = .true.
         end do
      end if

      if (.not. any_potential) deallocate (driver_config%fragment_potentials)
   end subroutine copy_fragment_potentials

   pure function fragmentation_allows_overlap(mqc_config) result(yes)
      !! Whether primaries may overlap, derived from `method` alone
      !!
      !! **Two consumers, one derivation.** `config_to_system_geometry` refuses
      !! overlapping fragments on its own path, reading this rather than a deck
      !! flag -- when `method: gmbe` first started deriving the switch, only
      !! the driver got it and a GMBE deck was refused at the door by the
      !! geometry, telling it to set the key the migration had just removed.
      use mqc_fragmentation_method, only: parse_fragmentation_method, &
                                          fragmentation_method_overlapping
      type(mqc_config_t), intent(in) :: mqc_config
      logical :: yes

      integer :: method
      logical :: ok

      yes = .false.
      if (.not. allocated(mqc_config%frag_method)) return
      call parse_fragmentation_method(mqc_config%frag_method, method, ok)
      if (ok) yes = fragmentation_method_overlapping(method)
   end function fragmentation_allows_overlap

   subroutine resolve_fragmentation_method(mqc_config, driver_config, error)
      !! Turn `keywords.fragmentation.method` into the two internal switches
      !!
      !! `expansion_kind` and `allow_overlapping_fragments` are what the driver
      !! and the term builders read; this is the only place either is decided.
      use mqc_fragmentation_method, only: parse_fragmentation_method, &
                                          fragmentation_method_name, &
                                          fragmentation_method_expansion, &
                                          fragmentation_method_overlapping, &
                                          fragmentation_method_implemented, &
                                          fragmentation_method_list
      type(mqc_config_t), intent(in) :: mqc_config
      type(driver_config_t), intent(inout) :: driver_config
      type(error_t), intent(inout), optional :: error
         !! Optional, matching `config_to_driver`. A caller that does not want
         !! to hear about a bad method still gets the derived switches for a
         !! good one, and silently keeps the defaults for a bad one -- which is
         !! what every other refusal in this file does.

      integer :: method
      logical :: ok
      character(len=16) :: derived_kind
      logical :: derived_overlap

      ! Unfragmented runs never name it, and the schema only requires it inside
      ! a `fragmentation` block, so an absent method is not an error here.
      if (.not. allocated(mqc_config%frag_method)) return

      call parse_fragmentation_method(mqc_config%frag_method, method, ok)
      if (.not. ok) then
         if (present(error)) call error%set(ERROR_VALIDATION, &
                                            "Unknown keywords.fragmentation.method: '"// &
                                            trim(mqc_config%frag_method)//"'. Use one of "// &
                                            fragmentation_method_list()//".")
         return
      end if
      if (.not. fragmentation_method_implemented(method)) then
         if (present(error)) call error%set(ERROR_VALIDATION, &
                                            "keywords.fragmentation.method '"// &
                                            fragmentation_method_name(method)//"' is a method this "// &
                                            "program knows the name of and does not implement yet. "// &
                                            "Use one of "//fragmentation_method_list()//".")
         return
      end if

      derived_kind = fragmentation_method_expansion(method)
      derived_overlap = fragmentation_method_overlapping(method)

      driver_config%expansion_kind = derived_kind
      driver_config%allow_overlapping_fragments = derived_overlap
   end subroutine resolve_fragmentation_method

   subroutine config_to_system_geometry(mqc_config, sys_geom, error, molecule_index)
      !! Convert mqc_config_t geometry to system_geometry_t
      !! For unfragmented calculations (nfrag=0), treats entire system as single unit
      !! For fragmented calculations, currently assumes monomer-based fragmentation
      !! If molecule_index is provided, uses that specific molecule from multi-molecule mode
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      type(error_t), intent(out) :: error
      integer, intent(in), optional :: molecule_index  !! Which molecule to use (for multi-molecule mode)

      ! TODO(mqc): `i` is dead here.
      integer :: i
      logical :: use_angstrom

      ! Determine units
      use_angstrom = .true.
      if (allocated(mqc_config%units)) then
         if (trim(mqc_config%units) == "bohr") then
            use_angstrom = .false.
         end if
      end if

      ! Handle multi-molecule vs single molecule mode
      if (present(molecule_index)) then
         ! Multi-molecule mode: extract specific molecule
         if (molecule_index < 1 .or. molecule_index > mqc_config%nmol) then
            call error%set(ERROR_VALIDATION, "Invalid molecule_index in multi-molecule mode")
            return
         end if
         call molecule_to_system_geometry(mqc_config%molecules(molecule_index), &
                                          sys_geom, use_angstrom, &
                                          fragmentation_allows_overlap(mqc_config), error)
      else
         ! Single molecule mode (backward compatible)
         ! Check if geometry is loaded
         if (mqc_config%geometry%natoms == 0) then
            call error%set(ERROR_VALIDATION, "No geometry loaded in mqc_config")
            return
         end if

         if (mqc_config%nfrag == 0) then
            ! Unfragmented calculation: entire system is one "monomer"
            call geometry_to_system_unfragmented(mqc_config%geometry, sys_geom, use_angstrom)
            sys_geom%charge = mqc_config%charge
            sys_geom%multiplicity = mqc_config%multiplicity
         else
            ! Fragmented calculation with explicit fragments
            call geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, error)
            if (error%has_error()) then
               call error%add_context("mqc_config_adapter:config_to_system_geometry")
               return
            end if
         end if
      end if

      ! Set after both branches so it survives either path. A capping convention
      ! belongs to the run rather than to a fragment: every fragment built from
      ! this geometry inherits it, which is what caps a cut bond the same way in
      ! a monomer and in the dimer containing it, so the cap contributions
      ! cancel through the expansion.
      sys_geom%cap_scale = mqc_config%cap_scale

   end subroutine config_to_system_geometry

   subroutine geometry_to_system_unfragmented(geom, sys_geom, use_angstrom)
      !! Convert geometry to system_geometry_t for unfragmented calculation
      !! Treats entire system as a single monomer
      use mqc_geometry, only: geometry_type

      type(geometry_type), intent(in) :: geom
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom

      integer :: i

      ! For unfragmented: n_monomers=1, atoms_per_monomer=natoms
      sys_geom%n_monomers = 1
      sys_geom%atoms_per_monomer = geom%natoms
      sys_geom%total_atoms = geom%natoms

      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(geom%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(geom%coords)
      else
         sys_geom%coordinates = geom%coords
      end if

   end subroutine geometry_to_system_unfragmented

   subroutine initialize_fragmented_system(nfrag, geom, fragments, charge, multiplicity, &
                                           allow_overlapping, use_angstrom, sys_geom, error)
      !! Shared helper to initialize system_geometry_t for fragmented calculations
      !! Handles fragment allocation, size checking, and overlap validation
      use mqc_geometry, only: geometry_type
      use mqc_config_types, only: input_fragment_t

      integer, intent(in) :: nfrag
      type(geometry_type), intent(in) :: geom
      type(input_fragment_t), intent(in) :: fragments(:)
      integer, intent(in) :: charge, multiplicity
      logical, intent(in) :: allow_overlapping
      logical, intent(in) :: use_angstrom
      type(system_geometry_t), intent(out) :: sys_geom
      type(error_t), intent(out) :: error

      integer :: i, j, atoms_in_first_frag, max_frag_size
      logical :: all_same_size

      ! Set up basic system geometry
      sys_geom%n_monomers = nfrag
      sys_geom%total_atoms = geom%natoms
      sys_geom%charge = charge
      sys_geom%multiplicity = multiplicity

      ! Allocate fragment info arrays
      allocate (sys_geom%fragment_sizes(nfrag))
      allocate (sys_geom%fragment_charges(nfrag))
      allocate (sys_geom%fragment_multiplicities(nfrag))

      ! Get fragment sizes
      max_frag_size = 0
      atoms_in_first_frag = size(fragments(1)%indices)
      all_same_size = .true.

      do i = 1, nfrag
         sys_geom%fragment_sizes(i) = size(fragments(i)%indices)
         sys_geom%fragment_charges(i) = fragments(i)%charge
         sys_geom%fragment_multiplicities(i) = fragments(i)%multiplicity
         max_frag_size = max(max_frag_size, sys_geom%fragment_sizes(i))
         if (sys_geom%fragment_sizes(i) /= atoms_in_first_frag) then
            all_same_size = .false.
         end if
      end do

      ! Allocate fragment_atoms array
      allocate (sys_geom%fragment_atoms(max_frag_size, nfrag))
      sys_geom%fragment_atoms = -1  ! Initialize with invalid index

      ! Store fragment atom indices (0-indexed from input file)
      do i = 1, nfrag
         do j = 1, sys_geom%fragment_sizes(i)
            sys_geom%fragment_atoms(j, i) = fragments(i)%indices(j)
         end do
      end do

      ! Check for overlapping fragments if not allowed
      if (.not. allow_overlapping) then
         call check_fragment_overlap(fragments, nfrag, error)
         if (error%has_error()) then
            call error%add_context("mqc_config_adapter:geometry_to_system_fragmented")
            return
         end if
      end if

      ! Set atoms_per_monomer: use common size if identical, else 0
      if (all_same_size) then
         sys_geom%atoms_per_monomer = atoms_in_first_frag
      else
         sys_geom%atoms_per_monomer = 0  ! Signal variable-sized fragments
      end if

      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(geom%elements(i))
      end do

      ! Store coordinates (convert to Bohr if needed)
      if (use_angstrom) then
         sys_geom%coordinates = to_bohr(geom%coords)
      else
         sys_geom%coordinates = geom%coords
      end if

   end subroutine initialize_fragmented_system

   subroutine geometry_to_system_fragmented(mqc_config, sys_geom, use_angstrom, error)
      !! Convert geometry to system_geometry_t for fragmented calculation
      !! Supports both identical and variable-sized fragments
      type(mqc_config_t), intent(in) :: mqc_config
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom
      type(error_t), intent(out) :: error

      call initialize_fragmented_system(mqc_config%nfrag, mqc_config%geometry, mqc_config%fragments, &
                                        mqc_config%charge, mqc_config%multiplicity, &
                                        fragmentation_allows_overlap(mqc_config), use_angstrom, &
                                        sys_geom, error)

   end subroutine geometry_to_system_fragmented

   subroutine molecule_to_system_geometry(mol, sys_geom, use_angstrom, allow_overlapping, error)
      !! Convert a molecule_t to system_geometry_t
      !! Handles both unfragmented (nfrag=0) and fragmented molecules
      use mqc_config_types, only: molecule_t

      type(molecule_t), intent(in) :: mol
      type(system_geometry_t), intent(out) :: sys_geom
      logical, intent(in) :: use_angstrom
      type(error_t), intent(out) :: error
      logical, intent(in) :: allow_overlapping

      ! Check if geometry is loaded
      if (mol%geometry%natoms == 0) then
         call error%set(ERROR_VALIDATION, "No geometry loaded in molecule")
         return
      end if

      if (mol%nfrag == 0) then
         ! Unfragmented molecule
         call geometry_to_system_unfragmented(mol%geometry, sys_geom, use_angstrom)
         sys_geom%charge = mol%charge
         sys_geom%multiplicity = mol%multiplicity
      else
         ! Fragmented molecule
         call initialize_fragmented_system(mol%nfrag, mol%geometry, mol%fragments, &
                                           mol%charge, mol%multiplicity, &
                                           allow_overlapping, use_angstrom, &
                                           sys_geom, error)
      end if

   end subroutine molecule_to_system_geometry

   function get_logger_level(level_string) result(level_int)
      !! Convert string log level to integer value
      !! This function uses the pic_logger constants
      use pic_logger, only: debug_level, verbose_level, large_info_level, info_level, &
                            performance_level, warning_level, error_level, knowledge_level
      character(len=*), intent(in) :: level_string
      integer :: level_int

      select case (trim(adjustl(level_string)))
      case ("debug", "Debug", "DEBUG")
         level_int = debug_level
      case ("verbose", "Verbose", "VERBOSE")
         level_int = verbose_level
      case ("large_info", "Large_Info", "LARGE_INFO")
         level_int = large_info_level
      case ("info", "Info", "INFO")
         level_int = info_level
      case ("performance", "Performance", "PERFORMANCE")
         level_int = performance_level
      case ("warning", "Warning", "WARNING")
         level_int = warning_level
      case ("error", "Error", "ERROR")
         level_int = error_level
      case ("knowledge", "Knowledge", "KNOWLEDGE")
         level_int = knowledge_level
      case default
         ! Default to info level if unknown
         call logger%warning("Unknown log level string: "//level_string//". Defaulting to INFO level.")
         level_int = info_level
      end select
   end function get_logger_level

   subroutine check_fragment_overlap(fragments, nfrag, error)
      !! Check if any atoms appear in multiple fragments
      ! TODO(mqc): every pair of fragments against every pair of their atoms,
      ! so `O(nfrag^2 * natoms_per_frag^2)` and not the `O(nfrag *
      ! natoms_per_frag^2)` this used to claim. A single pass marking each atom
      ! with the fragment that owns it is linear in the atom count.
      use mqc_config_types, only: input_fragment_t
      use pic_io, only: to_char

      type(input_fragment_t), intent(in) :: fragments(:)
      integer, intent(in) :: nfrag
      type(error_t), intent(out) :: error

      integer :: i, j, k, l
      integer :: atom_i, atom_j

      ! Compare each pair of fragments
      do i = 1, nfrag - 1
         do j = i + 1, nfrag
            ! Compare atoms in fragment i with atoms in fragment j
            do k = 1, size(fragments(i)%indices)
               atom_i = fragments(i)%indices(k)
               do l = 1, size(fragments(j)%indices)
                  atom_j = fragments(j)%indices(l)
                  if (atom_i == atom_j) then
                     ! Found overlapping atom
                    call error%set(ERROR_VALIDATION, "Overlapping fragments detected: fragments "//to_char(i)//" and "// &
                                    to_char(j)//" both contain atom "//to_char(atom_i)// &
                                    ". Set allow_overlapping_fragments = true to allow this.")
                     return
                  end if
               end do
            end do
         end do
      end do

   end subroutine check_fragment_overlap

end module mqc_config_adapter
