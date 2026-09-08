!! Reader for MQC JSON input files
module mqc_json_config_reader
   !! The only input reader. Parses a JSON deck into `mqc_config_t`.
   !!
   !! Three behaviours are not obvious from the schema alone:
   !!
   !!   * **A molecule gives either `xyz` or `symbols` + `geometry`, never
   !!     both.** The `xyz` path is resolved relative to the directory holding
   !!     the JSON file, not the working directory, so an input deck and its
   !!     geometries move together.
   !!   * **`is_broken` is derived, not declared.** A bond is broken when the
   !!     set of fragments containing atom i differs from the set containing
   !!     atom j. Nothing in the JSON says so; it falls out of the fragment
   !!     definitions. Comparing *sets* rather than asking whether both atoms
   !!     share some fragment is what makes it right when fragments overlap.
   !!   * **Absent keys leave defaults alone.** The `optional_*` accessors
   !!     write only when the key is present, so `mqc_config_types` holds the
   !!     defaults and nothing here restates them. Route two keys through one
   !!     local and the first one's value survives into the second.
   !!
   !! Atom indices are 0-based throughout, in the JSON and in the config, and
   !! are stored as read.
   !!
   !! Unknown keys, missing required ones, and the cross-checks between a
   !! molecule's parallel lists are `mqc_json_schema`'s job, run over the
   !! document before any of this. By the time the accessors below run, the
   !! deck is known to contain only keys this module knows about.
   use pic_io, only: to_char
   use pic_types, only: dp
   use mqc_program_limits, only: MAX_ELEMENT_SYMBOL_LEN, MAX_MBE_LEVEL, &
                                 MAX_ORBITAL_LABEL_LEN
   use mqc_geometry, only: geometry_type
   use mqc_error, only: error_t, ERROR_IO, ERROR_PARSE, ERROR_VALIDATION
   use mqc_calc_types, only: calc_type_from_string, CALC_TYPE_UNKNOWN
   use mqc_calculation_defaults, only: EFP_RESPONSE_AUTO, EFP_RESPONSE_DENSE, &
                                       EFP_RESPONSE_MATRIX_FREE
   use mqc_method_types, only: parse_method_string, method_spin_scaling, &
                               method_wants_density_fitting, method_type_to_string, &
                               method_is_casci, METHOD_TYPE_UNKNOWN, &
                               METHOD_TYPE_MP2_F12, METHOD_TYPE_CCSD_F12, &
                               METHOD_TYPE_CCSD_T_F12
   use mqc_cuest_iface, only: parse_backend_name, method_runs_on_cuest, BACKEND_CUEST, &
                              BACKEND_TERCO
   use mqc_terco_bridge, only: terco_backend_available
   use mqc_cuest_bridge, only: cuest_backend_available
   use mqc_config_types, only: mqc_config_t, input_fragment_t, bond_t
   use mqc_xyz_reader, only: read_xyz_file
   use mqc_json_schema, only: ensure_valid_json
   ! `json_integer` and `json_string` are imported here and not inside
   ! `read_neo`, where they are used: a routine-level `use json_module` in a
   ! module that already imports it makes ifx (2026.0, and 2025.3 before it)
   ! drop the default initialisation of every `json_file` local in the module
   ! -- the parser's procedure-pointer component comes up null and the first
   ! `load` jumps to address zero. gfortran is unaffected.
   use json_module, only: json_file, json_integer, json_string
   implicit none
   private

   public :: read_json_config_file  !! Parse a JSON input file into mqc_config_t
   public :: read_json_config_text  !! Parse settings from a JSON string, no molecules

contains

   subroutine read_json_config_file(filename, config, error)
      !! Read and parse a JSON format input file
      character(len=*), intent(in) :: filename
      type(mqc_config_t), intent(out) :: config
      type(error_t), intent(out) :: error

      type(json_file) :: json
      logical :: file_exists

      inquire (file=filename, exist=file_exists)
      if (.not. file_exists) then
         call error%set(ERROR_IO, "Input file not found: "//trim(filename))
         return
      end if

      call json%initialize()
      call json%load(filename=filename)
      if (json%failed()) then
         call error%set(ERROR_PARSE, "Could not parse JSON input file: "//trim(filename))
         call json%destroy()
         return
      end if

      ! Validate before reading, so a misspelled key is reported as itself
      ! rather than silently taking the default of the key it was meant to be.
      call ensure_valid_json(json, error)
      if (.not. error%has_error()) then
         call populate_config(json, directory_of(filename), config, error)
      end if
      call json%destroy()

      if (error%has_error()) call error%add_context("mqc_json_config_reader:"//trim(filename))
   end subroutine read_json_config_file

   subroutine read_json_config_text(text, config, error)
      !! Parse settings -- method, driver, keywords -- from a JSON string
      !!
      !! The same document a deck uses, less the molecules: this is how a
      !! caller that already holds a system says what to *do* with it, and it is
      !! what crosses to the workers, so they go through the validator rather
      !! than around it.
      !!
      !! A molecules block is refused rather than ignored. It would otherwise
      !! be the geometry a caller believed they were running, silently losing
      !! to the one on the handle.
      character(len=*), intent(in) :: text
      type(mqc_config_t), intent(out) :: config
      type(error_t), intent(out) :: error

      type(json_file) :: json

      call json%initialize()
      call json%deserialize(text)
      if (json%failed()) then
         call error%set(ERROR_PARSE, "Could not parse JSON settings")
         call json%destroy()
         return
      end if

      call ensure_valid_json(json, error, settings_only=.true.)
      if (.not. error%has_error()) then
         ! No base directory: an xyz path is only ever reached through a
         ! molecules block, which this mode does not have.
         call populate_config(json, "", config, error, settings_only=.true.)
      end if
      call json%destroy()

      if (error%has_error()) call error%add_context("mqc_json_config_reader:settings")
   end subroutine read_json_config_text

   subroutine populate_config(json, base_dir, config, error, settings_only)
      !! Fill every section of the config from an already-loaded document
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: base_dir  !! Directory holding the JSON, for xyz paths
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error
      logical, intent(in), optional :: settings_only  !! Stop before molecules

      character(len=:), allocatable :: text
      integer :: n_mol, imol
      logical :: found, settings, found_fukui, found_charges
      logical :: backend_named

      ! The two defaults that are not declared on the type itself.
      config%log_level = "info"
      config%fragment_breakdown = "csv"

      ! ---- schema ----------------------------------------------------------
      call require_string(json, "schema.name", config%schema_name, error)
      if (error%has_error()) return
      call require_string(json, "schema.version", config%schema_version, error)
      if (error%has_error()) return
      call optional_int(json, "schema.index_base", config%index_base)
      call optional_string(json, "schema.units", config%units)
      if (.not. allocated(config%units)) config%units = "angstrom"

      ! ---- model -----------------------------------------------------------
      call require_string(json, "model.method", text, error)
      if (error%has_error()) return
      config%method = parse_method_string(text)
      ! Refused here rather than by the factory, which can only `ERROR STOP` --
      ! survivable for a deck, and not for a driven run, where it takes the
      ! caller's interpreter down with it. The spelling still exists here, so
      ! this is also the only place that can quote it back.
      call check_method_supported(text, config%method, error)
      if (error%has_error()) return
      ! "scs-mp2" and "sos-mp2" parse to the same method type as "mp2", so the
      ! spin scaling has to be picked up here, while the spelling still exists.
      ! Left to the type alone they would run plain MP2 and report it as SCS.
      call method_spin_scaling(text, config%corr_scs, config%corr_scs_ss, config%corr_scs_os)
      ! Likewise "ri-mp2" and "df-mp2", which the method type cannot distinguish
      ! from "mp2". A later keyword can still turn it off.
      config%corr_density_fitting = method_wants_density_fitting(text)
      ! And likewise "casci", which the method type cannot distinguish from
      ! "casscf" -- both are METHOD_TYPE_MCSCF, differing only in whether the
      ! orbitals are allowed to move. `keywords.mcscf.optimize_orbitals` below
      ! can still override it.
      config%mcscf_optimize_orbitals = .not. method_is_casci(text)
      call optional_string(json, "model.basis", config%basis)
      call optional_string(json, "model.ecp", config%ecp)
      call optional_string(json, "model.aux_basis", config%aux_basis)
      call optional_logical(json, "model.cartesian", config%cartesian)
      call optional_string(json, "model.functional", config%functional)
      call optional_logical(json, "model.unrestricted", config%unrestricted)

      ! ---- properties ------------------------------------------------------
      ! The `fukui` OBJECT is what asks for the analysis; `population` only
      ! says which charges. Defaulted here so the scheme the run used is a value
      ! in the config -- and so in the report and the JSON.
      call json%info("properties.fukui", found=found_fukui)
      if (found_fukui) then
         config%fukui_population = "chelpg"
         call optional_string(json, "properties.fukui.population", &
                              config%fukui_population)
         config%fukui_guess = "neutral"
         call optional_string(json, "properties.fukui.guess", config%fukui_guess)
         ! `properties.fukui.scf` is NOT read here. It is seeded from the
         ! resolved `keywords.scf`, which this routine has not reached yet, so
         ! it is read by `read_fukui_scf` below once those values exist.
      end if
      ! Same shape as `fukui`: the OBJECT is the request, `scheme` only says
      ! how, and it is defaulted here for the same reason.
      call json%info("properties.charges", found=found_charges)
      if (found_charges) then
         config%charges_scheme = "mulliken"
         call optional_string(json, "properties.charges.scheme", &
                              config%charges_scheme)
      end if
      call optional_string(json, "properties.bonding_analysis.type", &
                           config%bonding_analysis)
      call optional_real(json, "properties.bonding_analysis.energy_threshold", &
                         config%bonding_threshold)
      call optional_logical(json, "properties.bonding_analysis.energy_decomposition", &
                            config%bonding_energy)
      call optional_logical(json, "properties.bonding_analysis.no_sharing", &
                            config%bonding_no_sharing)
      call optional_string(json, "properties.bonding_analysis.no_sharing_ci", &
                           config%bonding_no_sharing_ci)
      call optional_logical(json, "properties.bonding_analysis.restrict_localization", &
                            config%bonding_restrict_localization)

      ! ---- driver ----------------------------------------------------------
      call require_string(json, "driver", text, error)
      if (error%has_error()) return
      config%calc_type = calc_type_from_string(text)
      if (config%calc_type == CALC_TYPE_UNKNOWN) then
         call error%set(ERROR_VALIDATION, "unknown driver '"//trim(text)// &
                        "'. Accepted: energy, gradient, hessian, optimize, "// &
                        "makefp, conformers")
         return
      end if

      ! ---- system ----------------------------------------------------------
      ! Written straight into the config fields, which already hold the defaults
      ! set above. A shared local would carry the previous key's value forward
      ! whenever one is absent, since `optional_string` leaves its target alone
      ! rather than clearing it.
      call optional_string(json, "system.logger.level", config%log_level)
      call optional_real(json, "system.memory_gb", config%memory_gb)
      call optional_logical_seen(json, "system.gpu", config%gpu, config%gpu_set)
      call optional_logical(json, "system.skip_json_output", config%skip_json_output)
      call optional_logical(json, "system.unchecked_input", config%unchecked_input)
      call optional_string(json, "system.checkpoint", config%checkpoint_file)
      call optional_string(json, "system.fragment_breakdown", config%fragment_breakdown)

      ! ---- keywords --------------------------------------------------------
      call optional_int(json, "keywords.scf.maxiter", config%scf_maxiter)
      call named(json, "keywords.scf.maxiter", config%scf_maxiter_set)
      call optional_real(json, "keywords.scf.tolerance", config%scf_tolerance)
      call named(json, "keywords.scf.tolerance", config%scf_tolerance_set)
      call optional_real(json, "keywords.scf.density_tolerance", config%scf_density_tolerance)
      call optional_real(json, "keywords.scf.gradient_tolerance", config%scf_gradient_tolerance)
      call optional_real(json, "keywords.scf.level_shift", config%scf_level_shift)
      call optional_logical(json, "keywords.scf.diis", config%scf_use_diis)
      call optional_int(json, "keywords.scf.diis_size", config%scf_diis_size)
      call optional_real(json, "keywords.scf.linear_dependence_threshold", &
                         config%scf_linear_dependence)
      call named(json, "keywords.scf.density_tolerance", config%scf_density_tolerance_set)
      call optional_string(json, "keywords.scf.guess", config%scf_guess)
      call optional_string(json, "keywords.scf.accelerator", config%scf_accelerator)
      call read_scf_eri_path(json, config, error)
      if (error%has_error()) return
      call optional_string(json, "keywords.scf.convergence_metric", &
                           config%scf_convergence_metric)
      call optional_logical(json, "keywords.scf.incremental_fock", config%scf_incremental_fock)
      call optional_string(json, "keywords.guess.type", config%guess_type)
      call read_guess_steps(json, config, error)
      if (error%has_error()) return
      call optional_logical(json, "keywords.scf.allow_crap_scf", config%allow_crap_scf)

      ! After `keywords.scf`: the ion settings are *seeded* from the neutral's
      ! before the deck is consulted, so "the deck did not say" needs no
      ! sentinel -- the field still holds what it inherited.
      call read_fukui_scf(json, config)
      call optional_logical(json, "keywords.scf.density_fitting", &
                            config%scf_density_fitting)
      ! MAKEFP's own group: the response solve and the screening fit. Nothing
      ! here is an SCF setting -- that is `keywords.scf`, which it reads.
      call optional_real(json, "keywords.efp.dynamic_tolerance", &
                         config%efp_dynamic_tolerance)
      call optional_int(json, "keywords.efp.dynamic_maxiter", config%efp_dynamic_maxiter)
      call optional_logical(json, "keywords.efp.allow_crap_response", &
                            config%efp_allow_crap_response)
      call optional_int(json, "keywords.efp.response_batch", config%efp_response_batch)
      call optional_real(json, "keywords.efp.vdw_scale", config%efp_vdw_scale)
      call read_efp_response(json, config, error)
      if (error%has_error()) return
      call read_efp_dispersion(json, config, error)
      if (error%has_error()) return
      call read_neo(json, config, error)
      if (error%has_error()) return
      call optional_logical(json, "keywords.correlation.freeze_core", &
                            config%corr_freeze_core)
      call optional_int(json, "keywords.correlation.n_frozen_core", &
                        config%corr_n_frozen_core)
      ! After the method name, so an explicit keyword wins over the spelling.
      call optional_logical(json, "keywords.correlation.density_fitting", &
                            config%corr_density_fitting)

      ! A fitted SCF needs a basis to fit against, and must name its own rather
      ! than inherit a JK set that is wrong for some of what asks for fitting.
      ! `aux_basis` is unallocated exactly when `model.aux_basis` was absent.
      if (config%scf_density_fitting .and. .not. allocated(config%aux_basis)) then
         call error%set(ERROR_VALIDATION, "keywords.scf.density_fitting is on but no "// &
                        "auxiliary basis was given; set model.aux_basis")
         return
      end if

      call optional_logical(json, "keywords.correlation.scs", config%corr_scs)
      call optional_real(json, "keywords.correlation.scs_ss", config%corr_scs_ss)
      call optional_real(json, "keywords.correlation.scs_os", config%corr_scs_os)
      call optional_int(json, "keywords.cc.maxiter", config%cc_maxiter)
      call optional_real(json, "keywords.cc.tolerance", config%cc_tolerance)
      call optional_logical(json, "keywords.cc.diis", config%cc_diis)
      call optional_int(json, "keywords.cc.diis_size", config%cc_diis_size)
      call optional_logical(json, "keywords.cc.spin_adapted", config%cc_spin_adapted)
      ! Recorded as "was it named" as well as "what was it": the method name is
      ! the usual source of this and only an explicit keyword may override it,
      ! which `optional_logical` alone cannot tell from an absent key.
      call optional_logical_seen(json, "keywords.cc.triples", config%cc_triples, &
                                 config%cc_triples_set)
      call read_avas_orbitals(json, config, error)
      call read_ormas_partition(json, config, error)
      if (error%has_error()) return
      call optional_real(json, "keywords.mcscf.avas.threshold", &
                         config%mcscf_avas_threshold)
      call optional_logical(json, "keywords.mcscf.full_valence", config%mcscf_full_valence)
      call optional_int(json, "keywords.mcscf.n_active_electrons", &
                        config%mcscf_n_active_electrons)
      call optional_int(json, "keywords.mcscf.n_active_orbitals", &
                        config%mcscf_n_active_orbitals)
      call optional_int(json, "keywords.mcscf.n_inactive_orbitals", &
                        config%mcscf_n_inactive_orbitals)
      call optional_int(json, "keywords.mcscf.max_macro_iter", config%mcscf_max_macro_iter)
      call optional_real(json, "keywords.mcscf.orbital_convergence", &
                         config%mcscf_orbital_convergence)
      ! Plain `optional_logical` rather than the `_seen` variant the triples
      ! use: the field already holds what the method spelling implied, so
      ! "present wins over the name, absent loses to it" falls out without a
      ! second flag. The triples need `_seen` because their default is settled
      ! in the adapter, one hop after the spelling has been discarded.
      call optional_logical(json, "keywords.mcscf.optimize_orbitals", &
                            config%mcscf_optimize_orbitals)
      call optional_int(json, "keywords.dft.grid_level", config%dft_grid_level)
      call optional_int(json, "keywords.dft.nlc_grid_level", &
                        config%dft_nlc_grid_level)
      call optional_real(json, "keywords.dft.screening_tolerance", &
                         config%dft_screening_tolerance)
      call optional_int(json, "keywords.dft.block_size", config%dft_block_size)
      block
         ! Read into a deferred-length local, since the config field is fixed
         ! width and `optional_string` takes an allocatable.
         character(len=:), allocatable :: backend_text
         call json%get("backend", backend_text, backend_named)
         if (backend_named .and. allocated(backend_text)) config%backend = backend_text
      end block
      ! Both spellings are now in hand, so the GPU request can be settled and
      ! checked against what this build and this method can actually do.
      call resolve_gpu_request(config, backend_named, error)
      if (error%has_error()) return
      ! The continuum is switched on by the presence of the block rather than by
      ! a flag inside it: a deck that names a dielectric wants solvent, and a
      ! separate "enabled" would let the two disagree.
      call json%info("keywords.pcm", found=found)
      config%pcm_enabled = found
      block
         ! Read into a deferred-length local, since the config field is fixed
         ! width and `optional_string` takes an allocatable.
         character(len=:), allocatable :: pcm_method_text
         logical :: pcm_method_named
         call json%get("keywords.pcm.method", pcm_method_text, pcm_method_named)
         if (pcm_method_named .and. allocated(pcm_method_text)) then
            config%pcm_method = pcm_method_text
         end if
      end block
      call optional_real(json, "keywords.pcm.dielectric", config%pcm_dielectric)
      call optional_int(json, "keywords.pcm.angular_points", config%pcm_angular_points)
      call optional_real(json, "keywords.pcm.radii_scale", config%pcm_radii_scale)
      call optional_real(json, "keywords.pcm.zeta", config%pcm_zeta)
      call optional_real(json, "keywords.pcm.tolerance", config%pcm_tolerance)
      call optional_int(json, "keywords.pcm.max_iter", config%pcm_max_iter)
      call optional_int(json, "keywords.dft.radial_points", config%dft_radial_points)
      call optional_int(json, "keywords.dft.angular_points", config%dft_angular_points)

      ! Both spellings of the displacement key are accepted.
      call optional_real(json, "keywords.hessian.finite_difference_displacement", &
                         config%hessian_displacement)
      call optional_real(json, "keywords.hessian.displacement", config%hessian_displacement)
      call optional_real(json, "keywords.hessian.temperature", config%hessian_temperature)
      call optional_real(json, "keywords.hessian.pressure", config%hessian_pressure)
      call optional_real(json, "keywords.hessian.response_tolerance", &
                         config%hessian_response_tol)
      call optional_int(json, "keywords.hessian.response_max_iter", &
                        config%hessian_response_max_iter)
      call optional_int(json, "keywords.hessian.response_batch", config%hessian_response_batch)

      call optional_real(json, "keywords.aimd.dt", config%aimd_dt)
      call optional_real(json, "keywords.aimd.timestep", config%aimd_dt)
      call optional_int(json, "keywords.aimd.nsteps", config%aimd_nsteps)
      call optional_int(json, "keywords.aimd.steps", config%aimd_nsteps)
      call optional_real(json, "keywords.aimd.initial_temperature", config%aimd_initial_temperature)
      call optional_real(json, "keywords.aimd.temperature", config%aimd_initial_temperature)
      call optional_int(json, "keywords.aimd.output_frequency", config%aimd_output_frequency)
      call optional_int(json, "keywords.aimd.output_freq", config%aimd_output_frequency)

      ! Geometry optimization. Both spellings of the two that have an obvious
      ! synonym, the way the hessian and aimd blocks above accept theirs.
      call optional_int(json, "keywords.optimization.max_steps", config%opt_max_steps)
      call optional_int(json, "keywords.optimization.steps", config%opt_max_steps)
      call optional_real(json, "keywords.optimization.gradient_tolerance", &
                         config%opt_gradient_tolerance)
      call optional_real(json, "keywords.optimization.tolerance", config%opt_gradient_tolerance)
      call optional_real(json, "keywords.optimization.energy_tolerance", &
                         config%opt_energy_tolerance)
      call optional_real(json, "keywords.optimization.max_step", config%opt_max_step)
      call optional_string(json, "keywords.optimization.coordinates", config%opt_coordinates)
      call optional_string(json, "keywords.optimization.coordinate_system", config%opt_coordinates)
      call optional_string(json, "keywords.optimization.algorithm", config%opt_algorithm)
      call optional_string(json, "keywords.optimization.optimizer", config%opt_algorithm)
      call optional_int(json, "keywords.optimization.lbfgs_memory", config%opt_lbfgs_memory)
      call optional_int(json, "keywords.optimization.print_level", config%opt_print_level)
      call optional_logical(json, "keywords.optimization.trajectory", config%opt_trajectory)
      call optional_logical(json, "keywords.optimization.freeze_terms", config%opt_freeze_terms)
      call optional_logical(json, "keywords.optimization.hess_end", config%opt_hess_end)
      call optional_string(json, "keywords.optimization.hessian_update", &
                           config%opt_hessian_update)
      call optional_string(json, "keywords.optimization.target", config%opt_target)
      call optional_string(json, "keywords.optimization.endpoint", config%opt_endpoint)
      call optional_string(json, "keywords.optimization.neb_endpoints", config%opt_neb_ends)
      call optional_int(json, "keywords.optimization.images", config%opt_n_images)
      call optional_real(json, "keywords.optimization.neb_spring", config%opt_neb_spring)
      call optional_string(json, "keywords.optimization.saddle_method", &
                           config%opt_saddle_method)
      call optional_real(json, "keywords.optimization.dimer_separation", &
                         config%opt_dimer_separation)
      call optional_int(json, "keywords.optimization.dimer_max_rotations", &
                        config%opt_dimer_max_rotations)
      call optional_real(json, "keywords.optimization.dimer_rotation_tolerance", &
                         config%opt_dimer_rot_tol)
      call optional_int(json, "keywords.optimization.zero_modes", config%opt_zero_modes)
      call optional_real(json, "keywords.optimization.soft_mode_threshold", &
                         config%opt_soft_modes)
      call optional_logical(json, "keywords.optimization.connect", config%opt_connect)
      call optional_real(json, "keywords.optimization.connect_distort", &
                         config%opt_connect_distort)
      call optional_real(json, "keywords.optimization.timestep", config%opt_timestep)
      call optional_real(json, "keywords.optimization.friction", config%opt_friction)
      call optional_real(json, "keywords.optimization.friction_factor", &
                         config%opt_friction_factor)
      call optional_real(json, "keywords.optimization.friction_rising", &
                         config%opt_friction_rising)
      call read_frozen_atoms(json, config, error)
      if (error%has_error()) return
      call read_constraints(json, config, error)
      if (error%has_error()) return

      call optional_string(json, "keywords.xtb.solvent", config%solvent)
      call optional_string(json, "keywords.xtb.solvation_model", config%solvation_model)
      call optional_logical(json, "keywords.xtb.use_cds", config%use_cds)
      call optional_logical(json, "keywords.xtb.use_shift", config%use_shift)
      call optional_real(json, "keywords.xtb.dielectric", config%dielectric)
      call optional_int(json, "keywords.xtb.cpcm_nang", config%cpcm_nang)
      call optional_real(json, "keywords.xtb.cpcm_rscale", config%cpcm_rscale)

      call read_fragmentation(json, config, error)
      if (error%has_error()) return
      call read_efmo(json, config)

      ! ---- molecules -------------------------------------------------------
      settings = .false.
      if (present(settings_only)) settings = settings_only

      call json%info("molecules", found=found, n_children=n_mol)
      if (settings) then
         if (found) then
            call error%set(ERROR_VALIDATION, "settings must not define 'molecules': the "// &
                           "system comes from the handle, so a geometry given here would "// &
                           "be read, validated, and then quietly discarded")
         end if
         return
      end if
      if (.not. found .or. n_mol <= 0) then
         call error%set(ERROR_VALIDATION, "Missing or empty required key: molecules")
         return
      end if

      if (n_mol == 1) then
         ! Single-molecule mode: the fields live at the top of the config and
         ! nmol stays 0, which is what every consumer keys off.
         config%nmol = 0
         call read_molecule(json, "molecules(1)", base_dir, config%charge, &
                            config%multiplicity, config%geometry, config%nfrag, &
                            config%fragments, config%uniform_system, config%nbonds, &
                            config%nbroken, config%bonds, error)
         if (error%has_error()) return
      else
         config%nmol = n_mol
         allocate (config%molecules(n_mol))
         do imol = 1, n_mol
            call read_molecule(json, "molecules("//int_to_key(imol)//")", base_dir, &
                               config%molecules(imol)%charge, &
                               config%molecules(imol)%multiplicity, &
                               config%molecules(imol)%geometry, &
                               config%molecules(imol)%nfrag, &
                               config%molecules(imol)%fragments, &
                               config%molecules(imol)%uniform_system, &
                               config%molecules(imol)%nbonds, &
                               config%molecules(imol)%nbroken, &
                               config%molecules(imol)%bonds, error)
            if (error%has_error()) return
         end do
      end if
   end subroutine populate_config

   subroutine read_fragmentation(json, config, error)
      !! The keywords.fragmentation block, including per-level cutoffs
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      ! TODO(mqc): `level`, `highest` and `cutoff` are dead here -- the cutoff
      ! loop they belonged to now lives in `read_cutoffs`.
      integer :: level, highest
      real(dp) :: cutoff
      logical :: found

      call json%info("keywords.fragmentation", found=found)
      if (.not. found) return

      call require_string(json, "keywords.fragmentation.method", config%frag_method, error)
      if (error%has_error()) return
      call optional_int(json, "keywords.fragmentation.level", config%frag_level)
      call optional_int(json, "keywords.fragmentation.max_intersection_level", &
                        config%max_intersection_level)
      call optional_string(json, "keywords.fragmentation.counterpoise", config%counterpoise)
      call optional_string(json, "keywords.fragmentation.far_field", config%fmo_far_field)
      call optional_real(json, "keywords.fragmentation.resppc", config%fmo_resppc)
      call optional_real(json, "keywords.fragmentation.rcut", config%efmo_rcut)
      ! Unitless and a *ratio* of a distance to a van der Waals contact, so
      ! zero or negative is not "no cutoff" the way a negative `resppc` is: it
      ! would put every pair in the effective list, which is EFP with in-vacuo
      ! monomers and not the method the deck asked for. Refused rather than run.
      if (config%efmo_rcut <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "keywords.fragmentation.rcut must be "// &
                        "positive. It is a separation in units of van der Waals "// &
                        "contact, so 1.0 is touching and 2.0 (the default) is twice "// &
                        "that; a value at or below zero leaves no pair quantum "// &
                        "mechanical at all.")
         return
      end if
      call optional_int(json, "keywords.fragmentation.max_outer", config%fmo_max_outer)
      call optional_real(json, "keywords.fragmentation.outer_tolerance", config%fmo_tolerance)
      call optional_int(json, "keywords.fragmentation.scf_max_iter", config%fmo_scf_max_iter)
      call optional_real(json, "keywords.fragmentation.scf_energy_tolerance", config%fmo_scf_energy_tol)
      call optional_real(json, "keywords.fragmentation.scf_density_tolerance", config%fmo_scf_density_tol)
      call optional_string(json, "keywords.fragmentation.embedding", config%embedding)
      call optional_string(json, "keywords.fragmentation.bond_breaking", config%bond_breaking)
      call optional_real(json, "keywords.fragmentation.cap_scale", config%cap_scale)
      call optional_string(json, "keywords.fragmentation.cutoff_method", config%cutoff_method)
      call optional_string(json, "keywords.fragmentation.distance_metric", config%distance_metric)
      call optional_int(json, "keywords.fragmentation.global_groups", config%global_groups)
      call optional_int(json, "keywords.fragmentation.nodes_per_group", config%nodes_per_group)

      call read_cutoffs(json, config, error)
   end subroutine read_fragmentation

   subroutine read_efmo(json, config)
      !! The keywords.efmo block
      !!
      !! `rcut` is not here: it decides which pairs are solved quantum
      !! mechanically, which is a property of the partition, so it is read from
      !! `keywords.fragmentation` beside `resppc`.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config

      call optional_logical(json, "keywords.efmo.charge_transfer", &
                            config%efmo_charge_transfer)
   end subroutine read_efmo

   subroutine read_cutoffs(json, config, error)
      !! Per-level distance cutoffs from keywords.fragmentation.cutoffs
      !!
      !! Keys are either an n-mer name ("dimer", "trimer", ...) or the level as
      !! a decimal string ("2", "3", ...); both spellings mean the same level
      !! and the JSON may mix them. The result is an array indexed by level and
      !! sized to the highest one present, so `fragment_cutoffs(3)` is the
      !! trimer cutoff whichever spelling produced it.
      !!
      !! Both spellings are probed rather than enumerated: json-fortran offers
      !! no way to list an object's keys through the `json_file` interface.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=*), parameter :: PREFIX = "keywords.fragmentation.cutoffs."
      integer :: level, highest
      real(dp) :: cutoff, previous
      logical :: found
      logical :: present_level(MAX_MBE_LEVEL)

      call json%info("keywords.fragmentation.cutoffs", found=found)
      if (.not. found) return

      highest = 0
      present_level = .false.
      do level = 2, MAX_MBE_LEVEL
         found = .false.
         if (len(nmer_name(level)) > 0) call json%info(PREFIX//nmer_name(level), found=found)
         if (.not. found) call json%info(PREFIX//int_to_key(level), found=found)
         if (found) then
            present_level(level) = .true.
            highest = level
         end if
      end do
      if (highest < 2) return

      ! Sized to every supported level and filled with -1: downstream code reads
      ! "no cutoff at this level" off that sentinel and would take a 0 as
      ! "screen everything out".
      allocate (config%fragment_cutoffs(MAX_MBE_LEVEL))
      config%fragment_cutoffs = -1.0_dp

      previous = huge(1.0_dp)
      do level = 2, highest
         if (.not. present_level(level)) cycle
         cutoff = 0.0_dp
         if (len(nmer_name(level)) > 0) call optional_real(json, PREFIX//nmer_name(level), cutoff)
         call optional_real(json, PREFIX//int_to_key(level), cutoff)
         if (cutoff <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, PREFIX//int_to_key(level)// &
                           " must be a positive number")
            return
         end if
         ! A higher-order cutoff above a lower-order one would screen in more
         ! trimers than dimers, which cannot be what anyone meant.
         if (cutoff > previous) then
            call error%set(ERROR_VALIDATION, "keywords.fragmentation.cutoffs: the level-"// &
                           int_to_key(level)//" cutoff exceeds a lower-order one; "// &
                           "cutoffs must decrease with n-mer level")
            return
         end if
         previous = cutoff
         config%fragment_cutoffs(level) = cutoff
      end do
   end subroutine read_cutoffs

   pure function nmer_name(level) result(name)
      !! The spelled-out key for an n-mer level, or "" past the named range
      integer, intent(in) :: level
      character(len=:), allocatable :: name

      select case (level)
      case (2)
         name = "dimer"
      case (3)
         name = "trimer"
      case (4)
         name = "tetramer"
      case (5)
         name = "pentamer"
      case (6)
         name = "hexamer"
      case (7)
         name = "heptamer"
      case (8)
         name = "octamer"
      case default
         name = ""
      end select
   end function nmer_name

   subroutine check_method_supported(spelling, method_type, error)
      !! Refuse a method name nothing downstream can run
      !!
      !! Two failures, and they read the same to whoever typed the name. A
      !! spelling the parser does not know at all comes back UNKNOWN; a
      !! spelling it knows and no backend implements -- the F12 family, which
      !! has a method type and no code behind it -- parses cleanly and would
      !! otherwise stop in the factory, which can only `ERROR STOP`.
      character(len=*), intent(in) :: spelling
      integer, intent(in) :: method_type
      type(error_t), intent(inout) :: error

      character(len=*), parameter :: known = &
                                     "gfn1, gfn2, hf, dft, mp2 (also scs-, sos-, ri-, df-), ccsd, ccsd(t) "// &
                                     "(also ri-, df-), casscf, casci, mcscf, sapt0, sapt2, efp2"

      if (method_type == METHOD_TYPE_UNKNOWN) then
         call error%set(ERROR_VALIDATION, "unknown model.method '"//trim(spelling)// &
                        "'. Accepted: "//known)
         return
      end if

      select case (method_type)
      case (METHOD_TYPE_MP2_F12, METHOD_TYPE_CCSD_F12, METHOD_TYPE_CCSD_T_F12)
         call error%set(ERROR_VALIDATION, "model.method '"//trim(spelling)// &
                        "' is recognised but not implemented: there is no F12 "// &
                        "correlation factor in this program. Accepted: "//known)
         return
      case default
         ! Every other type is implemented somewhere. Whether *this build* has
         ! the backend for it is deliberately not asked here: reading a deck
         ! must mean the same thing whatever was linked. `mqc_session` asks it
         ! instead, for the callers that cannot survive an `ERROR STOP`.
      end select
   end subroutine check_method_supported

   subroutine read_efp_dispersion(json, config, error)
      !! `keywords.efp.dispersion`: "all" or "dipole", or a refusal
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: text
      character(len=:), allocatable :: lowered
      integer :: i

      call optional_string(json, "keywords.efp.dispersion", text)
      if (.not. allocated(text)) return

      lowered = trim(adjustl(text))
      do i = 1, len(lowered)
         if (lowered(i:i) >= "A" .and. lowered(i:i) <= "Z") then
            lowered(i:i) = achar(iachar(lowered(i:i)) + 32)
         end if
      end do

      select case (lowered)
      case ("all")
         config%efp_quadrupole_blocks = .true.
      case ("dipole")
         config%efp_quadrupole_blocks = .false.
      case default
         call error%set(ERROR_VALIDATION, "unknown keywords.efp.dispersion '"//trim(text)// &
                        "'. Accepted: all, dipole")
      end select
   end subroutine read_efp_dispersion

   subroutine read_scf_eri_path(json, config, error)
      !! `keywords.scf.eri_path`, checked for spelling, or a refusal
      !!
      !! Refused here rather than where the path is chosen, so a typo is
      !! reported before any integral is computed. Whether the build *has*
      !! the path named is not asked here: reading a deck must mean the same
      !! thing whatever was linked. The backend refuses that one.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: text
      character(len=:), allocatable :: lowered
      integer :: i

      call optional_string(json, "keywords.scf.eri_path", text)
      if (.not. allocated(text)) return

      lowered = trim(adjustl(text))
      do i = 1, len(lowered)
         if (lowered(i:i) >= "A" .and. lowered(i:i) <= "Z") then
            lowered(i:i) = achar(iachar(lowered(i:i)) + 32)
         end if
      end do

      select case (lowered)
      case ("rys", "rotaxis", "auto")
         config%scf_eri_path = lowered
      case default
         call error%set(ERROR_VALIDATION, "unknown keywords.scf.eri_path '"//trim(text)// &
                        "'. Accepted: rys, rotaxis, auto")
      end select
   end subroutine read_scf_eri_path

   subroutine read_efp_response(json, config, error)
      !! `keywords.efp.response`, as a code, or a refusal
      !!
      !! Refused here rather than where the choice is acted on, an SCF and a
      !! multipole expansion later, and refused rather than quietly resolved: a
      !! deck that asked for `matrix_free` and got the size rule instead would
      !! report a timing for a route it did not pick.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: text
      character(len=:), allocatable :: lowered
      integer :: i

      call optional_string(json, "keywords.efp.response", text)
      if (.not. allocated(text)) return

      lowered = trim(adjustl(text))
      do i = 1, len(lowered)
         if (lowered(i:i) >= "A" .and. lowered(i:i) <= "Z") then
            lowered(i:i) = achar(iachar(lowered(i:i)) + 32)
         end if
      end do

      select case (lowered)
      case ("auto")
         config%efp_response = EFP_RESPONSE_AUTO
      case ("dense")
         config%efp_response = EFP_RESPONSE_DENSE
      case ("matrix_free")
         config%efp_response = EFP_RESPONSE_MATRIX_FREE
      case default
         call error%set(ERROR_VALIDATION, "unknown keywords.efp.response '"//trim(text)// &
                        "'. Accepted: auto, dense, matrix_free")
      end select
   end subroutine read_efp_response

   subroutine read_neo(json, config, error)
      !! `keywords.neo`: which nuclei get orbitals of their own, and in what basis
      !!
      !! `quantum_nuclei` is either a list of 0-based atom indices or a list of
      !! element symbols; the first entry decides which, and the two are not
      !! mixed. A `neo` block without it is refused rather than read as "none":
      !! a deck that opened the block meant to quantise something.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error
      character(len=*), parameter :: path = "keywords.neo.quantum_nuclei"
      character(len=:), allocatable :: text
      integer, allocatable :: indices(:)
      logical :: found
      integer :: n, i, kind

      call json%info("keywords.neo", found=found)
      if (.not. found) return
      call optional_string(json, "keywords.neo.nuclear_basis", text)
      if (allocated(text)) then
         if (len_trim(adjustl(text)) > len(config%neo_nuclear_basis)) then
            call error%set(ERROR_VALIDATION, "keywords.neo.nuclear_basis is longer than "// &
                           int_to_key(len(config%neo_nuclear_basis))//" characters")
            return
         end if
         config%neo_nuclear_basis = trim(adjustl(text))
         deallocate (text)
      end if
      call optional_string(json, "keywords.neo.epc", text)
      if (allocated(text)) then
         select case (trim(adjustl(text)))
         case ("17-1", "epc17-1")
            config%neo_epc = "17-1"
         case ("17-2", "epc17-2")
            config%neo_epc = "17-2"
         case ("", "none")
            config%neo_epc = ""
         case default
            call error%set(ERROR_VALIDATION, "unknown keywords.neo.epc '"//trim(text)// &
                           "'. Accepted: 17-1, 17-2, none")
            return
         end select
         deallocate (text)
      end if
      call json%info(path, found=found, n_children=n)
      if (.not. found) then
         call error%set(ERROR_VALIDATION, "keywords.neo needs quantum_nuclei: a list of "// &
                        '0-based atom indices, or of element symbols such as ["H"]')
         return
      end if
      if (n < 1) then
         call error%set(ERROR_VALIDATION, "keywords.neo.quantum_nuclei is empty")
         return
      end if
      call json%info(path//"(1)", found=found, var_type=kind)
      select case (kind)
      case (json_integer)
         call json%get(path, indices, found)
         if (.not. found) then
            call error%set(ERROR_VALIDATION, "keywords.neo.quantum_nuclei could not be read "// &
                           "as a list of integers")
            return
         end if
         if (any(indices < 0)) then
            call error%set(ERROR_VALIDATION, "keywords.neo.quantum_nuclei: atom indices "// &
                           "are 0-based and cannot be negative")
            return
         end if
         config%neo_quantum_indices = indices + 1
      case (json_string)
         allocate (config%neo_quantum_symbols(n))
         do i = 1, n
            call json%get(path//"("//int_to_key(i)//")", text, found)
            if (.not. found .or. len_trim(text) == 0) then
               call error%set(ERROR_VALIDATION, "keywords.neo.quantum_nuclei: entry "// &
                              int_to_key(i)//" is not an element symbol")
               return
            end if
            config%neo_quantum_symbols(i) = trim(adjustl(text))
         end do
      case default
         call error%set(ERROR_VALIDATION, "keywords.neo.quantum_nuclei must be a list of "// &
                        "0-based atom indices or a list of element symbols")
         return
      end select
      config%neo_active = .true.
   end subroutine read_neo

   subroutine read_frozen_atoms(json, config, error)
      !! `keywords.optimization.frozen_atoms`, a flat list of atom indices
      !!
      !! Zero-based in the deck, as every atom index in this format is, and
      !! stored one-based. The conversion happens once, here.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      integer :: n_frozen, i, atom
      logical :: found

      call json%info("keywords.optimization.frozen_atoms", found=found, n_children=n_frozen)
      if (.not. found .or. n_frozen < 1) return

      allocate (config%opt_frozen_atoms(n_frozen))
      do i = 1, n_frozen
         atom = -1
         call optional_int(json, "keywords.optimization.frozen_atoms("// &
                           int_to_key(i)//")", atom)
         if (atom < 0) then
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization.frozen_atoms entry "//trim(to_char(i))// &
                           " is not an atom index. Indices are 0-based.")
            return
         end if
         config%opt_frozen_atoms(i) = atom + 1
      end do
   end subroutine read_frozen_atoms

   subroutine read_constraints(json, config, error)
      !! `keywords.optimization.constraints`, one object per held coordinate
      !!
      !! Each names a `type` and the `atoms` it is measured over, and the two
      !! have to agree: a bond takes two atoms and a torsion four. The count is
      !! the only thing that can catch the mistake -- three atoms under
      !! `"type": "torsion"` is a perfectly well-formed list of integers.
      use mqc_optimizer_types, only: constraint_from_string, constraint_atom_count
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: prefix, kind_text
      integer :: n_cons, n_atoms_here, expected, i, k, atom
      logical :: found

      call json%info("keywords.optimization.constraints", found=found, n_children=n_cons)
      if (.not. found .or. n_cons < 1) return

      allocate (config%opt_constraint_kinds(n_cons))
      allocate (config%opt_constraint_atoms(4, n_cons))
      config%opt_constraint_atoms = 0

      do i = 1, n_cons
         prefix = "keywords.optimization.constraints("//int_to_key(i)//")"

         kind_text = ""
         call optional_string(json, prefix//".type", kind_text)
         config%opt_constraint_kinds(i) = constraint_from_string(kind_text)
         expected = constraint_atom_count(config%opt_constraint_kinds(i))
         if (expected == 0) then
            call error%set(ERROR_VALIDATION, "constraint "//trim(to_char(i))// &
                           ' has type "'//trim(kind_text)//'". Use bond, angle, '// &
                           "torsion, cartesian or bond-difference.")
            return
         end if

         call json%info(prefix//".atoms", found=found, n_children=n_atoms_here)
         if (.not. found) then
            call error%set(ERROR_VALIDATION, "constraint "//trim(to_char(i))// &
                           " names no atoms")
            return
         end if
         if (n_atoms_here /= expected) then
            call error%set(ERROR_VALIDATION, "constraint "//trim(to_char(i))//' is a "'// &
                           trim(kind_text)//'", which is measured over '// &
                           trim(to_char(expected))//" atoms, but "// &
                           trim(to_char(n_atoms_here))//" were given.")
            return
         end if

         do k = 1, n_atoms_here
            atom = -1
            call optional_int(json, prefix//".atoms("//int_to_key(k)//")", atom)
            if (atom < 0) then
               call error%set(ERROR_VALIDATION, "constraint "//trim(to_char(i))// &
                              " has an atom that is not an index. Indices are 0-based.")
               return
            end if
            config%opt_constraint_atoms(k, i) = atom + 1
         end do
      end do
   end subroutine read_constraints

   subroutine read_fukui_scf(json, config)
      !! Resolve `properties.fukui.scf` against `keywords.scf`
      !!
      !! **Must run after the `keywords.scf` block.** Inheritance here is a
      !! seed followed by an overwrite: every field starts at the value that
      !! converged the neutral, and `optional_*` leaves it alone when the deck
      !! does not name the key. That is why nothing below tests for "absent"
      !! and why no field carries a sentinel -- the state "unset" is spelled by
      !! the field still holding the inherited value.
      !!
      !! `inherit_scf: false` skips the seed, leaving the type defaults for
      !! anything the object does not name.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config

      call optional_logical(json, "properties.fukui.scf.inherit_scf", &
                            config%fukui_scf%inherit_scf)

      if (config%fukui_scf%inherit_scf) then
         config%fukui_scf%max_iter = config%scf_maxiter
         config%fukui_scf%energy_tol = config%scf_tolerance
         config%fukui_scf%density_tol = config%scf_density_tolerance
         config%fukui_scf%grad_tol = config%scf_gradient_tolerance
         config%fukui_scf%linear_dependence = config%scf_linear_dependence
         config%fukui_scf%level_shift = config%scf_level_shift
         config%fukui_scf%use_diis = config%scf_use_diis
         config%fukui_scf%diis_size = config%scf_diis_size
         config%fukui_scf%incremental_fock = config%scf_incremental_fock
         config%fukui_scf%allow_crap_scf = config%allow_crap_scf
         if (allocated(config%scf_accelerator)) then
            config%fukui_scf%accelerator = config%scf_accelerator
         end if
         ! Both spellings, in the order `config_to_driver` resolves them for the
         ! neutral: `keywords.scf.guess` is the superseded one and
         ! `keywords.guess.type` the current one, so the latter wins where a
         ! deck somehow carries both.
         if (allocated(config%scf_guess)) config%fukui_scf%guess = config%scf_guess
         if (allocated(config%scf_convergence_metric)) then
            config%fukui_scf%convergence_metric = config%scf_convergence_metric
         end if
         if (allocated(config%guess_type)) config%fukui_scf%guess = config%guess_type
         ! The ladder comes across whole or not at all. A rung names a basis,
         ! and half a ladder is not a cheaper start -- it is a different one.
         if (allocated(config%guess_steps)) then
            config%fukui_scf%guess_steps = config%guess_steps
         end if
      end if

      call optional_int(json, "properties.fukui.scf.maxiter", config%fukui_scf%max_iter)
      call optional_real(json, "properties.fukui.scf.tolerance", config%fukui_scf%energy_tol)
      call optional_real(json, "properties.fukui.scf.density_tolerance", &
                         config%fukui_scf%density_tol)
      call optional_real(json, "properties.fukui.scf.gradient_tolerance", &
                         config%fukui_scf%grad_tol)
      call optional_real(json, "properties.fukui.scf.linear_dependence_threshold", &
                         config%fukui_scf%linear_dependence)
      call optional_real(json, "properties.fukui.scf.level_shift", &
                         config%fukui_scf%level_shift)
      call optional_logical(json, "properties.fukui.scf.diis", config%fukui_scf%use_diis)
      call optional_int(json, "properties.fukui.scf.diis_size", config%fukui_scf%diis_size)
      call optional_logical(json, "properties.fukui.scf.incremental_fock", &
                            config%fukui_scf%incremental_fock)
      call optional_logical(json, "properties.fukui.scf.allow_crap_scf", &
                            config%fukui_scf%allow_crap_scf)
      call read_fixed_string(json, "properties.fukui.scf.accelerator", &
                             config%fukui_scf%accelerator)
      call read_fixed_string(json, "properties.fukui.scf.convergence_metric", &
                             config%fukui_scf%convergence_metric)
      call read_fixed_string(json, "properties.fukui.scf.guess", config%fukui_scf%guess)
   end subroutine read_fukui_scf

   subroutine read_fixed_string(json, key, target)
      !! `optional_string` onto a fixed-length field, leaving it alone if absent
      !!
      !! `optional_string` takes a deferred-length allocatable, which the
      !! `character(len=32)` fields on `scf_numerics_t` are not. Reading into a
      !! local and assigning only when the key was there keeps the "absent
      !! means keep what you inherited" rule intact.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: key
      character(len=*), intent(inout) :: target
      character(len=:), allocatable :: buf

      call optional_string(json, key, buf)
      if (allocated(buf)) target = buf
   end subroutine read_fixed_string

   subroutine read_guess_steps(json, config, error)
      !! The basis-set-projection ladder, one entry per preliminary SCF
      !!
      !! Order is the order written: the first step is the cheapest basis and the
      !! last hands its density to the target one from `model.basis`, which is not
      !! listed here. Three SCFs -- STO-3G, 6-31G, cc-pVTZ -- is two steps.
      !!
      !! The schema has already refused a `subscf` that does not belong to a
      !! projection guess, and a projection guess with no steps.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: prefix
      integer :: n_steps, i
      logical :: found

      call json%info("keywords.guess.subscf.steps", found=found, n_children=n_steps)
      if (.not. found .or. n_steps < 1) return

      allocate (config%guess_steps(n_steps))
      do i = 1, n_steps
         prefix = "keywords.guess.subscf.steps("//trim(to_char(i))//")"
         call optional_string(json, prefix//".basis", config%guess_steps(i)%basis)
         if (.not. allocated(config%guess_steps(i)%basis)) then
            call error%set(ERROR_VALIDATION, "guess step "//trim(to_char(i))// &
                           " has no basis")
            return
         end if
         call optional_int(json, prefix//".maxiter", config%guess_steps(i)%maxiter)
         call optional_real(json, prefix//".tolerance", config%guess_steps(i)%tolerance)
      end do
   end subroutine read_guess_steps

   subroutine read_molecule(json, prefix, base_dir, charge, multiplicity, geom, &
                            nfrag, fragments, uniform, nbonds, nbroken, bonds, error)
      !! One entry of the molecules array
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix    !! e.g. "molecules(1)"
      character(len=*), intent(in) :: base_dir  !! For resolving a relative xyz path
      integer, intent(out) :: charge, multiplicity
      type(geometry_type), intent(inout) :: geom
      integer, intent(out) :: nfrag
      type(input_fragment_t), allocatable, intent(out) :: fragments(:)
      logical, intent(out) :: uniform
      integer, intent(out) :: nbonds, nbroken
      type(bond_t), allocatable, intent(out) :: bonds(:)
      type(error_t), intent(out) :: error

      charge = 0
      multiplicity = 1
      nfrag = 0
      nbonds = 0
      nbroken = 0

      call read_molecule_geometry(json, prefix, base_dir, geom, error)
      if (error%has_error()) return

      call require_int(json, prefix//".molecular_charge", charge, error)
      if (error%has_error()) return
      call require_int(json, prefix//".molecular_multiplicity", multiplicity, error)
      if (error%has_error()) return
      if (multiplicity <= 0) then
         call error%set(ERROR_VALIDATION, prefix//".molecular_multiplicity must be >= 1")
         return
      end if

      call read_fragments(json, prefix, base_dir, geom%natoms, nfrag, fragments, uniform, error)
      if (error%has_error()) return

      call read_connectivity(json, prefix, geom%natoms, nfrag, fragments, &
                             nbonds, nbroken, bonds, error)
   end subroutine read_molecule

   subroutine read_ormas_partition(json, config, error)
      !! `keywords.mcscf.ormas`, the subspaces and their occupation windows
      !!
      !! Absent leaves the arrays unallocated, which is what a complete active
      !! space looks like from here. The schema has already established that
      !! all three lists are present and the same length, so this only reads.
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      integer :: n_spaces, i
      logical :: found

      if (error%has_error()) return
      call json%info("keywords.mcscf.ormas.subspaces", found=found, n_children=n_spaces)
      if (.not. found .or. n_spaces <= 0) return

      allocate (config%mcscf_ormas_subspaces(n_spaces))
      allocate (config%mcscf_ormas_min_electrons(n_spaces))
      allocate (config%mcscf_ormas_max_electrons(n_spaces))

      do i = 1, n_spaces
         call require_int(json, "keywords.mcscf.ormas.subspaces("// &
                          int_to_key(i)//")", config%mcscf_ormas_subspaces(i), error)
         call require_int(json, "keywords.mcscf.ormas.min_electrons("// &
                          int_to_key(i)//")", config%mcscf_ormas_min_electrons(i), error)
         call require_int(json, "keywords.mcscf.ormas.max_electrons("// &
                          int_to_key(i)//")", config%mcscf_ormas_max_electrons(i), error)
         if (error%has_error()) return
      end do
   end subroutine read_ormas_partition

   subroutine read_avas_orbitals(json, config, error)
      !! `keywords.mcscf.avas.orbitals`, a list of atomic orbital labels
      type(json_file), intent(inout) :: json
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: label
      integer :: n_labels, i
      logical :: found

      if (error%has_error()) return
      call json%info("keywords.mcscf.avas.orbitals", found=found, n_children=n_labels)
      if (.not. found .or. n_labels <= 0) return

      allocate (character(len=MAX_ORBITAL_LABEL_LEN) :: config%mcscf_avas_orbitals(n_labels))
      do i = 1, n_labels
         call require_string(json, "keywords.mcscf.avas.orbitals("//int_to_key(i)//")", &
                             label, error)
         if (error%has_error()) return
         if (len_trim(label) > MAX_ORBITAL_LABEL_LEN) then
            call error%set(ERROR_VALIDATION, "the orbital label '"//trim(label)// &
                           "' is longer than "//int_to_key(MAX_ORBITAL_LABEL_LEN)// &
                           " characters. They read like 'N 2p' or 'Cr 3d'.")
            return
         end if
         config%mcscf_avas_orbitals(i) = trim(adjustl(label))
      end do
   end subroutine read_avas_orbitals

   subroutine read_molecule_geometry(json, prefix, base_dir, geom, error)
      !! Either an xyz file reference or inline symbols plus a flat coordinate list
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix, base_dir
      type(geometry_type), intent(inout) :: geom
      type(error_t), intent(out) :: error

      character(len=:), allocatable :: xyz_rel, symbol
      real(dp), allocatable :: flat(:)
      integer :: natoms, iatom
      logical :: has_xyz, has_symbols, found

      call json%info(prefix//".xyz", found=has_xyz)
      call json%info(prefix//".symbols", found=has_symbols, n_children=natoms)

      if (has_xyz .and. has_symbols) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": give either xyz or symbols+geometry, not both")
         return
      end if

      if (has_xyz) then
         call require_string(json, prefix//".xyz", xyz_rel, error)
         if (error%has_error()) return
         ! Relative to the input deck, so a deck and its geometries travel
         ! together regardless of where the job is launched from.
         call read_xyz_file(join_path(base_dir, xyz_rel), geom, error)
         if (error%has_error()) call error%add_context(prefix//".xyz")
         return
      end if

      if (.not. has_symbols .or. natoms <= 0) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": must give either xyz or symbols+geometry")
         return
      end if

      call json%get(prefix//".geometry", flat, found)
      if (.not. found .or. .not. allocated(flat)) then
         call error%set(ERROR_VALIDATION, prefix//".geometry is missing")
         return
      end if
      if (size(flat) /= 3*natoms) then
         call error%set(ERROR_VALIDATION, prefix//".geometry: expected "// &
                        int_to_key(3*natoms)//" values, got "//int_to_key(size(flat)))
         return
      end if

      geom%natoms = natoms
      allocate (character(len=MAX_ELEMENT_SYMBOL_LEN) :: geom%elements(natoms))
      allocate (geom%coords(3, natoms))
      do iatom = 1, natoms
         call require_string(json, prefix//".symbols("//int_to_key(iatom)//")", symbol, error)
         if (error%has_error()) return
         geom%elements(iatom) = trim(adjustl(symbol))
         geom%coords(1:3, iatom) = flat(3*(iatom - 1) + 1:3*iatom)
      end do
   end subroutine read_molecule_geometry

   subroutine read_fragments(json, prefix, base_dir, natoms, nfrag, fragments, uniform, error)
      !! The fragments array, with its parallel charge, multiplicity and potential lists
      !!
      !! `fragment_potentials` is the one of those that may be given as a **single
      !! string rather than a list**, when `uniform_system` says every fragment is
      !! the same species.
      !!
      !! Potential paths are resolved against the deck, exactly as `xyz` is.
      !!
      !! `uniform_system` is checked rather than believed: fragments claimed
      !! uniform must agree in atom count, so a deck that mixes species and says
      !! otherwise is refused rather than handed the wrong potential.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: base_dir  !! For resolving relative potential paths
      integer, intent(in) :: natoms
      integer, intent(out) :: nfrag
      type(input_fragment_t), allocatable, intent(out) :: fragments(:)
      logical, intent(out) :: uniform
      type(error_t), intent(out) :: error

      integer, allocatable :: indices(:)
      character(len=:), allocatable :: frag_path, shared, one
      integer :: ifrag, n_here, n_pot
      logical :: found, has_list

      nfrag = 0
      uniform = .false.
      call json%info(prefix//".fragments", found=found, n_children=nfrag)
      if (.not. found .or. nfrag <= 0) then
         nfrag = 0
         return
      end if

      call optional_logical(json, prefix//".uniform_system", uniform)

      ! A single string broadcasts; a list is read per fragment. Which one was written
      ! is settled by asking for the string first -- `json%get` on a list into a scalar
      ! does not find one, so this distinguishes the two without inspecting types.
      call json%get(prefix//".fragment_potentials", shared, found)
      if (.not. (found .and. allocated(shared))) then
         if (allocated(shared)) deallocate (shared)
      end if
      n_pot = 0
      call json%info(prefix//".fragment_potentials", found=has_list, n_children=n_pot)
      has_list = has_list .and. .not. allocated(shared)
      if (has_list .and. n_pot /= nfrag) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": 'fragment_potentials' has "//int_to_key(n_pot)// &
                        " entries for "//int_to_key(nfrag)//" fragments")
         return
      end if
      if (allocated(shared) .and. .not. uniform .and. nfrag > 1) then
         call error%set(ERROR_VALIDATION, prefix// &
                        ": 'fragment_potentials' is a single potential for "// &
                        int_to_key(nfrag)//" fragments, which only makes sense with "// &
                        "'uniform_system': true")
         return
      end if

      allocate (fragments(nfrag))
      do ifrag = 1, nfrag
         frag_path = prefix//".fragments("//int_to_key(ifrag)//")"

         call json%get(frag_path, indices, found)
         if (.not. found .or. .not. allocated(indices)) then
            call error%set(ERROR_PARSE, frag_path//" is not a list of atom indices")
            return
         end if
         n_here = size(indices)
         if (n_here == 0) then
            call error%set(ERROR_VALIDATION, frag_path//" is empty")
            return
         end if
         if (any(indices < 0) .or. any(indices >= natoms)) then
            call error%set(ERROR_VALIDATION, frag_path// &
                           ": atom index out of range for natoms="//int_to_key(natoms)// &
                           " (indices are 0-based)")
            return
         end if
         fragments(ifrag)%indices = indices

         ! Charges and multiplicities default per fragment, so a missing list
         ! is neutral singlets rather than an error.
         fragments(ifrag)%charge = 0
         fragments(ifrag)%multiplicity = 1
         call optional_int(json, prefix//".fragment_charges("//int_to_key(ifrag)//")", &
                           fragments(ifrag)%charge)
         call optional_int(json, prefix//".fragment_multiplicities("//int_to_key(ifrag)//")", &
                           fragments(ifrag)%multiplicity)

         ! A fragment with no potential named for it stays quantum, so an absent
         ! entry is not an error -- it is the mixed QM/EFP case.
         if (allocated(shared)) then
            fragments(ifrag)%potential = join_path(base_dir, shared)
         else if (has_list) then
            if (allocated(one)) deallocate (one)
            call json%get(prefix//".fragment_potentials("//int_to_key(ifrag)//")", &
                          one, found)
            if (found .and. allocated(one)) then
               if (len_trim(one) > 0) then
                  fragments(ifrag)%potential = join_path(base_dir, trim(one))
               end if
            end if
         end if
      end do

      ! The uniformity claim, checked. Atom count is the cheap necessary
      ! condition for a potential's atoms to be superposable onto a fragment's;
      ! the loader still checks composition when it reads the file, which is
      ! where the elements become known.
      if (uniform) then
         do ifrag = 2, nfrag
            if (size(fragments(ifrag)%indices) /= size(fragments(1)%indices)) then
               call error%set(ERROR_VALIDATION, prefix// &
                              ": 'uniform_system' is true but fragment "// &
                              int_to_key(ifrag)//" has "// &
                              int_to_key(size(fragments(ifrag)%indices))// &
                              " atoms where the first has "// &
                              int_to_key(size(fragments(1)%indices)))
               return
            end if
         end do
      end if
   end subroutine read_fragments

   subroutine read_connectivity(json, prefix, natoms, nfrag, fragments, &
                                nbonds, nbroken, bonds, error)
      !! The connectivity list, deriving which bonds fragmentation breaks
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: prefix
      integer, intent(in) :: natoms, nfrag
      type(input_fragment_t), allocatable, intent(in) :: fragments(:)
      integer, intent(out) :: nbonds, nbroken
      type(bond_t), allocatable, intent(out) :: bonds(:)
      type(error_t), intent(out) :: error

      integer, allocatable :: triple(:)
      character(len=:), allocatable :: bond_path
      integer :: ibond
      logical :: found

      nbonds = 0
      nbroken = 0
      call json%info(prefix//".connectivity", found=found, n_children=nbonds)
      if (.not. found .or. nbonds <= 0) then
         nbonds = 0
         return
      end if

      allocate (bonds(nbonds))
      do ibond = 1, nbonds
         bond_path = prefix//".connectivity("//int_to_key(ibond)//")"
         call json%get(bond_path, triple, found)
         if (.not. found .or. .not. allocated(triple)) then
            call error%set(ERROR_PARSE, bond_path//" is not a list")
            return
         end if
         if (size(triple) /= 3) then
            call error%set(ERROR_VALIDATION, bond_path//" must be [i, j, order]")
            return
         end if
         if (any(triple(1:2) < 0) .or. any(triple(1:2) >= natoms)) then
            call error%set(ERROR_VALIDATION, bond_path// &
                           ": atom index out of range for natoms="//int_to_key(natoms))
            return
         end if
         if (triple(3) <= 0) then
            call error%set(ERROR_VALIDATION, bond_path//": bond order must be positive")
            return
         end if

         bonds(ibond)%atom_i = triple(1)
         bonds(ibond)%atom_j = triple(2)
         bonds(ibond)%order = triple(3)
         bonds(ibond)%is_broken = bond_is_broken(triple(1), triple(2), nfrag, fragments)
         if (bonds(ibond)%is_broken) nbroken = nbroken + 1
      end do
   end subroutine read_connectivity

   pure function bond_is_broken(atom_i, atom_j, nfrag, fragments) result(broken)
      !! Whether fragmentation severs the bond between two atoms
      !!
      !! Broken means the two atoms do not belong to the same set of fragments.
      !! Comparing *sets* rather than asking "are they both in some fragment"
      !! is what makes this correct for overlapping fragments, where an atom
      !! can legitimately be in more than one.
      integer, intent(in) :: atom_i, atom_j
      integer, intent(in) :: nfrag
      type(input_fragment_t), allocatable, intent(in) :: fragments(:)
      logical :: broken

      integer :: ifrag
      logical :: in_i, in_j

      broken = .false.
      if (nfrag <= 0 .or. .not. allocated(fragments)) return

      do ifrag = 1, nfrag
         in_i = any(fragments(ifrag)%indices == atom_i)
         in_j = any(fragments(ifrag)%indices == atom_j)
         if (in_i .neqv. in_j) then
            broken = .true.
            return
         end if
      end do
   end function bond_is_broken

   ! ==========================================================================
   !  Accessors
   !
   !  json-fortran reports a missing key through `found` rather than failing,
   !  so "required" and "optional" differ only in what happens when it is
   !  false. The optional forms leave the target at whatever default the config
   !  type declared, which is how the defaults stay in one place.
   ! ==========================================================================

   subroutine require_string(json, path, value, error)
      !! Fetch a string, or record a missing-key error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(out) :: value
      type(error_t), intent(out) :: error

      logical :: found

      call json%get(path, value, found)
      if (.not. found .or. .not. allocated(value)) then
         call error%set(ERROR_VALIDATION, "Missing required key: "//path)
      end if
   end subroutine require_string

   subroutine require_int(json, path, value, error)
      !! Fetch an integer, or record a missing-key error
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(out) :: value
      type(error_t), intent(out) :: error

      logical :: found

      call json%get(path, value, found)
      if (.not. found) call error%set(ERROR_VALIDATION, "Missing required key: "//path)
   end subroutine require_int

   subroutine optional_string(json, path, value)
      !! Fetch a string if present, leaving `value` untouched otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(inout) :: value

      character(len=:), allocatable :: text
      logical :: found

      call json%get(path, text, found)
      if (found .and. allocated(text)) value = text
   end subroutine optional_string

   subroutine optional_int(json, path, value)
      !! Fetch an integer if present, leaving `value` at its default otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(inout) :: value

      integer :: found_value
      logical :: found

      call json%get(path, found_value, found)
      if (found) value = found_value
   end subroutine optional_int

   subroutine named(json, path, was_named)
      !! Whether the deck mentioned a key at all
      !!
      !! `optional_real` leaves a default in place when a key is absent, which
      !! loses the difference between "the user did not say" and "the user asked
      !! for exactly the default". MAKEFP needs it: it converges its SCF harder
      !! than the shared default, because the multipoles and the response are
      !! taken from that density, unless a deck actually asks otherwise.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      logical, intent(out) :: was_named

      real(dp) :: ignored

      call json%get(path, ignored, was_named)
   end subroutine named

   subroutine optional_real(json, path, value)
      !! Fetch a real if present, leaving `value` at its default otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      real(dp), intent(inout) :: value

      real(dp) :: found_value
      logical :: found

      call json%get(path, found_value, found)
      if (found) value = found_value
   end subroutine optional_real

   subroutine resolve_gpu_request(config, backend_named, error)
      !! Settle `system.gpu` against `backend`, then refuse a request nothing can honour
      !!
      !! Three things happen here, all of them before a calculation is set up
      !! rather than once per fragment inside one:
      !!
      !!   * `system.gpu` resolves into `backend`, so everything downstream reads
      !!     one field however the deck spelled the question. Naming both and
      !!     disagreeing -- `"gpu": true` beside `"backend": "cpu"` -- is refused
      !!     rather than resolved by precedence.
      !!   * asking for cuEST on a build without it is refused. This is asked of
      !!     the bridge that linked, not of a preprocessor symbol, so the answer
      !!     describes the binary.
      !!   * asking for cuEST for a method it has no implementation of is refused.
      !!     Substituting the CPU silently would report a provenance and a set of
      !!     timings that were not true.
      type(mqc_config_t), intent(inout) :: config
      logical, intent(in) :: backend_named  !! Whether the deck had a root `backend` key
      type(error_t), intent(out) :: error

      integer :: backend_kind

      if (config%gpu_set) then
         if (backend_named) then
            call parse_backend_name(config%backend, backend_kind, error)
            if (error%has_error()) return
            if ((backend_kind == BACKEND_CUEST .or. &
                 backend_kind == BACKEND_TERCO) .neqv. config%gpu) then
               call error%set(ERROR_VALIDATION, "system.gpu and backend disagree: "// &
                              "system.gpu is "//trim(merge("true ", "false", config%gpu))// &
                              " but backend is '"//trim(config%backend)// &
                              "'. Name one of them, not both.")
               return
            end if
         end if
         ! Only when the deck did NOT name a backend. Overwriting a named one
         ! here would take `"gpu": true` beside `"backend": "terco"` and run
         ! cuEST, reporting a backend the deck never asked for.
         if (.not. backend_named) then
            if (config%gpu) then
               config%backend = "cuest"
            else
               config%backend = "libcint"
            end if
         end if
      end if

      call parse_backend_name(config%backend, backend_kind, error)
      if (error%has_error()) return

      if (backend_kind == BACKEND_TERCO) then
         if (.not. terco_backend_available()) then
            call error%set(ERROR_VALIDATION, "backend 'terco' was asked for, but this "// &
                           "build has no terco backend. Rebuild with CMake and "// &
                           "-DMQC_ENABLE_TERCO=ON -DTERCO_ROOT=<terco checkout>, or "// &
                           "drop the request to run on the CPU.")
         end if
         return
      end if

      if (backend_kind /= BACKEND_CUEST) return

      if (.not. cuest_backend_available()) then
         call error%set(ERROR_VALIDATION, gpu_asked_by(config)//" but this build has "// &
                        "no cuEST backend. Rebuild with CMake and "// &
                        "-DMQC_ENABLE_CUEST=ON -DCUEST_ROOT=<cuest archive>, or drop "// &
                        "the request to run on the CPU.")
         return
      end if

      if (.not. method_runs_on_cuest(config%method)) then
         call error%set(ERROR_VALIDATION, gpu_asked_by(config)//" but '"// &
                        trim(method_type_to_string(config%method))//"' has no cuEST "// &
                        "implementation; only Hartree-Fock and DFT are offloaded. "// &
                        "Run it on the CPU, or pick a method the GPU has.")
         return
      end if
   end subroutine resolve_gpu_request

   pure function gpu_asked_by(config) result(text)
      !! Which key asked for the GPU, so the refusal names what the deck wrote
      type(mqc_config_t), intent(in) :: config
      character(len=:), allocatable :: text

      if (config%gpu_set) then
         text = "system.gpu is true"
      else
         text = "backend '"//trim(config%backend)//"' was asked for"
      end if
   end function gpu_asked_by

   subroutine optional_logical_seen(json, path, value, seen)
      !! Fetch a boolean, and say whether the deck named it
      !!
      !! Needed wherever a default comes from somewhere other than the field's
      !! initialiser -- the method name, say -- because "absent" and "false" have
      !! to be told apart before the name can be allowed to lose to a keyword.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      logical, intent(inout) :: value
      logical, intent(out) :: seen

      logical :: found_value

      call json%get(path, found_value, seen)
      if (seen) value = found_value
   end subroutine optional_logical_seen

   subroutine optional_logical(json, path, value)
      !! Fetch a boolean if present, leaving `value` at its default otherwise
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      logical, intent(inout) :: value

      logical :: found_value
      logical :: found

      call json%get(path, found_value, found)
      if (found) value = found_value
   end subroutine optional_logical

   ! ==========================================================================
   !  Small string utilities
   ! ==========================================================================

   pure function int_to_key(value) result(key)
      !! An integer as the string json-fortran's paths and keys want
      integer, intent(in) :: value
      character(len=:), allocatable :: key
      character(len=16) :: buffer

      write (buffer, "(I0)") value
      key = trim(buffer)
   end function int_to_key

   pure function directory_of(path) result(directory)
      !! The directory part of a path, or "." when there is none
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: directory

      integer :: slash

      slash = index(path, "/", back=.true.)
      if (slash > 0) then
         directory = path(1:slash - 1)
      else
         directory = "."
      end if
   end function directory_of

   pure function join_path(directory, relative) result(joined)
      !! Join a directory and a relative path
      !!
      !! An absolute `relative` wins, so a deck can point at a geometry
      !! somewhere else entirely.
      character(len=*), intent(in) :: directory, relative
      character(len=:), allocatable :: joined

      if (len_trim(relative) > 0) then
         if (relative(1:1) == "/") then
            joined = trim(relative)
            return
         end if
      end if
      joined = trim(directory)//"/"//trim(relative)
   end function join_path

end module mqc_json_config_reader
