module test_mqc_json_reader
   !! Behavioural tests for the JSON input reader.
   !!
   !! Ported from the `.mqc` parser's suite when that format was retired: the
   !! assertions are the same requirements, restated against the format that
   !! replaced it. They are what the reader is actually *for*, as distinct from
   !! `test_mqc_json_config`, which checks it against real decks.
   !!
   !! Each test writes a scratch deck and reads it back. The deck is assembled
   !! by `write_deck` from a handful of optional blocks, so a test says what it
   !! varies and nothing else -- the `.mqc` originals repeated forty lines of
   !! file-writing per case, which buried the one line that mattered.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_types, only: mqc_config_t
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, &
                               METHOD_TYPE_DFT, METHOD_TYPE_MP2, METHOD_TYPE_CCSD_T, &
                               METHOD_TYPE_MCSCF
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE, DEFAULT_RESPONSE_TOL, &
                                       DEFAULT_RESPONSE_MAX_ITER, DEFAULT_SCF_CONV, &
                                       DEFAULT_SCF_DENSITY_CONV, DEFAULT_VDW_SCALE, &
                                       DEFAULT_DYNAMIC_TOL, DEFAULT_DYNAMIC_MAXITER, &
                                       EFP_RESPONSE_AUTO, EFP_RESPONSE_DENSE, &
                                       EFP_RESPONSE_MATRIX_FREE
   use mqc_error, only: error_t
   use mqc_cuest_iface, only: parse_backend_name, BACKEND_AUTO, BACKEND_CUEST, &
                              BACKEND_CZT, method_runs_on_cuest
   use mqc_cuest_bridge, only: cuest_backend_available
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_json_reader_tests
   public :: remove_deck  !! Tidy the scratch file away

   character(len=*), parameter :: DECK = "test_json_reader_scratch.json"

contains

   subroutine collect_mqc_json_reader_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("minimal_deck", test_minimal), &
                  new_unittest("scf_eri_path_is_read_and_a_misspelling_refused", &
                               test_scf_eri_path), &
                  new_unittest("fragments", test_fragments), &
                  new_unittest("connectivity_marks_broken_bonds", test_connectivity), &
                  new_unittest("no_fragments", test_no_fragments), &
                  new_unittest("xtb_method_spelling", test_method_xtb), &
                  new_unittest("logger_level", test_log_level), &
                  new_unittest("memory_budget", test_memory_gb), &
                  new_unittest("efp_dispersion_blocks", test_efp_dispersion), &
                  new_unittest("hessian_settings", test_hessian), &
                  new_unittest("allow_crap_scf", test_allow_crap_scf), &
                  new_unittest("scf_tolerances_record_being_named", test_scf_tolerances), &
                  new_unittest("efp_keywords", test_efp_keywords), &
                  new_unittest("neo_keywords", test_neo_keywords), &
                  new_unittest("hessian_defaults", test_hessian_defaults), &
                  new_unittest("aimd_settings", test_aimd), &
                  new_unittest("fragmentation_settings", test_fragmentation), &
                  new_unittest("bond_breaking_defaults", test_bond_breaking_defaults), &
                  new_unittest("fmo_scf_keywords", test_fmo_scf_keywords), &
                  new_unittest("df_without_aux_fails", test_df_without_aux), &
                  new_unittest("fragmentation_cutoffs", test_cutoffs), &
                  new_unittest("efmo_keywords", test_efmo_keywords), &
                  new_unittest("efmo_rcut_must_be_positive", test_efmo_rcut_refused), &
                  new_unittest("cutoffs_must_decrease", test_cutoffs_monotonic), &
                  new_unittest("global_groups", test_global_groups), &
                  new_unittest("nodes_per_group", test_nodes_per_group), &
                  new_unittest("xtb_solvation", test_solvation), &
                  new_unittest("inline_symbols_and_geometry", test_inline_geometry), &
                  new_unittest("error_missing_schema", test_missing_schema), &
                  new_unittest("error_missing_molecules", test_missing_molecules), &
                  new_unittest("cc_keywords", test_cc_keywords), &
                  new_unittest("cc_spin_adapted_keyword", test_cc_spin_adapted), &
                  new_unittest("mcscf_keywords", test_mcscf_keywords), &
                  new_unittest("casci_spelling_fixes_the_orbitals", test_casci_spelling), &
                  new_unittest("backend_keyword", test_backend_keyword), &
                  new_unittest("system_gpu_keyword", test_gpu_keyword), &
                  new_unittest("system_gpu_conflicts_with_backend", test_gpu_conflict), &
                  new_unittest("cuest_method_allow_list", test_cuest_methods), &
                  new_unittest("pcm_keywords", test_pcm_keywords), &
                  new_unittest("dft_keywords", test_dft_keywords), &
                  new_unittest("fragment_potentials", test_fragment_potentials), &
                  new_unittest("uniform_system_broadcast", test_uniform_system), &
                  new_unittest("uniform_system_is_checked", test_uniform_rejected), &
                  new_unittest("bonding_analysis_property", test_bonding_analysis), &
                  new_unittest("charges_property", test_charges_property), &
                  new_unittest("avas_orbital_labels", test_avas_keywords), &
                  new_unittest("ormas_partition", test_ormas_keywords), &
                  new_unittest("full_valence_space", test_full_valence), &
                  new_unittest("optimization_hess_end", test_hess_end), &
                  new_unittest("optimization_target", test_opt_target), &
                  new_unittest("scf_diis_controls", test_diis_keywords), &
                  new_unittest("incremental_fock keyword", test_incremental_fock_keyword), &
                  new_unittest("fukui_scf_inherits_keywords_scf", test_fukui_scf_inherit), &
                  new_unittest("fukui_scf_overrides_per_key", test_fukui_scf_override), &
                  new_unittest("fukui_scf_inherit_false_drops_the_seed", test_fukui_scf_no_inherit), &
                  new_unittest("fukui_scf_inherits_the_current_guess_spelling", &
                               test_fukui_scf_guess_spelling), &
                  new_unittest("scf_maxiter_named_flag", test_maxiter_named), &
                  new_unittest("scf_convergence_metric", test_convergence_metric), &
                  new_unittest("error_malformed_json", test_malformed) &
                  ]
   end subroutine collect_mqc_json_reader_tests

   ! ==========================================================================
   !  Deck construction
   ! ==========================================================================

   subroutine write_deck(model, driver, keywords, system, molecule, root)
      !! Write a scratch deck, omitting any block given as an empty string
      !!
      !! Every argument is a JSON fragment without its trailing comma; the
      !! commas are placed here so a caller never has to think about them.
      character(len=*), intent(in) :: model     !! Body of "model"
      character(len=*), intent(in) :: driver    !! Value of "driver"
      character(len=*), intent(in) :: keywords  !! Body of "keywords", or ""
      character(len=*), intent(in) :: system    !! Body of "system", or ""
      character(len=*), intent(in) :: molecule  !! Body of molecules[0]
      character(len=*), intent(in), optional :: root
         !! Extra root-level keys, for the few that live there rather than in a
         !! block -- `backend` is the one this exists for.

      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "schema": {"name": "mqc-frag", "version": "1.0"},'
      if (present(root)) then
         if (len_trim(root) > 0) write (unit, "(A)") "  "//root//","
      end if
      write (unit, "(A)") '  "model": {'//model//'},'
      write (unit, "(A)") '  "driver": "'//driver//'",'
      if (len_trim(keywords) > 0) write (unit, "(A)") '  "keywords": {'//keywords//'},'
      if (len_trim(system) > 0) write (unit, "(A)") '  "system": {'//system//'},'
      write (unit, "(A)") '  "molecules": [{'//molecule//'}]'
      write (unit, "(A)") "}"
      close (unit)
   end subroutine write_deck

   pure function two_atoms() result(text)
      !! A two-atom hydrogen molecule given inline
      character(len=:), allocatable :: text
      text = '"symbols": ["H", "H"], "geometry": [0.0, 0.0, 0.0, 0.7, 0.0, 0.0], '// &
             '"molecular_charge": 0, "molecular_multiplicity": 1'
   end function two_atoms

   subroutine test_scf_eri_path(error)
      !! `keywords.scf.eri_path`: absent leaves it unset, a spelling is kept
      !! lowercase, and anything else is refused at read time
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. allocated(config%scf_eri_path), &
                 "a deck that says nothing must leave eri_path to the backend")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"scf": {"eri_path": "RotAxis"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%scf_eri_path), "scf.eri_path was not read")
      if (allocated(error)) return
      call check(error, config%scf_eri_path == "rotaxis", &
                 "scf.eri_path must be stored lowercase, got '"//config%scf_eri_path//"'")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"scf": {"eri_path": "rotated"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled eri_path must be refused when the deck is read")
      if (allocated(error)) return
      call check(error, index(parse_error%get_message(), "eri_path") > 0, &
                 "the refusal must name the key")
   end subroutine test_scf_eri_path

   subroutine read_deck(config, parse_error)
      !! Read the scratch deck back
      type(mqc_config_t), intent(out) :: config
      type(error_t), intent(out) :: parse_error

      call read_json_config_file(DECK, config, parse_error)
   end subroutine read_deck

   ! ==========================================================================
   !  Tests
   ! ==========================================================================

   subroutine test_minimal(error)
      !! The smallest deck the reader accepts, and every field it implies
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), &
                 "minimal deck should parse: "//parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%schema_name, "mqc-frag")
      if (allocated(error)) return
      call check(error, config%schema_version, "1.0")
      if (allocated(error)) return
      call check(error, config%index_base, 0, "index_base defaults to 0")
      if (allocated(error)) return
      call check(error, config%units, "angstrom", "units default to angstrom")
      if (allocated(error)) return
      call check(error, config%method, METHOD_TYPE_GFN2)
      if (allocated(error)) return
      call check(error, config%calc_type, CALC_TYPE_ENERGY)
      if (allocated(error)) return
      call check(error, config%charge, 0)
      if (allocated(error)) return
      call check(error, config%multiplicity, 1)
      if (allocated(error)) return
      call check(error, config%nmol, 0, "one molecule uses single-molecule mode")
      if (allocated(error)) return
      call check(error, config%geometry%natoms, 2)
      if (allocated(error)) return
      call check(error, trim(config%geometry%elements(1)), "H")
      if (allocated(error)) return
      call check(error, close_enough(config%geometry%coords(1, 2), 0.7_dp), &
                 "second atom sits at x = 0.7")
   end subroutine test_minimal

   subroutine test_fragments(error)
      !! Fragment definitions and their per-fragment charge and multiplicity
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Gradient", "", "", &
                      '"symbols": ["H", "H", "O", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0], '// &
                      '"molecular_charge": 1, "molecular_multiplicity": 2, '// &
                      '"fragments": [[0, 1], [2, 3]], '// &
                      '"fragment_charges": [0, 1], '// &
                      '"fragment_multiplicities": [1, 2]')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%method, METHOD_TYPE_GFN1)
      if (allocated(error)) return
      call check(error, config%calc_type, CALC_TYPE_GRADIENT)
      if (allocated(error)) return
      call check(error, config%charge, 1)
      if (allocated(error)) return
      call check(error, config%multiplicity, 2)
      if (allocated(error)) return
      call check(error, config%nfrag, 2)
      if (allocated(error)) return
      call check(error, size(config%fragments(1)%indices), 2)
      if (allocated(error)) return
      call check(error, config%fragments(1)%indices(1), 0, "indices stay 0-based")
      if (allocated(error)) return
      call check(error, config%fragments(1)%charge, 0)
      if (allocated(error)) return
      call check(error, config%fragments(2)%charge, 1)
      if (allocated(error)) return
      call check(error, config%fragments(2)%multiplicity, 2)
   end subroutine test_fragments

   subroutine test_connectivity(error)
      !! A bond between two fragments is broken; one inside a fragment is not
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", &
                      '"symbols": ["H", "H", "H", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1, '// &
                      '"fragments": [[0, 1], [2, 3]], '// &
                      '"connectivity": [[0, 1, 1], [1, 2, 1]]')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nbonds, 2)
      if (allocated(error)) return
      ! 0-1 lies inside fragment 1; 1-2 crosses from fragment 1 to fragment 2.
      call check(error,.not. config%bonds(1)%is_broken, &
                 "a bond inside one fragment is not broken")
      if (allocated(error)) return
      call check(error, config%bonds(2)%is_broken, &
                 "a bond between two fragments is broken")
      if (allocated(error)) return
      call check(error, config%nbroken, 1)
      if (allocated(error)) return
      call check(error, config%bonds(2)%atom_i, 1)
      if (allocated(error)) return
      call check(error, config%bonds(2)%atom_j, 2)
      if (allocated(error)) return
      call check(error, config%bonds(2)%order, 1)
   end subroutine test_connectivity

   subroutine test_no_fragments(error)
      !! A deck without a fragments key is unfragmented, not an error
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nfrag, 0)
      if (allocated(error)) return
      call check(error, config%nbonds, 0)
   end subroutine test_no_fragments

   subroutine test_method_xtb(error)
      !! "XTB-GFN1" and "gfn1" mean the same method
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%method, METHOD_TYPE_GFN1)
      if (allocated(error)) return

      call write_deck('"method": "gfn1"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%method, METHOD_TYPE_GFN1, &
                 "the bare spelling should reach the same method")
   end subroutine test_method_xtb

   subroutine test_efp_dispersion(error)
      !! keywords.efp.dispersion switches the quadrupole blocks, defaults to all
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf"', "MakeFP", '"efp": {"dispersion": "dipole"}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%efp_quadrupole_blocks, "dipole should switch the blocks off")
      if (allocated(error)) return

      call write_deck('"method": "hf"', "MakeFP", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%efp_quadrupole_blocks, "the blocks should default to on")
      if (allocated(error)) return

      call write_deck('"method": "hf"', "MakeFP", '"efp": {"dispersion": "octupole"}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "an unknown dispersion word must be refused")
   end subroutine test_efp_dispersion

   subroutine test_memory_gb(error)
      !! system.memory_gb reaches memory_gb, and stays negative when absent
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", &
                      '"memory_gb": 120.5', two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, abs(config%memory_gb - 120.5_dp) < 1.0e-12_dp, &
                 "memory_gb should carry the deck's figure")
      if (allocated(error)) return

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%memory_gb < 0.0_dp, &
                 "memory_gb should stay negative, meaning the machine decides")
   end subroutine test_memory_gb

   subroutine test_log_level(error)
      !! system.logger.level reaches log_level, and defaults when absent
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", &
                      '"logger": {"level": "verbose"}', two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%log_level), "log_level should be set")
      if (allocated(error)) return
      call check(error, config%log_level, "verbose")
      if (allocated(error)) return

      ! Absent, the default stands -- and does not pick up a neighbouring key.
      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%log_level, "info", "log_level defaults to info")
      if (allocated(error)) return
      call check(error, config%fragment_breakdown, "csv", &
                 "fragment_breakdown defaults to csv")
   end subroutine test_log_level

   subroutine test_allow_crap_scf(error)
      !! The escape hatch is off unless the deck says otherwise
      !!
      !! Both directions are checked. The default matters more than the set
      !! value: it decides whether a non-converged SCF stops the run or is
      !! folded into the total, and a default that drifted to true would put
      !! this program back to reporting numbers built from SCFs that never
      !! settled -- which is the failure it exists to prevent.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"scf": {"maxiter": 40}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%allow_crap_scf, &
                 "allow_crap_scf must default to false when the deck is silent")
      if (allocated(error)) return

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"scf": {"maxiter": 40, "allow_crap_scf": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%allow_crap_scf, "allow_crap_scf was not read")
   end subroutine test_allow_crap_scf

   subroutine test_scf_tolerances(error)
      !! Both SCF thresholds, and whether the deck named them
      !!
      !! The values alone are not enough here, which is why the flags exist. A
      !! caller with a stricter default of its own -- MAKEFP converges to
      !! 1e-10/1e-8 because the multipoles and the response come off that
      !! density -- has to be able to tell "the deck asked for 1e-6" from "the
      !! deck said nothing and 1e-6 is what the field was initialised to". Both
      !! arrive as 1e-6. If these flags ever come back true for a silent deck,
      !! MAKEFP quietly loses three digits of density; if they stay false for a
      !! deck that named a tolerance, the user is ignored, which is the bug this
      !! replaced -- a deck asking for 1e-6 watched the SCF grind past 37 steps.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 40}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%scf_tolerance_set, &
                 "a deck that never mentions scf.tolerance must not look like it did")
      if (allocated(error)) return
      call check(error,.not. config%scf_density_tolerance_set, &
                 "a deck that never mentions scf.density_tolerance must not look like it did")
      if (allocated(error)) return
      call check(error, close_enough(config%scf_tolerance, DEFAULT_SCF_CONV), &
                 "the silent deck must still carry the shared default")
      if (allocated(error)) return
      call check(error, close_enough(config%scf_density_tolerance, DEFAULT_SCF_DENSITY_CONV), &
                 "the silent deck must still carry the shared density default")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"tolerance": 1.0e-8, "density_tolerance": 1.0e-7}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%scf_tolerance, 1.0e-8_dp), &
                 "scf.tolerance was not read")
      if (allocated(error)) return
      call check(error, close_enough(config%scf_density_tolerance, 1.0e-7_dp), &
                 "scf.density_tolerance was not read")
      if (allocated(error)) return
      call check(error, config%scf_tolerance_set, &
                 "scf.tolerance was read but not recorded as named")
      if (allocated(error)) return
      call check(error, config%scf_density_tolerance_set, &
                 "scf.density_tolerance was read but not recorded as named")
      if (allocated(error)) return

      ! Naming the default is not the same request as saying nothing, and this
      ! is the case that a bare `optional_real` cannot tell apart.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"tolerance": 1.0e-6}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%scf_tolerance_set, &
                 "asking for exactly the default is still asking")
      if (allocated(error)) return
      call check(error,.not. config%scf_density_tolerance_set, &
                 "one key named must not set the other's flag")
   end subroutine test_scf_tolerances

   subroutine test_efp_keywords(error)
      !! `keywords.efp`: the silent deck, the deck that names all four, the typo
      !!
      !! The first case is the one that matters. Everything in this group is a
      !! number the MAKEFP path used to hold as a literal, and the whole of the
      !! change is that a deck can now name it -- so a deck that names none of
      !! them has to arrive carrying exactly what the solver would have used on
      !! its own. That is why these defaults are the constants themselves and not
      !! copies of their values: a copy stays right until someone moves the
      !! original, and then a silent deck starts overriding it with the old
      !! number, which is invisible from anywhere.
      !!
      !! The third case pins a refusal. `response` picks between two routes to the
      !! same equations, and a deck asking for one is usually about to time it; a
      !! misspelling that fell back to "auto" would produce a number the user
      !! would read as belonging to the route they asked for.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "6-31g"', "MakeFP", &
                      '"scf": {"maxiter": 40}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%efp_dynamic_tolerance, DEFAULT_DYNAMIC_TOL), &
                 "a deck with no efp block must carry the solver's own tolerance")
      if (allocated(error)) return
      call check(error, config%efp_dynamic_maxiter == DEFAULT_DYNAMIC_MAXITER, &
                 "a deck with no efp block must carry the solver's own iteration cap")
      if (allocated(error)) return
      call check(error, config%efp_response == EFP_RESPONSE_AUTO, &
                 "a deck with no efp block must leave the route to the size rule")
      if (allocated(error)) return
      call check(error, close_enough(config%efp_vdw_scale, DEFAULT_VDW_SCALE), &
                 "a deck with no efp block must carry the screening grid's own scale")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "MakeFP", &
                      '"efp": {"dynamic_tolerance": 5.0e-5, "dynamic_maxiter": 40, '// &
                      '"response": "matrix_free", "vdw_scale": 0.8}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%efp_dynamic_tolerance, 5.0e-5_dp), &
                 "efp.dynamic_tolerance was not read")
      if (allocated(error)) return
      call check(error, config%efp_dynamic_maxiter == 40, &
                 "efp.dynamic_maxiter was not read")
      if (allocated(error)) return
      call check(error, config%efp_response == EFP_RESPONSE_MATRIX_FREE, &
                 "efp.response was not read")
      if (allocated(error)) return
      call check(error, close_enough(config%efp_vdw_scale, 0.8_dp), &
                 "efp.vdw_scale was not read")
      if (allocated(error)) return

      ! The other two spellings, since only one of the three can be the default.
      call write_deck('"method": "hf", "basis": "6-31g"', "MakeFP", &
                      '"efp": {"response": "Dense"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%efp_response == EFP_RESPONSE_DENSE, &
                 "efp.response 'Dense' was not read, or was read case-sensitively")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "MakeFP", &
                      '"efp": {"response": "matrix-free"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "an unknown efp.response must be refused, not resolved to auto")
      if (allocated(error)) return
      call check(error, index(parse_error%get_message(), "matrix_free") > 0, &
                 "the refusal must name the spellings that would have worked: "// &
                 parse_error%get_message())
   end subroutine test_efp_keywords

   subroutine test_neo_keywords(error)
      !! `keywords.neo`: absent, by index, by symbol, and the block that names nothing
      !!
      !! The indices are 0-based in the deck and 1-based in the config, which is
      !! the convention every other atom list here follows; the test pins the
      !! shift so a deck's `[0]` can never quietly mean the second atom.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"scf": {"maxiter": 40}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%neo_active, "a deck with no neo block quantised something")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"neo": {"quantum_nuclei": [1], "nuclear_basis": "pb5-d"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%neo_active, "neo.quantum_nuclei by index was not read")
      if (allocated(error)) return
      call check(error, allocated(config%neo_quantum_indices), "the index list is missing")
      if (allocated(error)) return
      call check(error, size(config%neo_quantum_indices) == 1 .and. &
                 config%neo_quantum_indices(1) == 2, "deck index 1 must be the second atom")
      if (allocated(error)) return
      call check(error, trim(config%neo_nuclear_basis) == "pb5-d", "neo.nuclear_basis was not read")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"neo": {"quantum_nuclei": ["H"]}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%neo_quantum_symbols), "the symbol list is missing")
      if (allocated(error)) return
      call check(error, trim(config%neo_quantum_symbols(1)) == "H", "the symbol was not read")
      if (allocated(error)) return
      call check(error, trim(config%neo_nuclear_basis) == "pb4-d", &
                 "a silent nuclear_basis must be PB4-D")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"neo": {"nuclear_basis": "pb4-d"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a neo block naming no nucleus was accepted")
      if (allocated(error)) return

      call write_deck('"method": "dft", "basis": "6-31g", "functional": "b3lyp"', "Energy", &
                      '"neo": {"quantum_nuclei": ["H"], "epc": "17-2"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, trim(config%neo_epc) == "17-2", "neo.epc was not read")
      if (allocated(error)) return

      call write_deck('"method": "dft", "basis": "6-31g", "functional": "b3lyp"', "Energy", &
                      '"neo": {"quantum_nuclei": ["H"], "epc": "19"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "an unknown epc functional was accepted")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "6-31g"', "Energy", &
                      '"neo": {"quantum_nuclei": [0], "nuclear_bases": "pb4-d"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a misspelt neo key was accepted")
   end subroutine test_neo_keywords

   subroutine test_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Hessian", &
                      '"hessian": {"finite_difference_displacement": 0.005, '// &
                      '"temperature": 300.0, "pressure": 1.5, '// &
                      '"response_tolerance": 1.0e-6, "response_max_iter": 20, '// &
                      '"response_batch": 7}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%calc_type, CALC_TYPE_HESSIAN)
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_displacement, 0.005_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_temperature, 300.0_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_pressure, 1.5_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_response_tol, 1.0e-6_dp))
      if (allocated(error)) return
      call check(error, config%hessian_response_max_iter, 20)
      if (allocated(error)) return
      call check(error, config%hessian_response_batch, 7)
   end subroutine test_hessian

   subroutine test_hessian_defaults(error)
      !! An empty hessian block leaves every default in place
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Hessian", '"hessian": {}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_displacement, DEFAULT_DISPLACEMENT))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_temperature, DEFAULT_TEMPERATURE))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_pressure, DEFAULT_PRESSURE))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_response_tol, DEFAULT_RESPONSE_TOL))
      if (allocated(error)) return
      call check(error, config%hessian_response_max_iter, DEFAULT_RESPONSE_MAX_ITER)
   end subroutine test_hessian_defaults

   subroutine test_aimd(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"aimd": {"dt": 0.5, "nsteps": 1000, '// &
                      '"initial_temperature": 350.0, "output_frequency": 10}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%aimd_dt, 0.5_dp))
      if (allocated(error)) return
      call check(error, config%aimd_nsteps, 1000)
      if (allocated(error)) return
      call check(error, close_enough(config%aimd_initial_temperature, 350.0_dp))
      if (allocated(error)) return
      call check(error, config%aimd_output_frequency, 10)
   end subroutine test_aimd

   subroutine test_fragmentation(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "gmbe", "level": 3, '// &
                      '"max_intersection_level": 5, "embedding": "none", '// &
                      '"bond_breaking": "caps", "cap_scale": 0.71, '// &
                      '"cutoff_method": "distance", "distance_metric": "min"}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      ! `gmbe` rather than `MBE` + allow_overlapping_fragments: the flag is
      ! gone and the method names the expansion on its own now.
      call check(error, config%frag_method, "gmbe")
      if (allocated(error)) return
      call check(error, config%frag_level, 3)
      if (allocated(error)) return
      call check(error, config%max_intersection_level, 5)
      if (allocated(error)) return
      call check(error, config%embedding, "none")
      if (allocated(error)) return
      call check(error, config%bond_breaking, "caps")
      if (allocated(error)) return
      call check(error, abs(config%cap_scale - 0.71_dp) < 1.0e-12_dp, &
                 "cap_scale should read back as given")
      if (allocated(error)) return
      call check(error, config%cutoff_method, "distance")
      if (allocated(error)) return
      call check(error, config%distance_metric, "min")
   end subroutine test_fragmentation

   subroutine test_bond_breaking_defaults(error)
      !! A deck that names neither key keeps the behaviour this program had
      !! before they existed: cut bonds refused, and a cap -- when some other
      !! path places one -- exactly where the removed atom was.
      !!
      !! Worth a case of its own rather than trusting the declaration. These two
      !! defaults are what makes adding the keys a no-op for every deck already
      !! written, and `cap_scale` in particular is the difference between the
      !! cap gradient being the whole chain rule and being an approximation.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 2}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. allocated(config%bond_breaking), &
                 "bond_breaking should be unset when the deck is silent")
      if (allocated(error)) return
      call check(error, abs(config%cap_scale - 1.0_dp) < 1.0e-12_dp, &
                 "cap_scale should default to placing the cap at the removed atom")
   end subroutine test_bond_breaking_defaults

   subroutine test_df_without_aux(error)
      !! A density-fitted SCF with no auxiliary basis is refused, not defaulted
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      ! density_fitting on, no model.aux_basis -> must fail
      call write_deck('"method": "HF", "basis": "6-31g"', "Energy", &
                      '"scf": {"density_fitting": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "density_fitting with no aux basis should be refused")
      if (allocated(error)) return

      ! the same deck naming an aux basis is accepted
      call write_deck('"method": "HF", "basis": "6-31g", "aux_basis": "def2-universal-jkfit"', &
                      "Energy", '"scf": {"density_fitting": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
   end subroutine test_df_without_aux

   subroutine test_fmo_scf_keywords(error)
      !! FMO's inner per-fragment SCF controls reach the config from the deck
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "fmo", "level": 2, '// &
                      '"scf_max_iter": 250, '// &
                      '"scf_energy_tolerance": 1.0e-10, '// &
                      '"scf_density_tolerance": 1.0e-8}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fmo_scf_max_iter, 250)
      if (allocated(error)) return
      call check(error, close_enough(config%fmo_scf_energy_tol, 1.0e-10_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%fmo_scf_density_tol, 1.0e-8_dp))
      if (allocated(error)) return

      ! Absent keys leave the defaults untouched.
      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "fmo", "level": 2}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fmo_scf_max_iter, 100)
      if (allocated(error)) return
      call check(error, close_enough(config%fmo_scf_energy_tol, 1.0e-9_dp))
   end subroutine test_fmo_scf_keywords

   subroutine test_efmo_keywords(error)
      !! EFMO's two settings reach the config, and their defaults survive silence
      !!
      !! They come from two blocks on purpose: `rcut` decides which pairs get a
      !! quantum dimer, which is a property of the partition and so sits in
      !! `keywords.fragmentation` beside FMO's `resppc`; `charge_transfer` says
      !! what EFMO does with the far pairs and sits in `keywords.efmo`. A case
      !! rather than trust, because a key added to the wrong allow-list is
      !! refused by the schema and a key read from the wrong path is silently
      !! ignored.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "HF", "basis": "6-31g"', "Energy", &
                      '"fragmentation": {"method": "efmo", "level": 2, '// &
                      '"rcut": 1.25}, "efmo": {"charge_transfer": false}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%frag_method, "efmo")
      if (allocated(error)) return
      call check(error, close_enough(config%efmo_rcut, 1.25_dp))
      if (allocated(error)) return
      call check(error,.not. config%efmo_charge_transfer, &
                 "charge_transfer: false should switch the far-pair CT term off")
      if (allocated(error)) return

      ! Silence leaves the paper's defaults: R_cut = 2.0, charge transfer on as
      ! GAMESS's EFMO has it.
      call write_deck('"method": "HF", "basis": "6-31g"', "Energy", &
                      '"fragmentation": {"method": "efmo", "level": 2}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%efmo_rcut, 2.0_dp))
      if (allocated(error)) return
      call check(error, config%efmo_charge_transfer, &
                 "charge transfer should default to on")
   end subroutine test_efmo_keywords

   subroutine test_efmo_rcut_refused(error)
      !! `rcut` at or below zero is refused rather than run
      !!
      !! It is a ratio of a separation to a van der Waals contact, so unlike
      !! FMO's `resppc` -- where negative means "no approximation" -- there is
      !! no reading under which a non-positive value is a request. It would put
      !! every pair in the effective list, which is EFP with in-vacuo monomers
      !! and a different method from the one the deck named.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "HF", "basis": "6-31g"', "Energy", &
                      '"fragmentation": {"method": "efmo", "level": 2, '// &
                      '"rcut": 0.0}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "rcut: 0.0 should be refused")
      if (allocated(error)) return

      call write_deck('"method": "HF", "basis": "6-31g"', "Energy", &
                      '"fragmentation": {"method": "efmo", "level": 2, '// &
                      '"rcut": -1.0}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a negative rcut should be refused")
   end subroutine test_efmo_rcut_refused

   subroutine test_cutoffs(error)
      !! Named and numeric n-mer keys both land at the right level
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"cutoffs": {"dimer": 10.0, "trimer": 8.0}}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%fragment_cutoffs), "cutoffs should be allocated")
      if (allocated(error)) return
      call check(error, close_enough(config%fragment_cutoffs(2), 10.0_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%fragment_cutoffs(3), 8.0_dp))
      if (allocated(error)) return
      ! A level nobody gave keeps the sentinel, which downstream reads as
      ! "no cutoff here" -- a zero would screen every tetramer out instead.
      call check(error, config%fragment_cutoffs(4) < 0.0_dp, &
                 "an unspecified level keeps the negative sentinel")
      if (allocated(error)) return

      ! The numeric spelling means the same thing.
      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"cutoffs": {"2": 10.0, "3": 8.0}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%fragment_cutoffs(3), 8.0_dp), &
                 "numeric cutoff keys should reach the same level")
   end subroutine test_cutoffs

   subroutine test_cutoffs_monotonic(error)
      !! A trimer cutoff above the dimer one screens in more trimers than
      !! dimers, which cannot be intended
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"cutoffs": {"dimer": 4.0, "trimer": 9.0}}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error, parse_error%has_error(), &
                 "cutoffs that increase with level should be rejected")
   end subroutine test_cutoffs_monotonic

   subroutine test_global_groups(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "global_groups": 3}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%global_groups, 3)
   end subroutine test_global_groups

   subroutine test_nodes_per_group(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "nodes_per_group": 4}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nodes_per_group, 4)
   end subroutine test_nodes_per_group

   subroutine test_backend_keyword(error)
      !! The root-level `backend` key, and that an unknown one is refused
      !!
      !! Refused rather than defaulted, because a typo that fell back to "auto"
      !! would be a deck asking for one backend and silently getting another --
      !! and the two are independent implementations, so that is a provenance
      !! error, not a preference ignored.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: kind

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, trim(config%backend) == "auto", "the default backend is auto")
      if (allocated(error)) return

      call parse_backend_name("gpu", kind, parse_error)
      call check(error, kind == BACKEND_CUEST, "'gpu' must mean cuEST")
      if (allocated(error)) return
      call parse_backend_name("cpu", kind, parse_error)
      call check(error, kind == BACKEND_CZT, "'cpu' must mean libcint")
      if (allocated(error)) return
      call parse_backend_name("", kind, parse_error)
      call check(error, kind == BACKEND_AUTO, "an empty name must mean auto")
      if (allocated(error)) return

      call parse_error%clear()
      call parse_backend_name("cuda", kind, parse_error)
      call check(error, parse_error%has_error(), &
                 "an unknown backend name must be refused, not defaulted")
   end subroutine test_backend_keyword

   subroutine test_gpu_keyword(error)
      !! `system.gpu`, and that a GPU this build cannot give is refused at read time
      !!
      !! The refusal has to happen while the deck is still a deck. Left to the
      !! calculation it would arrive once per fragment, after fragmentation, from
      !! a bridge that computes nothing -- which is a much longer way round to
      !! the same sentence.
      !!
      !! Written to hold on either build. Whether a GPU can be had is a property
      !! of what linked, so the test asks the same question the reader does
      !! rather than assuming an answer.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      ! Absent: the backend request is left alone.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%gpu_set, "an absent system.gpu must not be 'seen'")
      if (allocated(error)) return
      call check(error, trim(config%backend) == "auto", "an absent system.gpu leaves auto")
      if (allocated(error)) return

      ! False: the CPU, named, on any build.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", &
                      '"gpu": false', two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, trim(config%backend) == "libcint", &
                 "system.gpu false must pin the CPU backend")
      if (allocated(error)) return

      ! True: honoured where cuEST is built, refused where it is not.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", &
                      '"gpu": true', two_atoms())
      call read_deck(config, parse_error)
      if (cuest_backend_available()) then
         call check(error,.not. parse_error%has_error(), parse_error%get_message())
         if (allocated(error)) return
         call check(error, trim(config%backend) == "cuest", &
                    "system.gpu true must ask for cuEST")
      else
         call check(error, parse_error%has_error(), &
                    "system.gpu true must be refused on a build without cuEST")
      end if
      if (allocated(error)) return

      ! A method with no GPU implementation is refused too. On a build without
      ! cuEST the missing backend is what stops it first; either way the deck
      ! does not run, which is the property being fixed here.
      call write_deck('"method": "gfn2"', "Energy", "", '"gpu": true', two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "system.gpu true must be refused for a method cuEST cannot run")
   end subroutine test_gpu_keyword

   subroutine test_gpu_conflict(error)
      !! `system.gpu` and `backend` naming different things is refused
      !!
      !! Refused rather than given a precedence rule: they are two spellings of
      !! one choice, and a deck that says both while meaning opposite things has
      !! no reading that is not a guess. Agreeing is fine -- redundant, not
      !! contradictory.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", &
                      '"gpu": true', two_atoms(), root='"backend": "cpu"')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "system.gpu true beside backend 'cpu' must be refused")
      if (allocated(error)) return

      ! The other direction, and on a build without cuEST it must still be the
      ! disagreement that is reported rather than the missing backend.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", &
                      '"gpu": false', two_atoms(), root='"backend": "gpu"')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "system.gpu false beside backend 'gpu' must be refused")
      if (allocated(error)) return
      call check(error, index(parse_error%get_message(), "disagree") > 0, &
                 "the message must name the disagreement: "//parse_error%get_message())
      if (allocated(error)) return

      ! Agreeing is not a conflict.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", &
                      '"gpu": false', two_atoms(), root='"backend": "cpu"')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
   end subroutine test_gpu_conflict

   subroutine test_cuest_methods(error)
      !! Which methods cuEST is allowed to be asked for
      !!
      !! An allow-list, so a method added later is CPU-only until someone writes
      !! it on the GPU. Checked here because it is the half of the GPU gate a
      !! build without cuEST cannot exercise through a deck -- the missing
      !! backend stops those decks first.
      type(error_type), allocatable, intent(out) :: error

      call check(error, method_runs_on_cuest(METHOD_TYPE_HF), "cuEST runs Hartree-Fock")
      if (allocated(error)) return
      call check(error, method_runs_on_cuest(METHOD_TYPE_DFT), "cuEST runs Kohn-Sham")
      if (allocated(error)) return
      call check(error,.not. method_runs_on_cuest(METHOD_TYPE_MP2), &
                 "MP2 has no cuEST implementation")
      if (allocated(error)) return
      call check(error,.not. method_runs_on_cuest(METHOD_TYPE_CCSD_T), &
                 "CCSD(T) has no cuEST implementation")
      if (allocated(error)) return
      call check(error,.not. method_runs_on_cuest(METHOD_TYPE_GFN2), &
                 "the semi-empirical methods never reach cuEST")
   end subroutine test_cuest_methods

   subroutine test_pcm_keywords(error)
      !! keywords.pcm, and that the block's presence is what switches it on
      !!
      !! There is no `enabled` flag inside the block on purpose: a deck that names
      !! a dielectric wants solvent, and a separate switch would let the two
      !! disagree. So an absent block must leave the calculation in the gas phase
      !! and a present one must turn the continuum on, which is the pair below.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "dft", "functional": "pbe0", "basis": "sto-3g"', &
                      "Energy", '"pcm": {"dielectric": 78.39, "angular_points": 302, '// &
                      '"radii_scale": 1.1, "max_iter": 50}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%pcm_enabled, "naming the block must enable the continuum")
      if (allocated(error)) return
      call check(error, abs(config%pcm_dielectric - 78.39_dp) < 1.0e-10_dp)
      if (allocated(error)) return
      call check(error, config%pcm_angular_points, 302)
      if (allocated(error)) return
      call check(error, abs(config%pcm_radii_scale - 1.1_dp) < 1.0e-10_dp)
      if (allocated(error)) return
      call check(error, config%pcm_max_iter, 50)
      if (allocated(error)) return
      ! Untouched keys keep their defaults rather than being zeroed.
      call check(error, abs(config%pcm_zeta - 2.0_dp) < 1.0e-10_dp, &
                 "an unmentioned key must keep its default")
      if (allocated(error)) return

      ! And no block at all is the gas phase.
      call write_deck('"method": "dft", "functional": "pbe0", "basis": "sto-3g"', &
                      "Energy", '"dft": {"grid_level": 3}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%pcm_enabled, &
                 "without the block the calculation must stay in the gas phase")
   end subroutine test_pcm_keywords

   subroutine test_dft_keywords(error)
      !! keywords.dft, and that explicit counts take the level out of charge
      !!
      !! The two ways of asking for a grid are exclusive: a level picks per-element
      !! counts from the standard tables, and explicit counts override it for every
      !! atom. Supplying counts has to switch the level off, or a run would report
      !! one grid and integrate on another -- which is what it used to do.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "dft", "functional": "b3lyp", "basis": "sto-3g"', &
                      "Energy", '"dft": {"grid_level": 5}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%dft_grid_level, 5)
      if (allocated(error)) return
      ! Absent counts stay negative, which is how the adapter knows the level is
      ! the thing that was asked for.
      call check(error, config%dft_radial_points < 0, "absent radial_points stays unset")
      if (allocated(error)) return

      call write_deck('"method": "dft", "functional": "pbe", "basis": "sto-3g"', &
                      "Energy", '"dft": {"radial_points": 99, "angular_points": 590}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%dft_radial_points, 99)
      if (allocated(error)) return
      call check(error, config%dft_angular_points, 590)
   end subroutine test_dft_keywords

   subroutine test_cc_spin_adapted(error)
      !! `keywords.cc.spin_adapted`, and that its default is on
      !!
      !! The default is the whole point of the case. Both formulations are exact
      !! for the closed-shell reference coupled cluster supports and agree to
      !! machine precision, so nothing in a deck's output would reveal which one
      !! ran -- only the time and the memory would move. A default flipped by
      !! accident is therefore a silent change of what every CCSD deck does.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      ! Absent: spin-adapted, which is the faster and smaller of the two.
      call write_deck('"method": "ccsd", "basis": "sto-3g"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%cc_spin_adapted, &
                 "coupled cluster must default to the spin-adapted formulation")
      if (allocated(error)) return

      ! Named false: the spin-orbital path, which is what a doubtful number is
      ! checked against.
      call write_deck('"method": "ccsd", "basis": "sto-3g"', "Energy", &
                      '"cc": {"spin_adapted": false}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%cc_spin_adapted, &
                 "keywords.cc.spin_adapted false was not read")
      if (allocated(error)) return

      ! And the allow-list is what let it through, not a validator that is off.
      call write_deck('"method": "ccsd", "basis": "sto-3g"', "Energy", &
                      '"cc": {"spin_adapated": false}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled cc key was accepted")
   end subroutine test_cc_spin_adapted

   subroutine test_cc_keywords(error)
      !! keywords.cc, and that "triples" records whether it was named
      !!
      !! The presence flag is the part worth testing. The triples default comes
      !! from the method name rather than from the field's initialiser -- "ccsd"
      !! and "ccsd(t)" being separate method types -- so a deck writing
      !! `"triples": false` has to be distinguishable from one that said nothing,
      !! or the name could never be overridden downward.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "ccsd(t)", "basis": "sto-3g"', "Energy", &
                      '"cc": {"maxiter": 40, "tolerance": 1.0e-10, '// &
                      '"triples": false, "diis": false, "diis_size": 4}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%cc_maxiter, 40)
      if (allocated(error)) return
      call check(error, close_enough(config%cc_tolerance, 1.0e-10_dp))
      if (allocated(error)) return
      call check(error, config%cc_diis, .false., "cc.diis should be settable")
      if (allocated(error)) return
      call check(error, config%cc_diis_size, 4)
      if (allocated(error)) return
      call check(error, config%cc_triples, .false., "an explicit triples:false is read")
      if (allocated(error)) return
      call check(error, config%cc_triples_set, .true., &
                 "naming triples must be recorded as named")
      if (allocated(error)) return

      ! Absent, so the method name decides and the adapter must not be overridden.
      call write_deck('"method": "ccsd(t)", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%cc_triples_set, .false., &
                 "an absent triples key must not read as named")
      if (allocated(error)) return
      call check(error, config%cc_maxiter, 100, "cc.maxiter keeps its default")
   end subroutine test_cc_keywords

   subroutine test_mcscf_keywords(error)
      !! keywords.mcscf, every key of it
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "casscf", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"n_active_electrons": 6, "n_active_orbitals": 6, '// &
                      '"n_inactive_orbitals": 4, "max_macro_iter": 120, '// &
                      '"orbital_convergence": 1.0e-7, "optimize_orbitals": false}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%mcscf_n_active_electrons, 6)
      if (allocated(error)) return
      call check(error, config%mcscf_n_active_orbitals, 6)
      if (allocated(error)) return
      call check(error, config%mcscf_n_inactive_orbitals, 4)
      if (allocated(error)) return
      call check(error, config%mcscf_max_macro_iter, 120)
      if (allocated(error)) return
      call check(error, close_enough(config%mcscf_orbital_convergence, 1.0e-7_dp))
      if (allocated(error)) return
      ! The keyword must beat the method name, which said "casscf".
      call check(error, config%mcscf_optimize_orbitals, .false., &
                 "an explicit optimize_orbitals:false must override the method name")
      if (allocated(error)) return

      ! Absent, so every field keeps its default and the name decides again.
      call write_deck('"method": "casscf", "basis": "cc-pvdz"', "Energy", "", "", &
                      two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%mcscf_n_active_electrons, 0, &
                 "an unset active space stays unset rather than being guessed")
      if (allocated(error)) return
      call check(error, config%mcscf_n_inactive_orbitals, -1, &
                 "n_inactive_orbitals defaults to 'derive it'")
      if (allocated(error)) return
      call check(error, config%mcscf_max_macro_iter, 100)
      if (allocated(error)) return
      call check(error, config%mcscf_optimize_orbitals, .true., &
                 "casscf without the keyword optimises orbitals")
   end subroutine test_mcscf_keywords

   subroutine test_casci_spelling(error)
      !! "casci" and "casscf" are one method type and differ by this boolean
      !!
      !! `parse_method_string` maps both to METHOD_TYPE_MCSCF, so the spelling
      !! is the only place the distinction survives the parse -- the same
      !! situation the "ri-" prefix is in. If this default were lost there would
      !! be no way to ask for a CASCI at all except by writing the keyword.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "casci", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"n_active_electrons": 4, "n_active_orbitals": 4}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%method, METHOD_TYPE_MCSCF, &
                 "casci must parse to the MCSCF method type")
      if (allocated(error)) return
      call check(error, config%mcscf_optimize_orbitals, .false., &
                 "the casci spelling must leave the orbitals fixed")
      if (allocated(error)) return

      ! And the keyword still wins over the name in the other direction.
      call write_deck('"method": "casci", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"n_active_electrons": 4, "n_active_orbitals": 4, '// &
                      '"optimize_orbitals": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%mcscf_optimize_orbitals, .true., &
                 "an explicit optimize_orbitals:true must override the casci spelling")
   end subroutine test_casci_spelling

   subroutine test_solvation(error)
      !! Every xTB solvation knob, including the two only .mqc could reach
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"xtb": {"solvent": "water", "solvation_model": "cpcm", '// &
                      '"dielectric": 80.0, "cpcm_nang": 230, "cpcm_rscale": 1.2, '// &
                      '"use_cds": false, "use_shift": false}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%solvent, "water")
      if (allocated(error)) return
      call check(error, config%solvation_model, "cpcm")
      if (allocated(error)) return
      call check(error, close_enough(config%dielectric, 80.0_dp))
      if (allocated(error)) return
      call check(error, config%cpcm_nang, 230)
      if (allocated(error)) return
      call check(error, close_enough(config%cpcm_rscale, 1.2_dp))
      if (allocated(error)) return
      call check(error, config%use_cds, .false., "use_cds should be settable from JSON")
      if (allocated(error)) return
      call check(error, config%use_shift, .false., "use_shift should be settable from JSON")
      if (allocated(error)) return

      ! Both default to on, so a deck that says nothing keeps the terms.
      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"xtb": {"solvent": "water"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%use_cds, .true., "use_cds defaults to true")
      if (allocated(error)) return
      call check(error, config%use_shift, .true., "use_shift defaults to true")
   end subroutine test_solvation

   subroutine test_inline_geometry(error)
      !! symbols + a flat coordinate list, the alternative to an xyz reference
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H"], '// &
                      '"geometry": [0.0, 0.0, 0.117, 0.0, 0.757, -0.469, '// &
                      '0.0, -0.757, -0.469], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%geometry%natoms, 3)
      if (allocated(error)) return
      call check(error, trim(config%geometry%elements(1)), "O")
      if (allocated(error)) return
      call check(error, close_enough(config%geometry%coords(3, 1), 0.117_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%geometry%coords(2, 3), -0.757_dp))
      if (allocated(error)) return

      ! A coordinate list that does not match the symbol count is an error,
      ! not a silently truncated molecule.
      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H"], "geometry": [0.0, 0.0, 0.1], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a short geometry list should be rejected")
   end subroutine test_inline_geometry

   subroutine test_missing_schema(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "model": {"method": "XTB-GFN2"},'
      write (unit, "(A)") '  "driver": "Energy",'
      write (unit, "(A)") '  "molecules": [{'//two_atoms()//'}]'
      write (unit, "(A)") "}"
      close (unit)

      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a deck without schema should fail")
   end subroutine test_missing_schema

   subroutine test_missing_molecules(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "schema": {"name": "mqc-frag", "version": "1.0"},'
      write (unit, "(A)") '  "model": {"method": "XTB-GFN2"},'
      write (unit, "(A)") '  "driver": "Energy"'
      write (unit, "(A)") "}"
      close (unit)

      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a deck without molecules should fail")
   end subroutine test_missing_molecules

   subroutine test_fragment_potentials(error)
      !! One potential per fragment, and a fragment left quantum by omitting one
      !!
      !! The list is parallel to `fragments`, like the charges and multiplicities
      !! beside it. An entry that is absent or empty leaves that fragment without a
      !! potential, which is not an error: it is how a mixed QM/EFP system is written,
      !! and there is deliberately no second way to declare a fragment.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H", "O", "H", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, -0.7,0,0, '// &
                      '3,0,0, 3.7,0,0, 2.3,0,0], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1, '// &
                      '"fragments": [[0, 1, 2], [3, 4, 5]], '// &
                      '"fragment_potentials": ["water.efp", ""]')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nfrag, 2)
      if (allocated(error)) return
      call check(error, allocated(config%fragments(1)%potential), &
                 "the first fragment should carry a potential")
      if (allocated(error)) return
      ! Resolved against the deck's own directory, exactly as `xyz` is, so a deck
      ! names the potential beside it and runs from anywhere. The test deck lives in
      ! the working directory, hence the "./".
      call check(error, config%fragments(1)%potential, "./water.efp")
      if (allocated(error)) return
      call check(error,.not. allocated(config%fragments(2)%potential), &
                 "an empty entry should leave the fragment quantum")
   end subroutine test_fragment_potentials

   subroutine test_uniform_system(error)
      !! `uniform_system` lets one potential stand for every fragment
      !!
      !! The point of it is a cluster with thousands of identical fragments, where a
      !! parallel list would be one filename repeated and nothing else. What is checked
      !! here is that the single string reaches every fragment, not just the first.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H", "O", "H", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, -0.7,0,0, '// &
                      '3,0,0, 3.7,0,0, 2.3,0,0], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1, '// &
                      '"fragments": [[0, 1, 2], [3, 4, 5]], '// &
                      '"uniform_system": true, '// &
                      '"fragment_potentials": "water.efp"')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%uniform_system, "uniform_system should be set")
      if (allocated(error)) return
      ! Both resolved against the deck, not just the first -- broadcasting must not
      ! bypass the path handling.
      call check(error, config%fragments(1)%potential, "./water.efp")
      if (allocated(error)) return
      call check(error, config%fragments(2)%potential, "./water.efp")
   end subroutine test_uniform_system

   subroutine test_uniform_rejected(error)
      !! A deck that claims uniformity and does not have it is refused
      !!
      !! `uniform_system` is an assertion about the system, not just a spelling
      !! convenience, so it is checked rather than believed. Fragments of different
      !! sizes cannot all be the same species, and the failure this prevents is the
      !! quiet one: every fragment handed a potential describing a different molecule,
      !! producing an energy rather than an error.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H", "O", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, -0.7,0,0, 3,0,0, 3.7,0,0], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1, '// &
                      '"fragments": [[0, 1, 2], [3, 4]], '// &
                      '"uniform_system": true, '// &
                      '"fragment_potentials": "water.efp"')
      call read_deck(config, parse_error)

      call check(error, parse_error%has_error(), &
                 "fragments of different sizes cannot be a uniform system")
   end subroutine test_uniform_rejected

   subroutine test_charges_property(error)
      !! `properties.charges`, and that the object is the request
      !!
      !! The shape matters more than the value read. `scheme` is not what asks
      !! for charges -- the object is -- so a deck naming no scheme still gets
      !! them, and a deck with no `charges` object gets none however else it is
      !! configured. Defaulting the scheme downstream instead would leave the
      !! config unable to say which partition a run used, and the report and
      !! the JSON both quote it.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      ! Absent object: no charges, and nothing to echo.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), "")
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. allocated(config%charges_scheme), &
                 "no charges object means no charges")
      if (allocated(error)) return

      ! Empty object: charges, at the documented default.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"charges": {}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%charges_scheme), &
                 "the object alone should ask for charges")
      if (allocated(error)) return
      call check(error, trim(config%charges_scheme), "mulliken", &
                 "and should default to the scheme the schema documents")
      if (allocated(error)) return

      ! Named scheme: taken as written.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"charges": {"scheme": "chelpg"}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, trim(config%charges_scheme), "chelpg", &
                 "a named scheme should be what the deck asked for")
      if (allocated(error)) return

      ! An unknown key inside the object is refused, like every other object.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"charges": {"schema": "chelpg"}}')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled key inside properties.charges should be refused")
   end subroutine test_charges_property

   subroutine test_bonding_analysis(error)
      !! `properties.bonding_analysis`, and that an unknown one is refused
      !!
      !! `properties` sits beside `keywords` rather than inside it, and the
      !! distinction is the reason this test exists at all: `keywords` change
      !! the number that comes out, `properties` ask for something further to be
      !! done with a wave function that is already determined. So a bonding
      !! analysis leaves `driver` as "energy" and changes no energy.
      !!
      !! The refusal is the more important half. A deck naming an analysis
      !! nobody implements would otherwise run a perfectly good energy and print
      !! nothing, which reads as the analysis having found nothing to say rather
      !! than as never having run.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "energy_threshold": 2.5}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%bonding_analysis), &
                 "the analysis name should have been read")
      if (allocated(error)) return
      call check(error, trim(config%bonding_analysis), "gms_quao", &
                 "and should be what the deck asked for")
      if (allocated(error)) return
      call check(error, abs(config%bonding_threshold - 2.5_dp) < 1.0e-12_dp, &
                 "the reporting threshold should have been read")
      if (allocated(error)) return

      call check(error, config%bonding_no_sharing .eqv. .false., &
                 "the no-sharing analysis should be off unless asked for")
      if (allocated(error)) return
      call check(error, config%bonding_energy .eqv. .false., &
                 "the energy decomposition should be off unless asked for")
      if (allocated(error)) return

      ! Opt-in because it needs the dense two-electron integrals, which the
      ! bonding tables on their own do not.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "energy_decomposition": true}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%bonding_energy, &
                 "the deck asked for the energy decomposition")
      if (allocated(error)) return

      ! The no-sharing analysis is opt-in because it costs a full valence CI,
      ! so the deck asking for it and the deck not mentioning it must differ.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "no_sharing": true}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%bonding_no_sharing, &
                 "the deck asked for the no-sharing analysis")
      if (allocated(error)) return
      call check(error,.not. allocated(config%bonding_no_sharing_ci), &
                 "an unmentioned CI route should be left for the default to fill")
      if (allocated(error)) return

      ! How that CI is obtained is a separate choice from whether to run it. The
      ! two routes describe the same wave function, so this picks how a number
      ! is computed rather than which number it is.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "no_sharing": true, '// &
                      '"no_sharing_ci": "resolve"}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%bonding_no_sharing_ci), &
                 "the deck named a CI route")
      if (allocated(error)) return
      call check(error, trim(config%bonding_no_sharing_ci), "resolve", &
                 "the CI route did not survive the parse")
      if (allocated(error)) return

      ! Constraining the localization is a separate opt-in, and off by default,
      ! because it changes the orbitals rather than only how they are reached.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "no_sharing": true, '// &
                      '"restrict_localization": true}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%bonding_restrict_localization, &
                 "the deck asked to constrain the localization")
      if (allocated(error)) return

      ! A route nobody implements is refused with both spellings named, rather
      ! than falling back to the default and computing something the deck did
      ! not ask for.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "no_sharing": true, '// &
                      '"no_sharing_ci": "davidson"}}')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "an unknown no-sharing CI route should be refused")
      if (allocated(error)) return

      ! Naming the analysis and nothing else keeps the default threshold.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao"}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, abs(config%bonding_threshold - 1.0_dp) < 1.0e-12_dp, &
                 "an unmentioned threshold should keep its default")
      if (allocated(error)) return
      call check(error, config%bonding_no_sharing .eqv. .false., &
                 "an unmentioned no-sharing flag should keep its default")
      if (allocated(error)) return

      ! Absent leaves it unset, which is what "no analysis" looks like.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      if (allocated(config%bonding_analysis)) then
         call check(error, trim(config%bonding_analysis), "none", &
                    "an absent properties block must not request an analysis")
         if (allocated(error)) return
      end if

      ! An analysis nobody implements is refused, by name, with the list.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "nbo"}}')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "an unknown bonding analysis should be refused rather than "// &
                 "silently producing no output")
      if (allocated(error)) return

      ! And so is a key inside the block that nothing reads.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bond_analysis": '// &
                      '{"type": "gms_quao"}}')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled key inside properties should be refused")
      if (allocated(error)) return

      ! A misspelling one level deeper too. This is the case the subgroup was
      ! made for: with a bare string there was nowhere to put a setting, and
      ! nowhere for the validator to catch one that was put wrongly.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"type": "gms_quao", "threshold": 2.0}}')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled setting inside the analysis block should be refused")
      if (allocated(error)) return

      ! Naming the block without saying which analysis is a deck that means
      ! nothing, rather than one that means "the default analysis".
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", "", "", &
                      two_atoms(), '"properties": {"bonding_analysis": '// &
                      '{"energy_threshold": 2.0}}')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a settings-only analysis block should be refused")
   end subroutine test_bonding_analysis

   subroutine test_hess_end(error)
      !! `keywords.optimization.hess_end`, read and defaulted off
      !!
      !! The keyword asks for a Hessian at the converged geometry, so a deck
      !! that never mentions it must come back false -- defaulting the other
      !! way would add a Hessian to every optimization ever run, which is the
      !! expensive direction to be wrong in.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Optimize", &
                      '"optimization": {"hess_end": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%opt_hess_end, "the request should be read")
      if (allocated(error)) return

      ! Said explicitly, the other way.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Optimize", &
                      '"optimization": {"hess_end": false}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%opt_hess_end, "false should be read as false")
      if (allocated(error)) return

      ! Not mentioned at all. `optional_logical` leaves the default alone, so
      ! this checks the default is the one in the type and not a second copy
      ! written into the reader.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Optimize", &
                      '"optimization": {"max_steps": 50}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%opt_hess_end, &
                 "a deck that never asked should not get a Hessian")
   end subroutine test_hess_end

   subroutine test_opt_target(error)
      !! `keywords.optimization.target`, through the reader and the schema
      !!
      !! The key has to survive three separate gates before it means anything:
      !! the schema allow-list, which refuses what it does not know; the
      !! reader, which puts the string on the config; and the default, which
      !! stays unallocated when nothing asked, so the type's own
      !! `OPT_TARGET_MINIMUM` is what an existing deck still gets. The
      !! vocabulary itself is checked in `test_mqc_optimizer_types`; what is
      !! checked here is that the string arrives at all.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Optimize", &
                      '"optimization": {"target": "saddle", "algorithm": "prfo", '// &
                      '"coordinates": "dlc"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%opt_target), "the target should be read")
      if (allocated(error)) return
      call check(error, trim(config%opt_target) == "saddle", &
                 "the target should arrive as it was spelled")
      if (allocated(error)) return

      ! Not mentioned. `optional_string` leaves the component unallocated, and
      ! the adapter only assigns when it is allocated -- which is what keeps
      ! every optimization written before this keyword existed a minimisation.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Optimize", &
                      '"optimization": {"max_steps": 50}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. allocated(config%opt_target), &
                 "a deck that never asked should leave the target to its default")
      if (allocated(error)) return

      ! The schema gate. A misspelled key inside the optimization block is
      ! refused before the reader sees it, so a typo cannot silently become a
      ! minimisation reported as a saddle search.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Optimize", &
                      '"optimization": {"targett": "saddle"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled optimization key should be refused")
   end subroutine test_opt_target

   subroutine test_full_valence(error)
      !! `keywords.mcscf.full_valence`, and that it is the only space named
      !!
      !! The third way of saying which orbitals are active, beside counts and an
      !! AVAS block. Each of the three describes the whole space, so any two of
      !! them together leave the deck's meaning depending on which the reader
      !! reached first -- refused rather than resolved, as AVAS already is
      !! against counts.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "casci", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"full_valence": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%mcscf_full_valence, "the valence request should be read")
      if (allocated(error)) return

      ! Absent leaves it off, and the space has to be named some other way.
      call write_deck('"method": "casci", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"n_active_electrons": 6, "n_active_orbitals": 6}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%mcscf_full_valence, &
                 "no mention should leave it off")
      if (allocated(error)) return

      ! With counts as well.
      call write_deck('"method": "casci", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"full_valence": true, "n_active_electrons": 6, '// &
                      '"n_active_orbitals": 6}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "naming the space twice should be refused")
      if (allocated(error)) return

      ! With AVAS as well.
      call write_deck('"method": "casci", "basis": "cc-pvdz"', "Energy", &
                      '"mcscf": {"full_valence": true, '// &
                      '"avas": {"orbitals": ["N 2p"]}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "two automatic selections at once should be refused")
      if (allocated(error)) return
   end subroutine test_full_valence

   subroutine test_ormas_keywords(error)
      !! `keywords.mcscf.ormas`, the subspaces and their occupation windows
      !!
      !! Three lists of the same length: where each subspace starts, and the
      !! fewest and most electrons it may hold counting both spins. Lengths that
      !! disagree are a typo in the deck, so they are caught here with the key
      !! names rather than four layers down as a complaint about enumerating
      !! occupation classes.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"n_active_electrons": 6, "n_active_orbitals": 7, '// &
                      '"ormas": {"subspaces": [1, 4], "min_electrons": [4, 0], '// &
                      '"max_electrons": [6, 2]}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%mcscf_ormas_subspaces), &
                 "the partition should have been read")
      if (allocated(error)) return
      call check(error, size(config%mcscf_ormas_subspaces), 2, "two subspaces")
      if (allocated(error)) return
      call check(error, all(config%mcscf_ormas_subspaces == [1, 4]), "where they start")
      if (allocated(error)) return
      call check(error, all(config%mcscf_ormas_min_electrons == [4, 0]), "the minima")
      if (allocated(error)) return
      call check(error, all(config%mcscf_ormas_max_electrons == [6, 2]), "the maxima")
      if (allocated(error)) return

      ! Absent leaves the arrays unallocated, which is what a complete active
      ! space looks like from here.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"n_active_electrons": 6, "n_active_orbitals": 6}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. allocated(config%mcscf_ormas_subspaces), &
                 "no ormas block should leave the partition unset")
      if (allocated(error)) return

      ! Lists of different lengths.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"n_active_electrons": 6, "n_active_orbitals": 7, '// &
                      '"ormas": {"subspaces": [1, 4], "min_electrons": [4], '// &
                      '"max_electrons": [6, 2]}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "three lists of different lengths are not a partition")
      if (allocated(error)) return

      ! A window list with no subspaces to go with it.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"n_active_electrons": 6, "n_active_orbitals": 7, '// &
                      '"ormas": {"min_electrons": [4, 0], "max_electrons": [6, 2]}}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a partition with no subspaces should be refused")
      if (allocated(error)) return
   end subroutine test_ormas_keywords

   subroutine test_avas_keywords(error)
      !! `keywords.mcscf.avas`, and that it cannot coexist with explicit counts
      !!
      !! An active space can be named two ways -- by counts, or by the atomic
      !! orbitals it should be built from -- and giving both is refused rather
      !! than resolved by precedence. AVAS decides the space from the labels, so
      !! whichever counts a deck also wrote would be silently discarded, and a
      !! deck whose meaning depends on which key the reader reached first is
      !! worse than one that is rejected.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"avas": {"orbitals": ["N 2s", "N 2p"], '// &
                      '"threshold": 0.3}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%mcscf_avas_orbitals), &
                 "the orbital labels should have been read")
      if (allocated(error)) return
      call check(error, size(config%mcscf_avas_orbitals), 2, "two labels")
      if (allocated(error)) return
      call check(error, trim(config%mcscf_avas_orbitals(1)), "N 2s", "the first label")
      if (allocated(error)) return
      call check(error, trim(config%mcscf_avas_orbitals(2)), "N 2p", "the second")
      if (allocated(error)) return
      call check(error, abs(config%mcscf_avas_threshold - 0.3_dp) < 1.0e-12_dp, &
                 "and the threshold")
      if (allocated(error)) return

      ! Absent leaves the threshold at the published default.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"avas": {"orbitals": ["N 2p"]}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, abs(config%mcscf_avas_threshold - 0.2_dp) < 1.0e-12_dp, &
                 "an unmentioned threshold should keep the published default")
      if (allocated(error)) return

      ! Both ways of naming a space at once.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"avas": {"orbitals": ["N 2p"]}, '// &
                      '"n_active_electrons": 6, "n_active_orbitals": 6}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "naming the space twice should be refused rather than resolved")
      if (allocated(error)) return

      ! A block that selects on nothing.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"avas": {"orbitals": []}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "an empty orbital list should be refused")
      if (allocated(error)) return

      ! And a misspelled key inside the block.
      call write_deck('"method": "casscf", "basis": "sto-3g"', "Energy", &
                      '"mcscf": {"avas": {"orbitals": ["N 2p"], "cutoff": 0.2}}', &
                      "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a misspelled setting inside the avas block should be refused")
   end subroutine test_avas_keywords

   subroutine test_diis_keywords(error)
      !! `keywords.scf.diis` and `keywords.scf.diis_size`
      !!
      !! Both fields existed all the way down to `scf_config_t` and nothing
      !! read them, so every SCF ever run used a subspace of eight whatever the
      !! deck said -- and a deck that said something got no error either, since
      !! the schema simply did not know the keys. That is the failure this
      !! guards: a keyword that is accepted, plumbed, and ignored.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"diis_size": 16}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%scf_diis_size == 16, "the subspace size should be read")
      if (allocated(error)) return
      call check(error, config%scf_use_diis, "DIIS should still be on")
      if (allocated(error)) return

      ! Off, which is a diagnostic rather than a setting.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"diis": false}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%scf_use_diis, "false should be read as false")
      if (allocated(error)) return

      ! Never mentioned: the defaults are the ones in the type, not a second
      ! copy written into the reader.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 50}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%scf_use_diis, "DIIS defaults on")
      if (allocated(error)) return
      call check(error, config%scf_diis_size == 8, "the default subspace is eight")
   end subroutine test_diis_keywords

   subroutine test_incremental_fock_keyword(error)
      !! `keywords.scf.incremental_fock`
      !!
      !! The switch exists so that "is the incremental build doing this?" can be
      !! answered from a deck rather than by rebuilding the binary. That makes
      !! the default the load-bearing part: an SCF that quietly stopped building
      !! incrementally would lose most of its speed with nothing to say so, and
      !! one that could not be turned off leaves a stalling run with no way to
      !! rule the correction out.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"incremental_fock": false}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%scf_incremental_fock, "false should be read as false")
      if (allocated(error)) return

      ! Never mentioned: on, from the default in the type rather than a second
      ! copy written into the reader.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 50}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%scf_incremental_fock, "incremental building defaults on")
      if (allocated(error)) return

      ! Spelled true explicitly, which a deck pinning the default would do.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"incremental_fock": true}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%scf_incremental_fock, "true should be read as true")
   end subroutine test_incremental_fock_keyword

   subroutine test_fukui_scf_inherit(error)
      !! An absent `properties.fukui.scf` leaves the ions on `keywords.scf`
      !!
      !! This is the whole reason the group replaced three sentinel fields. The
      !! ions are a second SCF beside the neutral, and a deck that says nothing
      !! about them must converge them the way it asked for the neutral -- not
      !! the way the type defaults happen to read.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 77, "level_shift": 0.25, "diis_size": 3}', &
                      "", two_atoms(), root='"properties": {"fukui": {"population": "chelpg"}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fukui_scf%max_iter == 77, "maxiter did not inherit")
      if (allocated(error)) return
      call check(error, config%fukui_scf%level_shift == 0.25_dp, &
                 "level_shift did not inherit")
      if (allocated(error)) return
      call check(error, config%fukui_scf%diis_size == 3, "diis_size did not inherit")
      if (allocated(error)) return
      call check(error, config%fukui_scf%inherit_scf, "inherit_scf defaults true")
   end subroutine test_fukui_scf_inherit

   subroutine test_fukui_scf_override(error)
      !! A named key wins; its neighbours still inherit
      !!
      !! The half that a plain extended type could not express. Naming one
      !! setting for the ions must not silently reset the rest to type
      !! defaults, which is what would happen if "unset" were decided by
      !! comparing against a default rather than by seeding first.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 77, "level_shift": 0.25, "diis_size": 3}', &
                      "", two_atoms(), &
                      root='"properties": {"fukui": {"population": "chelpg", '// &
                      '"scf": {"maxiter": 250, "accelerator": "ediis"}}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fukui_scf%max_iter == 250, "the named maxiter must win")
      if (allocated(error)) return
      call check(error, config%fukui_scf%accelerator == "ediis", "accelerator")
      if (allocated(error)) return
      ! Not named, so still the neutral's -- not the 0.0 and 8 of the type.
      call check(error, config%fukui_scf%level_shift == 0.25_dp, &
                 "an unnamed neighbour must keep the inherited value")
      if (allocated(error)) return
      call check(error, config%fukui_scf%diis_size == 3, &
                 "an unnamed neighbour must keep the inherited value")
      if (allocated(error)) return
      ! The neutral is untouched by any of it.
      call check(error, config%scf_maxiter == 77, "the neutral must not move")
   end subroutine test_fukui_scf_override

   subroutine test_fukui_scf_guess_spelling(error)
      !! The ions inherit `keywords.guess.type`, not only `keywords.scf.guess`
      !!
      !! There are two spellings for the neutral's guess -- `keywords.scf.guess`
      !! is superseded and `keywords.guess.type` is current -- and
      !! `config_to_driver` resolves them in that order. Seeding from only the
      !! superseded one left a deck written the current way inheriting twelve
      !! keys and silently not the thirteenth: the ladder came across while the
      !! guess that selects it stayed `auto`.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      ! The current spelling, on its own.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"guess": {"type": "gwh"}', "", two_atoms(), &
                      root='"properties": {"fukui": {"population": "chelpg"}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fukui_scf%guess == "gwh", &
                 "keywords.guess.type must reach the ions")
      if (allocated(error)) return

      ! The superseded spelling still works.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"guess": "core"}', "", two_atoms(), &
                      root='"properties": {"fukui": {"population": "chelpg"}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fukui_scf%guess == "core", &
                 "keywords.scf.guess must still reach the ions")
      if (allocated(error)) return

      ! And the ions may still name their own, which beats either.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"guess": {"type": "gwh"}', "", two_atoms(), &
                      root='"properties": {"fukui": {"population": "chelpg", '// &
                      '"scf": {"guess": "sad"}}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%fukui_scf%guess == "sad", &
                 "the ions' own guess must win over the inherited one")
   end subroutine test_fukui_scf_guess_spelling

   subroutine test_fukui_scf_no_inherit(error)
      !! `inherit_scf: false` falls back to the defaults, not to the neutral
      !!
      !! For a neutral converged tightly, or with an accelerator chosen for it,
      !! where carrying that onto a charged species is the wrong start. What it
      !! changes is only the fallback: a key the object names still wins.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 77, "level_shift": 0.25, "diis_size": 3}', &
                      "", two_atoms(), &
                      root='"properties": {"fukui": {"population": "chelpg", '// &
                      '"scf": {"inherit_scf": false, "level_shift": 0.5}}}')
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%fukui_scf%inherit_scf, "inherit_scf must read false")
      if (allocated(error)) return
      call check(error, config%fukui_scf%max_iter == 100, &
                 "an unnamed key must fall back to the default, not the neutral")
      if (allocated(error)) return
      call check(error, config%fukui_scf%diis_size == 8, &
                 "an unnamed key must fall back to the default, not the neutral")
      if (allocated(error)) return
      call check(error, config%fukui_scf%level_shift == 0.5_dp, &
                 "a named key wins whatever inherit_scf says")
   end subroutine test_fukui_scf_no_inherit

   subroutine test_convergence_metric(error)
      !! `keywords.scf.convergence_metric` reaches the config
      !!
      !! The fifth step AGENTS.md asks for. A key added everywhere but here is
      !! one nothing would notice losing.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"convergence_metric": "commutator"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%scf_convergence_metric), &
                 "the metric must reach the config")
      if (allocated(error)) return
      call check(error, config%scf_convergence_metric == "commutator", "commutator")
      if (allocated(error)) return

      ! Unmentioned stays unallocated, so the default in the method layer is
      ! the single copy rather than one written again here.
      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 50}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. allocated(config%scf_convergence_metric), &
                 "an unmentioned metric must not be invented")
   end subroutine test_convergence_metric

   subroutine test_maxiter_named(error)
      !! Whether the deck named `keywords.scf.maxiter`, not just its value
      !!
      !! MAKEFP defaults to 200 iterations where the shared default is 100, so
      !! it has to tell "the deck asked for 100" from "the deck said nothing".
      !! Only the flag carries that, and the value alone cannot.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"maxiter": 100}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%scf_maxiter_set, &
                 "a deck naming maxiter must set the flag, even at the default value")
      if (allocated(error)) return

      call write_deck('"method": "hf", "basis": "sto-3g"', "Energy", &
                      '"scf": {"tolerance": 1e-8}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error,.not. config%scf_maxiter_set, &
                 "a deck that never mentions maxiter must not set the flag")
   end subroutine test_maxiter_named

   subroutine test_malformed(error)
      !! Broken JSON is a parse error, not a crash or a half-filled config
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") '{"schema": {"name": "mqc-frag",'
      close (unit)

      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "malformed JSON should fail cleanly")
   end subroutine test_malformed

   pure function close_enough(a, b) result(same)
      real(dp), intent(in) :: a, b
      logical :: same
      same = abs(a - b) <= 1.0e-10_dp*max(1.0_dp, abs(b))
   end function close_enough

   subroutine remove_deck()
      !! Delete the scratch deck
      !!
      !! Called from the driver rather than from each test, so it runs whether
      !! or not a test returned early. It matters more here than for most
      !! scratch files: one case deliberately writes malformed JSON, and a
      !! leftover copy trips the repository's check-json hook.
      integer :: unit
      logical :: exists

      inquire (file=DECK, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=DECK, status="old")
      close (unit, status="delete")
   end subroutine remove_deck

end module test_mqc_json_reader

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_reader, only: collect_mqc_json_reader_tests, remove_deck
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_json_reader", collect_mqc_json_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   call remove_deck()

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
