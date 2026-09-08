module test_mqc_config_adapter
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_config_types, only: input_fragment_t, mqc_config_t
   use mqc_config_adapter, only: check_fragment_overlap, config_to_driver, driver_config_t
   use mqc_optimizer_types, only: OPT_TARGET_MINIMUM, OPT_TARGET_SADDLE
   use mqc_method_types, only: METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_OPTIMIZE
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_config_adapter_tests

contains

   !! Collect all exported unit tests
   subroutine collect_mqc_config_adapter_tests(testsuite)
      !! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_overlap_detection", test_no_overlap), &
                  new_unittest("overlap_detection", test_overlap_detected), &
                  new_unittest("single_fragment_no_overlap", test_single_fragment), &
                  new_unittest("driver_global_groups", test_driver_global_groups), &
                  new_unittest("driver_nodes_per_group", test_driver_nodes_per_group), &
                  new_unittest("hessian_displacement_reaches_the_method", &
                               test_hessian_displacement), &
                  new_unittest("saddle_target_reaches_driver", test_saddle_target), &
                  new_unittest("saddle_needs_a_saddle_algorithm", test_saddle_algorithm), &
                  new_unittest("unknown_target_refused", test_unknown_target), &
                  new_unittest("driver_cartesian", test_driver_cartesian), &
                  new_unittest("driver_cartesian_default", test_driver_cartesian_default), &
                  new_unittest("driver_carries_maxiter_named_flag", test_driver_maxiter_named), &
                  new_unittest("fragmentation_method_selects_the_expansion", test_frag_method), &
                  new_unittest("fragmentation_method_refuses_nonsense", test_frag_method_bad) &
                  ]
   end subroutine collect_mqc_config_adapter_tests

   subroutine test_no_overlap(error)
      !! Test that non-overlapping fragments pass validation
      type(error_type), allocatable, intent(out) :: error
      type(input_fragment_t), allocatable :: fragments(:)
      type(error_t) :: parse_error

      ! Create two non-overlapping fragments
      allocate (fragments(2))

      ! Fragment 1: atoms 0, 1, 2
      allocate (fragments(1)%indices(3))
      fragments(1)%indices = [0, 1, 2]
      fragments(1)%charge = 0
      fragments(1)%multiplicity = 1

      ! Fragment 2: atoms 3, 4, 5
      allocate (fragments(2)%indices(3))
      fragments(2)%indices = [3, 4, 5]
      fragments(2)%charge = 0
      fragments(2)%multiplicity = 1

      ! Check for overlap - should find none
      call check_fragment_overlap(fragments, 2, parse_error)

      call check(error,.not. parse_error%has_error(), "Non-overlapping fragments should pass validation: "//parse_error%get_full_trace())
      if (allocated(error)) return

      ! Clean up
      call fragments(1)%destroy()
      call fragments(2)%destroy()
      deallocate (fragments)
   end subroutine test_no_overlap

   subroutine test_overlap_detected(error)
      !! Test that overlapping fragments are detected
      type(error_type), allocatable, intent(out) :: error
      type(input_fragment_t), allocatable :: fragments(:)
      type(error_t) :: parse_error

      ! Create two overlapping fragments
      allocate (fragments(2))

      ! Fragment 1: atoms 0, 1, 2
      allocate (fragments(1)%indices(3))
      fragments(1)%indices = [0, 1, 2]
      fragments(1)%charge = 0
      fragments(1)%multiplicity = 1

      ! Fragment 2: atoms 2, 3, 4 (overlaps on atom 2)
      allocate (fragments(2)%indices(3))
      fragments(2)%indices = [2, 3, 4]
      fragments(2)%charge = 0
      fragments(2)%multiplicity = 1

      ! Check for overlap - should detect it
      call check_fragment_overlap(fragments, 2, parse_error)

      call check(error, parse_error%has_error(), "Overlapping fragments should be detected")
      if (allocated(error)) return

      ! Check that error message mentions the overlapping atom
      call check(error, index(parse_error%get_message(), "atom 2") > 0, &
                 "Error message should mention overlapping atom")
      if (allocated(error)) return

      ! Clean up
      call fragments(1)%destroy()
      call fragments(2)%destroy()
      deallocate (fragments)
   end subroutine test_overlap_detected

   subroutine test_single_fragment(error)
      !! Test that a single fragment has no overlap (edge case)
      type(error_type), allocatable, intent(out) :: error
      type(input_fragment_t), allocatable :: fragments(:)
      type(error_t) :: parse_error

      ! Create a single fragment
      allocate (fragments(1))

      ! Fragment 1: atoms 0, 1, 2
      allocate (fragments(1)%indices(3))
      fragments(1)%indices = [0, 1, 2]
      fragments(1)%charge = 0
      fragments(1)%multiplicity = 1

      ! Check for overlap - single fragment should pass
      call check_fragment_overlap(fragments, 1, parse_error)

 call check(error,.not. parse_error%has_error(), "Single fragment should pass validation: "//parse_error%get_full_trace())
      if (allocated(error)) return

      ! Clean up
      call fragments(1)%destroy()
      deallocate (fragments)
   end subroutine test_single_fragment

   subroutine test_driver_global_groups(error)
      !! Test global_groups is copied into driver_config
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_ENERGY
      config%nfrag = 0
      config%global_groups = 2
      config%nodes_per_group = 0

      call config_to_driver(config, driver_config)

      call check(error, driver_config%global_groups, 2, "global_groups should be copied")
      if (allocated(error)) return

      call check(error, driver_config%nodes_per_group, 0, "nodes_per_group should default to 0")
   end subroutine test_driver_global_groups

   subroutine frag_driver(method, driver_config, err)
      !! A fragmented deck naming `method`, through the adapter
      character(len=*), intent(in) :: method
      type(driver_config_t), intent(out) :: driver_config
      type(error_t), intent(inout) :: err

      type(mqc_config_t) :: config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_ENERGY
      config%nfrag = 0
      config%frag_method = method
      call config_to_driver(config, driver_config, error=err)
   end subroutine frag_driver

   subroutine test_frag_method(error)
      !! Each method reaches the two switches the driver actually branches on
      !!
      !! `method` was required, validated for presence and then never read --
      !! the expansion came from `expansion` and `allow_overlapping_fragments`
      !! instead -- two more ways to say the same thing, both now deleted.
      !! Every fragmented deck in the tree said "MBE", including the ones
      !! running FMO. These assert the key now decides, and is the only thing
      !! that does.
      type(error_type), allocatable, intent(out) :: error
      type(driver_config_t) :: dc
      type(error_t) :: err

      call frag_driver("mbe", dc, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, trim(dc%expansion_kind) == "mbe" .and. &
                 .not. dc%allow_overlapping_fragments, "mbe")
      if (allocated(error)) return

      call frag_driver("fmo", dc, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, trim(dc%expansion_kind) == "fmo", "fmo")
      if (allocated(error)) return

      call frag_driver("ee-mbe", dc, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, trim(dc%expansion_kind) == "ee-mbe", "ee-mbe")
      if (allocated(error)) return

      ! GMBE is the overlapping flag with the plain term builder, which is how
      ! the driver has always told the two apart.
      call frag_driver("gmbe", dc, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, dc%allow_overlapping_fragments .and. &
                 trim(dc%expansion_kind) == "mbe", "gmbe")
      if (allocated(error)) return

      ! Case and separator are not significant: decks in this tree write both
      ! "MBE" and "mbe", and the hyphen is the kind of thing typed either way.
      call frag_driver("EE_MBE", dc, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, trim(dc%expansion_kind) == "ee-mbe", "EE_MBE should read as ee-mbe")
   end subroutine test_frag_method

   subroutine test_frag_method_bad(error)
      !! A name with no expansion behind it is refused
      !!
      !! And, beside it, that a name which *does* have one is not: the two
      !! outcomes come from the same table, so a method added to the parser and
      !! forgotten in `fragmentation_method_implemented` fails here.
      type(error_type), allocatable, intent(out) :: error
      type(driver_config_t) :: dc
      type(error_t) :: err

      call frag_driver("wibble", dc, err)
      call check(error, err%has_error(), "an unknown method must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "wibble") > 0, &
                 "the refusal must name the spelling used")
      if (allocated(error)) return
      call err%clear()

      ! `efmo` used to be refused here as reserved-but-unimplemented, which is
      ! a different fact from "no such method" and was checked as such. Phase 2
      ! built it, so what is asserted now is that it is accepted and that it
      ! resolves to its own expansion rather than falling back to `mbe` -- which
      ! is how a half-wired method would look: accepted, and running something
      ! else.
      call frag_driver("efmo", dc, err)
      call check(error,.not. err%has_error(), "efmo is implemented and must be "// &
                 "accepted: "//err%get_message())
      if (allocated(error)) return
      call check(error, trim(dc%expansion_kind) == "efmo", &
                 "efmo must resolve to its own expansion, not to "// &
                 trim(dc%expansion_kind))
      if (allocated(error)) return
      call check(error,.not. dc%allow_overlapping_fragments, &
                 "efmo fragments do not overlap")
      if (allocated(error)) return
      call err%clear()
   end subroutine test_frag_method_bad

   subroutine test_driver_maxiter_named(error)
      !! `keywords.scf.maxiter` and the flag saying the deck named it
      !!
      !! The adapter is the layer the MakeFP tests cannot reach: they call
      !! `make_efp_potential` directly, so removing this copy would leave every
      !! one of them passing while a deck's `maxiter` stopped arriving. That is
      !! the exact silent-drop path the work is meant to close, so the layer
      !! gets its own assertion.
      !!
      !! Both states matter. MAKEFP runs 200 iterations by default, larger than
      !! the shared 100 on purpose, so it needs to tell a deck that asked for
      !! 100 from one that said nothing -- and that distinction is the flag,
      !! not the value.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_ENERGY
      config%nfrag = 0

      ! Named.
      config%scf_maxiter = 137
      config%scf_maxiter_set = .true.
      call config_to_driver(config, driver_config)
      call check(error, driver_config%method_config%scf%max_iter == 137, &
                 "the value must reach the method config")
      if (allocated(error)) return
      call check(error, driver_config%method_config%scf%max_iter_set, &
                 "the named flag must reach the method config")
      if (allocated(error)) return

      ! Silent: the value still travels, but nothing may claim the deck asked.
      config%scf_maxiter_set = .false.
      call config_to_driver(config, driver_config)
      call check(error,.not. driver_config%method_config%scf%max_iter_set, &
                 "an unnamed maxiter must not arrive marked as named")
   end subroutine test_driver_maxiter_named

   subroutine test_driver_cartesian(error)
      !! `model.cartesian` reaches the SCF settings the backend reads
      !!
      !! Worth a test of its own because the consequence of it *not* arriving
      !! is silent: the run succeeds, converges and reports an energy in the
      !! basis the file declared rather than the one the deck asked for. For a
      !! shell above p those are different spaces, so it is a different answer
      !! and not a different spelling of one.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_ENERGY
      config%nfrag = 0
      config%cartesian = .true.

      call config_to_driver(config, driver_config)

      call check(error, driver_config%method_config%scf%cartesian, &
                 "model.cartesian should reach the SCF config")
   end subroutine test_driver_cartesian

   subroutine test_driver_cartesian_default(error)
      !! Absent, it is false, so a deck that says nothing honours the file
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_ENERGY
      config%nfrag = 0

      call config_to_driver(config, driver_config)

      call check(error,.not. driver_config%method_config%scf%cartesian, &
                 "cartesian should default to false")
   end subroutine test_driver_cartesian_default

   subroutine test_driver_nodes_per_group(error)
      !! Test nodes_per_group is copied into driver_config
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_ENERGY
      config%nfrag = 0
      config%global_groups = 0
      config%nodes_per_group = 4

      call config_to_driver(config, driver_config)

      call check(error, driver_config%nodes_per_group, 4, "nodes_per_group should be copied")
      if (allocated(error)) return

      call check(error, driver_config%global_groups, 0, "global_groups should default to 0")
   end subroutine test_driver_nodes_per_group

   subroutine optimize_config(config)
      !! The smallest deck that reaches the optimization vocabulary
      type(mqc_config_t), intent(out) :: config

      config%method = METHOD_TYPE_GFN2
      config%calc_type = CALC_TYPE_OPTIMIZE
      config%nfrag = 0
   end subroutine optimize_config

   subroutine test_hessian_displacement(error)
      !! The deck's finite-difference step has to reach the *method*
      !!
      !! `driver_config%hessian%displacement` was already set and already
      !! covered -- `test_mqc_json_reader` asserts the reader picks the value
      !! up. Nothing asserted the next rung, and that is exactly where it was
      !! dropped: `finite_difference_hessian` took no displacement argument at
      !! all, so every finite-difference Hessian ran at `DEFAULT_DISPLACEMENT`
      !! while the log printed that default back as though it had been asked
      !! for. A keyword accepted, validated, read, and then silently ignored.
      !!
      !! Both fields are checked. The driver one is what the vibrational path
      !! reads and the method one is what actually displaces the geometry, and
      !! they are set from the same source in two separate statements -- so one
      !! can be right while the other is not, which is the state this test was
      !! written to end.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config
      type(error_t) :: parse_error

      call optimize_config(config)
      config%hessian_displacement = 0.0125_dp

      call config_to_driver(config, driver_config, error=parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, abs(driver_config%hessian%displacement - 0.0125_dp) < 1.0e-12_dp, &
                 "the displacement should reach the driver")
      if (allocated(error)) return
      call check(error, &
                 abs(driver_config%method_config%hessian_displacement - 0.0125_dp) < 1.0e-12_dp, &
                 "the displacement should reach the method, which is what displaces")
   end subroutine test_hessian_displacement

   subroutine test_saddle_target(error)
      !! A saddle target with an algorithm that can look for one
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config
      type(error_t) :: parse_error

      call optimize_config(config)
      config%opt_target = "saddle"
      config%opt_algorithm = "prfo"
      config%opt_coordinates = "dlc"

      call config_to_driver(config, driver_config, error=parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, driver_config%optimization%target, OPT_TARGET_SADDLE, &
                 "the target should reach the driver")
   end subroutine test_saddle_target

   subroutine test_saddle_algorithm(error)
      !! Asked for a saddle, given a downhill optimizer
      !!
      !! Both spellings of the deck: one that names L-BFGS and one that names
      !! no algorithm at all and gets it as the default. The second is the case
      !! worth pinning -- it is the deck a first-time user writes, and it has to
      !! be refused rather than run to a minimum and reported as a saddle.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config
      type(error_t) :: parse_error

      call optimize_config(config)
      config%opt_target = "saddle"
      config%opt_algorithm = "lbfgs"

      call config_to_driver(config, driver_config, error=parse_error)
      call check(error, parse_error%has_error(), &
                 "a saddle search with a downhill algorithm should be refused")
      if (allocated(error)) return

      call optimize_config(config)
      config%opt_target = "saddle"

      call config_to_driver(config, driver_config, error=parse_error)
      call check(error, parse_error%has_error(), &
                 "a saddle search with the default algorithm should be refused too")
   end subroutine test_saddle_algorithm

   subroutine test_unknown_target(error)
      !! A target nobody recognises, with and without somewhere to report it
      !!
      !! The second half is the point. `error` is optional and the session
      !! omits it, so a parse failure that was written into the driver anyway
      !! would sit there as a sentinel below OPT_TARGET_MINIMUM and read as
      !! "not a saddle" -- a minimisation, run silently, from a deck that asked
      !! for something else.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver_config
      type(error_t) :: parse_error

      call optimize_config(config)
      config%opt_target = "banana"

      call config_to_driver(config, driver_config, error=parse_error)
      call check(error, parse_error%has_error(), &
                 "an unrecognised target should be refused by name")
      if (allocated(error)) return

      call optimize_config(config)
      config%opt_target = "banana"

      call config_to_driver(config, driver_config)
      call check(error, driver_config%optimization%target, OPT_TARGET_MINIMUM, &
                 "a target that failed to parse must not reach the driver")
   end subroutine test_unknown_target

end module test_mqc_config_adapter

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_config_adapter, only: collect_mqc_config_adapter_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_config_adapter", collect_mqc_config_adapter_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
