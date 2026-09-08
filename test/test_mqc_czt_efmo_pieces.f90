!! The pieces EFMO is assembled from, each checked on its own
module test_mqc_czt_efmo_pieces
   !! EFMO's energy (Sattasathuchana et al., JCTC 20, 2445 (2024), eq 6) is
   !!
   !!     E = sum_I E_I^0
   !!       + sum_{R_IJ <= R_cut} (E_IJ^0 - E_I^0 - E_J^0 - E_IJ^pol)
   !!       + sum_{R_IJ >  R_cut} (E_IJ^Coul + E_IJ^disp + E_IJ^ExRep + E_IJ^CT)
   !!       + E_pol^total
   !!
   !! and every piece of it already exists here in some form: MAKEFP runs the
   !! monomer SCF that `E_I^0` is, the EFP terms are validated against GAMESS,
   !! and FMO measures the `R_IJ` of eq 2. What was missing is the *interfaces*
   !! -- one pair at a time, one potential in memory, one SCF instead of two --
   !! and this file tests those five interfaces before anything is wired
   !! together, because each of them can be wrong in a way that shifts the total
   !! by a plausible amount and reports nothing.
   !!
   !! Water in 6-31G throughout, and the potential is built once and cached: a
   !! MAKEFP run is an SCF, a localization and twelve frequency-dependent
   !! response solves, and it is the whole cost of this file.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential, &
                                    write_efp_potential
   use mqc_czt_efp_read, only: efp_fragment_t, read_efp_potential
   use mqc_czt_efp_convert, only: potential_to_fragment
   use mqc_czt_efp_energy, only: efp_energy_t, efp_interaction_energy, &
                                 efp_pair_energy_t, efp_pair_terms, &
                                 pair_polarization_energy
   use mqc_czt_efmo_pairs, only: efmo_pair_distance, efmo_split_pairs, &
                                 vdw_scaled_distance
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_atomic_guess, only: build_restricted_guess
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_scf_types, only: scf_numerics_t
   use mqc_diis, only: ACCEL_DIIS
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_efmo_pieces_tests

   real(dp), parameter :: ANG = ANGSTROM_TO_BOHR
      !! Angstrom to Bohr, so the geometries below can be written the way they
      !! are measured. **The code's own constant, not a literal**: the cutoff
      !! test places a pair at a stated `R_IJ`, and `unitless_distance` divides
      !! by radii converted with this one, so a spelling that differed in the
      !! tenth digit would put the hand calculation out by that much.

   character(len=*), parameter :: BASIS = "6-31g"
      !! Small on purpose. Nothing here is compared against an external number,
      !! so what the basis has to be is cheap and Cartesian-clean.

   !! The potential and the fragment it converts to, built once. Every test
   !! wants one, and they are identical every time.
   type(efp_potential_t), save :: cached_pot
   type(efp_fragment_t), save :: cached_frag
   logical, save :: cached_ready = .false.
   real(dp), save :: cached_scf_energy = 0.0_dp

contains

   subroutine collect_mqc_czt_efmo_pieces_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efmo_makefp_hands_back_its_scf", test_makefp_scf), &
                  new_unittest("efmo_potential_to_fragment_round_trip", test_convert), &
                  new_unittest("efmo_pair_polarization_is_the_total_for_two", test_pair_pol), &
                  new_unittest("efmo_three_body_induction_vanishes_with_distance", &
                               test_three_body), &
                  new_unittest("efmo_pair_distance_and_split", test_split), &
                  new_unittest("efmo_pair_terms_sum_to_the_all_pairs_energy", test_pair_terms) &
                  ]
   end subroutine collect_mqc_czt_efmo_pieces_tests

   subroutine water_geometry(z, symbols, coords)
      !! One water, in Bohr, in the yz plane
      integer, intent(out) :: z(3)
      character(len=2), intent(out) :: symbols(3)
      real(dp), intent(out) :: coords(3, 3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      coords = reshape([0.00000000000000_dp, 0.00000000009155_dp, 0.10077199490609_dp, &
                        0.00000000000000_dp, 0.77250895271063_dp, -0.46780199741728_dp, &
                        0.00000000000000_dp, -0.77250895280218_dp, -0.46780199748881_dp], &
                       [3, 3])*ANG
   end subroutine water_geometry

   subroutine water_potential(pot, err)
      !! The cached MAKEFP potential for that water
      type(efp_potential_t), intent(out) :: pot
      type(error_t), intent(inout) :: err

      integer :: z(3)
      character(len=2) :: symbols(3)
      real(dp) :: coords(3, 3)
      type(rhf_result_t) :: scf

      if (.not. cached_ready) then
         call water_geometry(z, symbols, coords)
         call make_efp_potential(z, symbols, coords, BASIS, "WATER", cached_pot, err, &
                                 scf_out=scf)
         if (err%has_error()) return
         cached_scf_energy = scf%energy
         call potential_to_fragment(cached_pot, cached_frag, err)
         if (err%has_error()) return
         cached_ready = .true.
      end if
      pot = cached_pot
   end subroutine water_potential

   subroutine water_fragment(frag, err)
      !! The same water as a fragment, converted in memory
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: err

      type(efp_potential_t) :: pot

      call water_potential(pot, err)
      if (err%has_error()) return
      frag = cached_frag
   end subroutine water_fragment

   subroutine test_makefp_scf(error)
      !! `E_I^0` from MAKEFP is the monomer's own RHF energy
      !!
      !! **Piece one of Phase 1.** MAKEFP already runs a closed-shell SCF over
      !! the monomer; EFMO's fragment sum needs exactly that energy. Before
      !! `scf_out` existed the only way to have it was to run the SCF again,
      !! which for a real cluster is a second SCF per fragment -- and worse, an
      !! SCF that could silently be a *different* one, converged elsewhere or in
      !! another angular form, so the potential and the energy would describe
      !! two different determinants.
      !!
      !! So this reproduces the SCF MAKEFP performs, argument for argument, and
      !! requires the two energies to agree to 1e-10 -- the energy threshold the
      !! SCF is converged to, which is as tight as the comparison can be.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(scf_numerics_t) :: settings
      real(dp), allocatable :: guess_density(:, :)
      integer :: z(3), guess_kind
      character(len=2) :: symbols(3)
      real(dp) :: coords(3, 3)

      call water_potential(pot, err)
      call check(error,.not. err%has_error(), "building the potential failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      ! The same molecule MAKEFP builds: Cartesian, because a `.efp` is read by
      ! GAMESS. Spherical would be a different energy, and a comparison passing
      ! at 1e-10 is what says this is the same determinant.
      call water_geometry(z, symbols, coords)
      call build_czt_molecule(z, symbols, coords, BASIS, mol, err, force_cartesian=.true.)
      call check(error,.not. err%has_error(), "building the molecule failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call build_restricted_guess(mol, "auto", guess_kind, guess_density, err)
      call check(error,.not. err%has_error(), "building the guess failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      ! MAKEFP's own defaults, listed at its SCF call: 1e-10 energy, 1e-8
      ! density, the commutator following the density, 200 iterations, DIIS over
      ! eight vectors, no level shift, incremental Fock building.
      settings = scf_numerics_t()
      settings%max_iter = 200
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp
      settings%grad_tol = 1.0e-8_dp
      settings%diis_size = 8
      settings%level_shift = 0.0_dp
      settings%linear_dependence = 0.0_dp
      settings%incremental_fock = .true.
      call run_czt_rhf(mol, 10, 200, 1.0e-10_dp, 1.0e-8_dp, .false., scf, err, &
                       guess=guess_kind, guess_density=guess_density, &
                       grad_tol=1.0e-8_dp, diis_vectors=8, level_shift=0.0_dp, &
                       linear_dependence=0.0_dp, accelerator=ACCEL_DIIS, &
                       incremental_fock=.true., scf=settings)
      call check(error,.not. err%has_error(), "the reference SCF failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, scf%converged, "the reference SCF did not converge")
      if (allocated(error)) return

      call check(error, cached_scf_energy, scf%energy, thr=1.0e-10_dp, &
                 message="MAKEFP's SCF energy is not the monomer's RHF energy")
      if (allocated(error)) return
      ! And that the potential carries the same number, since an orchestrator
      ! holding potentials rather than SCF results reads it from there.
      call check(error, pot%scf_energy, scf%energy, thr=1.0e-10_dp, &
                 message="the potential's scf_energy is not the monomer's RHF energy")
      call mol%destroy()
   end subroutine test_makefp_scf

   subroutine test_convert(error)
      !! `potential_to_fragment` equals write-then-read, field by field
      !!
      !! **Piece two of Phase 1**, and the one that can be wrong most quietly:
      !! every consumer of an `efp_fragment_t` today was written against what
      !! the *reader* produces, so a converter that packs a tensor in another
      !! order or leaves a `has_*` flag down does not fail -- it computes a term
      !! from data in a layout the term does not expect, or drops the term.
      !!
      !! Compared at the file's own precision rather than exactly: the round
      !! trip goes through ten printed decimals, and the orbital blocks through
      !! eight significant figures, so those are what the comparison can ask
      !! for. What is compared exactly is every count and every flag.
      type(error_type), allocatable, intent(out) :: error

      type(efp_potential_t) :: pot
      type(efp_fragment_t) :: written, converted
      type(error_t) :: err
      character(len=*), parameter :: path = "test_efmo_convert.efp"
      real(dp), parameter :: FILE_TOL = 1.0e-8_dp
         !! What ten printed decimals leave, with room for the exponent-form
         !! blocks.
      real(dp), parameter :: ORBITAL_TOL = 1.0e-7_dp
         !! The wavefunction blocks print `ES15.8E2`, so eight significant
         !! figures on coefficients of order one.
      integer :: unit, stat

      call water_potential(pot, err)
      call check(error,.not. err%has_error(), "building the potential failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call write_efp_potential(pot, path, err)
      call check(error,.not. err%has_error(), "writing the potential failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call read_efp_potential(path, written, err)
      open (newunit=unit, file=path, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
      call check(error,.not. err%has_error(), "reading the potential back failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call potential_to_fragment(pot, converted, err)
      call check(error,.not. err%has_error(), "converting the potential failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      ! Nothing below means anything if the potential carries no blocks, and a
      ! comparison of two absent arrays passes. So: this potential has all of
      ! them.
      call check(error, written%has_dynamic .and. written%has_dipquad &
                 .and. written%has_quadquad .and. written%has_static_pol &
                 .and. written%has_lmo .and. written%has_fock .and. written%has_ctvec &
                 .and. written%has_ctfok .and. written%has_basis &
                 .and. written%has_screen .and. written%has_screen2, &
                 "the round trip produced a fragment missing blocks, so the "// &
                 "comparison below would be vacuous")
      if (allocated(error)) return

      ! --- the counts and the name ----------------------------------------------
      call check(error, converted%name == written%name, "the fragment name differs")
      if (allocated(error)) return
      call check(error, converted%n_points, written%n_points, message="n_points")
      if (allocated(error)) return
      call check(error, converted%n_atoms, written%n_atoms, message="n_atoms")
      if (allocated(error)) return
      call check(error, converted%multiplicity, written%multiplicity, message="multiplicity")
      if (allocated(error)) return
      call check(error, converted%n_lmo, written%n_lmo, message="n_lmo")
      if (allocated(error)) return
      call check(error, converted%n_freq, written%n_freq, message="n_freq")
      if (allocated(error)) return
      call check(error, converted%n_pol, written%n_pol, message="n_pol")
      if (allocated(error)) return
      call check(error, converted%n_dipquad, written%n_dipquad, message="n_dipquad")
      if (allocated(error)) return
      call check(error, converted%n_quadquad, written%n_quadquad, message="n_quadquad")
      if (allocated(error)) return
      call check(error, converted%n_shells, written%n_shells, message="n_shells")
      if (allocated(error)) return
      call check(error, converted%n_lmo_proj, written%n_lmo_proj, message="n_lmo_proj")
      if (allocated(error)) return
      call check(error, converted%nao_proj, written%nao_proj, message="nao_proj")
      if (allocated(error)) return
      call check(error, converted%n_occ_ct, written%n_occ_ct, message="n_occ_ct")
      if (allocated(error)) return
      call check(error, converted%n_mo_ct, written%n_mo_ct, message="n_mo_ct")
      if (allocated(error)) return

      ! --- every flag, exactly ---------------------------------------------------
      call check(error, converted%has_screen .eqv. written%has_screen, "has_screen")
      if (allocated(error)) return
      call check(error, converted%has_screen2 .eqv. written%has_screen2, "has_screen2")
      if (allocated(error)) return
      call check(error, converted%has_dynamic .eqv. written%has_dynamic, "has_dynamic")
      if (allocated(error)) return
      call check(error, converted%has_static_pol .eqv. written%has_static_pol, &
                 "has_static_pol")
      if (allocated(error)) return
      call check(error, converted%has_dipquad .eqv. written%has_dipquad, "has_dipquad")
      if (allocated(error)) return
      call check(error, converted%has_quadquad .eqv. written%has_quadquad, "has_quadquad")
      if (allocated(error)) return
      call check(error, converted%has_basis .eqv. written%has_basis, "has_basis")
      if (allocated(error)) return
      call check(error, converted%has_lmo .eqv. written%has_lmo, "has_lmo")
      if (allocated(error)) return
      call check(error, converted%has_fock .eqv. written%has_fock, "has_fock")
      if (allocated(error)) return
      call check(error, converted%has_ctvec .eqv. written%has_ctvec, "has_ctvec")
      if (allocated(error)) return
      call check(error, converted%has_ctfok .eqv. written%has_ctfok, "has_ctfok")
      if (allocated(error)) return

      ! --- every allocated block --------------------------------------------------
      call check(error, all(converted%labels == written%labels), "the point labels differ")
      if (allocated(error)) return
      call compare2(error, converted%points, written%points, FILE_TOL, "points")
      if (allocated(error)) return
      call compare(error, converted%mass, written%mass, 1.0e-6_dp, "mass")
      if (allocated(error)) return
      call compare(error, converted%charge, written%charge, FILE_TOL, "charge")
      if (allocated(error)) return
      call compare(error, converted%q_elec, written%q_elec, FILE_TOL, "q_elec")
      if (allocated(error)) return
      call compare(error, converted%q_nuc, written%q_nuc, FILE_TOL, "q_nuc")
      if (allocated(error)) return
      call compare2(error, converted%dipole, written%dipole, FILE_TOL, "dipole")
      if (allocated(error)) return
      call compare2(error, converted%quadrupole, written%quadrupole, FILE_TOL, "quadrupole")
      if (allocated(error)) return
      call compare2(error, converted%octopole, written%octopole, FILE_TOL, "octopole")
      if (allocated(error)) return
      call compare(error, converted%screen, written%screen, FILE_TOL, "screen")
      if (allocated(error)) return
      call compare(error, converted%screen2, written%screen2, FILE_TOL, "screen2")
      if (allocated(error)) return
      call compare(error, reshape(converted%static_pol, [size(converted%static_pol)]), &
                   reshape(written%static_pol, [size(written%static_pol)]), FILE_TOL, &
                   "static_pol")
      if (allocated(error)) return
      call compare2(error, converted%pol_points, written%pol_points, FILE_TOL, "pol_points")
      if (allocated(error)) return
      call compare(error, reshape(converted%dyn_pol, [size(converted%dyn_pol)]), &
                   reshape(written%dyn_pol, [size(written%dyn_pol)]), FILE_TOL, "dyn_pol")
      if (allocated(error)) return
      call compare2(error, converted%centroids, written%centroids, FILE_TOL, "centroids")
      if (allocated(error)) return
      call compare(error, converted%frequencies, written%frequencies, 1.0e-6_dp, &
                   "frequencies")
      if (allocated(error)) return
      call compare(error, reshape(converted%dipquad, [size(converted%dipquad)]), &
                   reshape(written%dipquad, [size(written%dipquad)]), FILE_TOL, "dipquad")
      if (allocated(error)) return
      call compare(error, reshape(converted%quadquad, [size(converted%quadquad)]), &
                   reshape(written%quadquad, [size(written%quadquad)]), FILE_TOL, "quadquad")
      if (allocated(error)) return
      call compare2(error, converted%lmo_gamess, written%lmo_gamess, ORBITAL_TOL, &
                    "lmo_gamess")
      if (allocated(error)) return
      call compare2(error, converted%fock_lmo, written%fock_lmo, FILE_TOL, "fock_lmo")
      if (allocated(error)) return
      call compare2(error, converted%ctvec_gamess, written%ctvec_gamess, ORBITAL_TOL, &
                    "ctvec_gamess")
      if (allocated(error)) return
      call compare(error, converted%eps_occ, written%eps_occ, FILE_TOL, "eps_occ")
      if (allocated(error)) return

      ! The projection basis goes through the reader in both cases, so this is
      ! exact -- and what it guards is that the converter hands the reader the
      ! *same lines*, `L` shells and atom headers included.
      call check(error, all(converted%shell_atom == written%shell_atom), "shell_atom")
      if (allocated(error)) return
      call check(error, all(converted%shell_l == written%shell_l), "shell_l")
      if (allocated(error)) return
      call check(error, all(converted%shell_first == written%shell_first), "shell_first")
      if (allocated(error)) return
      call check(error, all(converted%shell_nprim == written%shell_nprim), "shell_nprim")
      if (allocated(error)) return
      call compare(error, converted%prim_expo, written%prim_expo, 0.0_dp, "prim_expo")
      if (allocated(error)) return
      call compare(error, converted%prim_coef, written%prim_coef, 0.0_dp, "prim_coef")
      if (allocated(error)) return

      ! And that the converted fragment computes the same energy the read one
      ! does, which is the thing every consumer actually asks of it.
      call check(error, energy_of_pair(converted, err), energy_of_pair(written, err), &
                 thr=1.0e-8_dp, message="the two fragments do not interact alike")
      call written%destroy()
      call converted%destroy()
   end subroutine test_convert

   function energy_of_pair(frag, err) result(total)
      !! Two copies of one fragment, three Angstrom apart, as one number
      type(efp_fragment_t), intent(in) :: frag
      type(error_t), intent(inout) :: err
      real(dp) :: total

      type(efp_fragment_t) :: two(2)
      type(efp_energy_t) :: e
      real(dp) :: shifts(3, 2)

      two(1) = frag
      two(2) = frag
      shifts = 0.0_dp
      shifts(1, 2) = 3.0_dp*ANG
      e = efp_interaction_energy(two, shifts, err)
      total = e%total
   end function energy_of_pair

   subroutine test_pair_pol(error)
      !! For two fragments the pair induction *is* the total induction
      !!
      !! **Piece three of Phase 1.** `E_IJ^pol` is subtracted from every
      !! quantum-mechanical dimer so that `E_pol^total` does not count that pair's
      !! induction twice, and the subtraction is only clean if the two are the
      !! same quantity: the same static field, the same screening, the same
      !! iteration to the same tolerance. On a system of exactly two fragments
      !! that identity is checkable outright, and at 1e-12 -- the tolerance the
      !! induced dipoles are solved to -- because it is not an approximation but
      !! the same solver on the same input.
      !!
      !! Symmetry in `I` and `J` is asserted too: the induced dipoles are solved
      !! together, so a pair term that depended on which fragment was listed
      !! first would mean the solve is not symmetric in them.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: frags(2)
      type(error_t) :: err
      type(efp_energy_t) :: e
      real(dp) :: shifts(3, 2), pair, swapped

      call water_fragment(frags(1), err)
      call water_fragment(frags(2), err)
      call check(error,.not. err%has_error(), "building the fragments failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      shifts = 0.0_dp
      shifts(1, 2) = 3.0_dp*ANG
      e = efp_interaction_energy(frags, shifts, err)
      call check(error,.not. err%has_error(), "the interaction energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      pair = pair_polarization_energy(frags(1), frags(2), shifts(:, 1), shifts(:, 2), err)
      call check(error,.not. err%has_error(), "the pair induction failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, pair, e%polarization, thr=1.0e-12_dp, &
                 message="the pair induction is not the two-fragment total")
      if (allocated(error)) return

      swapped = pair_polarization_energy(frags(2), frags(1), shifts(:, 2), shifts(:, 1), err)
      call check(error,.not. err%has_error(), "the swapped pair induction failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, swapped, pair, thr=1.0e-14_dp, &
                 message="the pair induction is not symmetric in I and J")
   end subroutine test_pair_pol

   subroutine test_three_body(error)
      !! What the pair subtraction leaves behind is three-body induction
      !!
      !! With three fragments, `E_pol^total` minus the three pair terms is the
      !! part of the induction that no pair can hold: the field one fragment's
      !! induced dipoles add at a third. **This is the term EFMO keeps**, and it
      !! is the reason the pair subtraction exists rather than simply not
      !! solving the total.
      !!
      !! Two things are asserted, and the second is the one with teeth. That the
      !! remainder is a correction to the pair sum rather than larger than it,
      !! which says the pair terms are the same quantity the total is built
      !! from; and that pulling the third water away drives it to zero -- an induction that is genuinely
      !! three-body must vanish when one body is not there, where a
      !! bookkeeping mistake in the subtraction would leave a residue of the
      !! remaining pair's own induction.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: frags(3)
      type(error_t) :: err
      type(efp_energy_t) :: e
      real(dp) :: shifts(3, 3), pairs, remainder, remote
      integer :: k

      do k = 1, 3
         call water_fragment(frags(k), err)
      end do
      call check(error,.not. err%has_error(), "building the fragments failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      ! A compact triangle, so the three-body term is not numerically zero.
      shifts = 0.0_dp
      shifts(1, 2) = 4.0_dp*ANG
      shifts(1, 3) = 2.0_dp*ANG
      shifts(2, 3) = 3.5_dp*ANG

      e = efp_interaction_energy(frags, shifts, err)
      call check(error,.not. err%has_error(), "the interaction energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      pairs = three_pair_sum(frags, shifts, err)
      call check(error,.not. err%has_error(), "the pair inductions failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      remainder = e%polarization - pairs

      ! **Induction is strongly non-additive**, and this is where that shows:
      ! at this geometry the three pairs carry -9.63e-5 Ha of the total
      ! -1.72e-4, leaving -7.59e-5 -- 44% of the total -- that no pair holds.
      ! The energy is quadratic in the field each point sees, so the square of a
      ! sum of two fields keeps a cross term that neither pair has. That is
      ! precisely why EFMO solves one induction over every fragment and
      ! subtracts pair terms, rather than summing pair inductions.
      !
      ! So what is asserted is the ordering the subtraction needs -- the
      ! remainder is a correction to the pair sum rather than larger than it --
      ! and that it is not numerically zero, without which the vanishing test
      ! below would show nothing.
      call check(error, abs(remainder) < abs(pairs), &
                 "what the pair subtraction leaves is larger than the pair sum "// &
                 "itself, so the pair terms are not the same quantity the total is "// &
                 "built from")
      if (allocated(error)) return
      call check(error, abs(remainder) > 1.0e-9_dp, &
                 "this geometry produces no three-body induction at all, so the "// &
                 "test below cannot show it vanishing")
      if (allocated(error)) return

      ! The third water taken far away. Induction falls off as a high power of
      ! the separation, so at 60 Angstrom the three-body part is gone while the
      ! remaining pair's is untouched.
      shifts(:, 3) = [60.0_dp*ANG, 0.0_dp, 0.0_dp]
      e = efp_interaction_energy(frags, shifts, err)
      call check(error,.not. err%has_error(), "the separated energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      pairs = three_pair_sum(frags, shifts, err)
      call check(error,.not. err%has_error(), "the separated pair inductions failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      remote = e%polarization - pairs
      ! Measured 1.7e-8 Ha, against 7.6e-5 with the third water in the cluster:
      ! four orders down, and what is left is the tail of a field falling as the
      ! cube of a separation that is only twenty times larger.
      call check(error, abs(remote) < 1.0e-6_dp, &
                 "the three-body induction does not vanish when the third fragment "// &
                 "is taken away, so what the subtraction leaves is not three-body")
      if (allocated(error)) return
      call check(error, abs(remote) < 1.0e-3_dp*abs(remainder), &
                 "moving the third fragment away barely changed the three-body term")
   end subroutine test_three_body

   function three_pair_sum(frags, shifts, err) result(total)
      !! The three pair inductions of a three-fragment system
      type(efp_fragment_t), intent(in) :: frags(:)
      real(dp), intent(in) :: shifts(:, :)
      type(error_t), intent(inout) :: err
      real(dp) :: total

      integer :: a, b

      total = 0.0_dp
      do a = 1, size(frags) - 1
         do b = a + 1, size(frags)
            total = total + pair_polarization_energy(frags(a), frags(b), &
                                                     shifts(:, a), shifts(:, b), err)
            if (err%has_error()) return
         end do
      end do
   end function three_pair_sum

   subroutine test_split(error)
      !! `R_IJ` and the cutoff it decides, at a geometry worked out by hand
      !!
      !! **Piece four of Phase 1.** Two identical waters, both lying in the `yz`
      !! plane, translated along `x` by `d`. Every atom pair is then at least `d`
      !! apart and the oxygen pair is at exactly `d`, so the pair with the
      !! largest radii sum -- oxygen against oxygen, `2 x 1.52` Angstrom in the
      !! Bondi set this code carries -- is the minimum of eq 2:
      !!
      !!     R_IJ = d / (2 * 1.52 Angstrom)
      !!
      !! The two hydrogens sit `2 x 1.10` apart and the O-H pairs `2.62`, both
      !! at distances no shorter than `d`, so neither can undercut it. Placing
      !! `d` at `1.9` and `2.1` of that unit therefore puts the pair on either
      !! side of the default `R_cut = 2.0`, and which list it lands in is the
      !! whole QM/EFP decision of eq 6.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: O_VDW = 1.40_dp, H_VDW = 1.20_dp
         !! GAMESS's `$FMO VDWRAD` table, as `mqc_atomic_radii` carries it.
         !! Written out so the geometry below is derived here rather than read
         !! back from the code. Bondi's 1.52 and 1.10 are *not* what the cutoff
         !! is measured in -- the FMO literature quotes `R_cut` against this
         !! table, and the two disagree by six per cent on a water pair.
      real(dp), parameter :: RCUT = 2.0_dp
      type(error_t) :: err
      integer :: z(3), owner(6), z_pair(6)
      character(len=2) :: symbols(3)
      real(dp) :: coords(3, 3), xyz(3, 6), d, r
      integer, allocatable :: qm(:, :), efp(:, :)

      call water_geometry(z, symbols, coords)

      ! Inside the cutoff: R_IJ = 1.9.
      d = 1.9_dp*(2.0_dp*O_VDW)*ANG
      call two_waters(coords, z, d, xyz, z_pair, owner)
      call efmo_pair_distance(z_pair(1:3), xyz(:, 1:3), z_pair(4:6), xyz(:, 4:6), r, err)
      call check(error,.not. err%has_error(), "the pair distance failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, r, 1.9_dp, thr=1.0e-12_dp, &
                 message="R_IJ is not the oxygen pair over twice the oxygen radius")
      if (allocated(error)) return
      ! And that the scalar kernel agrees on that one atom pair, since the
      ! minimum is meant to be over pairs of exactly this quantity.
      call check(error, vdw_scaled_distance(8, xyz(:, 1), 8, xyz(:, 4)), 1.9_dp, &
                 thr=1.0e-12_dp, message="the scaled distance of the oxygen pair")
      if (allocated(error)) return

      call efmo_split_pairs(owner, z_pair, xyz, RCUT, qm, efp, err)
      call check(error,.not. err%has_error(), "the split failed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, size(qm, 2), 1, message="a pair inside the cutoff is not QM")
      if (allocated(error)) return
      call check(error, size(efp, 2), 0, message="a pair inside the cutoff is also EFP")
      if (allocated(error)) return
      call check(error, qm(1, 1) == 1 .and. qm(2, 1) == 2, "the QM pair is not (1, 2)")
      if (allocated(error)) return

      ! Outside it: R_IJ = 2.1, the same geometry pushed apart.
      d = 2.1_dp*(2.0_dp*O_VDW)*ANG
      call two_waters(coords, z, d, xyz, z_pair, owner)
      call efmo_pair_distance(z_pair(1:3), xyz(:, 1:3), z_pair(4:6), xyz(:, 4:6), r, err)
      call check(error,.not. err%has_error(), "the far pair distance failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, r, 2.1_dp, thr=1.0e-12_dp, message="R_IJ of the far pair")
      if (allocated(error)) return
      ! The hydrogen pair is the one that could undercut it -- the smallest
      ! radii sum in the system -- so its value is asserted rather than assumed.
      call check(error, vdw_scaled_distance(1, xyz(:, 2), 1, xyz(:, 5)) &
                 > 2.1_dp, "the hydrogen pair, not the oxygen pair, is the minimum "// &
                 "at this geometry, so the hand calculation above does not hold")
      if (allocated(error)) return
      call check(error, 2.0_dp*O_VDW > 2.0_dp*H_VDW, "the radii are not what this "// &
                 "test assumes")
      if (allocated(error)) return

      call efmo_split_pairs(owner, z_pair, xyz, RCUT, qm, efp, err)
      call check(error,.not. err%has_error(), "the far split failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, size(qm, 2), 0, message="a pair beyond the cutoff is still QM")
      if (allocated(error)) return
      call check(error, size(efp, 2), 1, message="a pair beyond the cutoff is not EFP")
      if (allocated(error)) return
      call check(error, efp(1, 1) == 1 .and. efp(2, 1) == 2, "the EFP pair is not (1, 2)")
      if (allocated(error)) return

      ! Three fragments: every pair in exactly one list, whatever the cutoff.
      call check(error, every_pair_once(1.0_dp), "with R_cut = 1.0 the three pairs "// &
                 "are not covered exactly once")
      if (allocated(error)) return
      call check(error, every_pair_once(2.0_dp), "with R_cut = 2.0 the three pairs "// &
                 "are not covered exactly once")
      if (allocated(error)) return
      call check(error, every_pair_once(1.0e6_dp), "with a huge R_cut the three pairs "// &
                 "are not covered exactly once")
   end subroutine test_split

   subroutine two_waters(coords, z, d, xyz, z_pair, owner)
      !! Two copies of one water, the second translated along `x` by `d`
      real(dp), intent(in) :: coords(3, 3)
      integer, intent(in) :: z(3)
      real(dp), intent(in) :: d
      real(dp), intent(out) :: xyz(3, 6)
      integer, intent(out) :: z_pair(6), owner(6)

      integer :: k

      do k = 1, 3
         xyz(:, k) = coords(:, k)
         xyz(:, k + 3) = coords(:, k) + [d, 0.0_dp, 0.0_dp]
         z_pair(k) = z(k)
         z_pair(k + 3) = z(k)
         owner(k) = 1
         owner(k + 3) = 2
      end do
   end subroutine two_waters

   function every_pair_once(rcut) result(ok)
      !! Three fragments, and each of the three pairs in exactly one list
      real(dp), intent(in) :: rcut
      logical :: ok

      type(error_t) :: err
      integer :: z(3), owner(9), z_all(9), seen(3, 3), k, p
      character(len=2) :: symbols(3)
      real(dp) :: coords(3, 3), xyz(3, 9)
      integer, allocatable :: qm(:, :), efp(:, :)

      call water_geometry(z, symbols, coords)
      do k = 1, 3
         xyz(:, k) = coords(:, k)
         xyz(:, k + 3) = coords(:, k) + [5.0_dp*ANG, 0.0_dp, 0.0_dp]
         xyz(:, k + 6) = coords(:, k) + [0.0_dp, 14.0_dp*ANG, 0.0_dp]
         z_all(k) = z(k)
         z_all(k + 3) = z(k)
         z_all(k + 6) = z(k)
         owner(k) = 1
         owner(k + 3) = 2
         owner(k + 6) = 3
      end do

      call efmo_split_pairs(owner, z_all, xyz, rcut, qm, efp, err)
      ok = .not. err%has_error()
      if (.not. ok) return
      seen = 0
      do p = 1, size(qm, 2)
         seen(qm(1, p), qm(2, p)) = seen(qm(1, p), qm(2, p)) + 1
      end do
      do p = 1, size(efp, 2)
         seen(efp(1, p), efp(2, p)) = seen(efp(1, p), efp(2, p)) + 1
      end do
      ok = seen(1, 2) == 1 .and. seen(1, 3) == 1 .and. seen(2, 3) == 1 &
           .and. size(qm, 2) + size(efp, 2) == 3
   end function every_pair_once

   subroutine test_pair_terms(error)
      !! The per-pair EFP terms sum to what the all-pairs energy reports
      !!
      !! **Piece five of Phase 1.** `efp_interaction_energy` computes one number
      !! per term over a whole system; EFMO needs the four non-induction terms
      !! for a *chosen list* of pairs, because a near pair gets quantum
      !! mechanics instead and must contribute none of them.
      !!
      !! Over the all-pairs list the two must agree, and that is a real check
      !! rather than a tautology for electrostatics: the system-wide routine
      !! flattens every fragment's points into one set and screens pairwise
      !! within it, where the per-pair route builds a two-fragment system each
      !! time. Exchange repulsion, dispersion and charge transfer go through the
      !! same pair routines in both, so what they check is the bookkeeping.
      !!
      !! Symmetry under swapping `I` and `J` is asserted separately: every term
      !! is a property of the unordered pair, and the pair routines take the two
      !! fragments in different roles.
      type(error_type), allocatable, intent(out) :: error

      type(efp_fragment_t) :: frags(3)
      type(error_t) :: err
      type(efp_energy_t) :: e
      type(efp_pair_energy_t), allocatable :: terms(:), swapped(:)
      real(dp) :: shifts(3, 3)
      integer :: pairs(2, 3), back(2, 3), k
      real(dp) :: coul, disp, exrep, ct

      do k = 1, 3
         call water_fragment(frags(k), err)
      end do
      call check(error,.not. err%has_error(), "building the fragments failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      shifts = 0.0_dp
      shifts(1, 2) = 3.0_dp*ANG
      shifts(1, 3) = 1.5_dp*ANG
      shifts(2, 3) = 2.8_dp*ANG

      e = efp_interaction_energy(frags, shifts, err)
      call check(error,.not. err%has_error(), "the interaction energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      pairs = reshape([1, 2, 1, 3, 2, 3], [2, 3])
      terms = efp_pair_terms(frags, shifts, pairs, err)
      call check(error,.not. err%has_error(), "the pair terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      coul = sum(terms%electrostatics)
      disp = sum(terms%dispersion)
      exrep = sum(terms%exchange_repulsion)
      ct = sum(terms%charge_transfer)
      call check(error, coul, e%electrostatics, thr=1.0e-10_dp, &
                 message="the per-pair electrostatics do not sum to the system's")
      if (allocated(error)) return
      call check(error, exrep, e%exchange_repulsion, thr=1.0e-10_dp, &
                 message="the per-pair exchange repulsion does not sum to the system's")
      if (allocated(error)) return
      call check(error, disp, e%dispersion, thr=1.0e-10_dp, &
                 message="the per-pair dispersion does not sum to the system's")
      if (allocated(error)) return
      call check(error, sum(terms%dispersion_e7), e%dispersion_e7, thr=1.0e-10_dp, &
                 message="the per-pair E7 does not sum to the system's")
      if (allocated(error)) return
      call check(error, ct, e%charge_transfer, thr=1.0e-10_dp, &
                 message="the per-pair charge transfer does not sum to the system's")
      if (allocated(error)) return
      ! No polarization anywhere in the pair terms: it is many-body and belongs
      ! to `E_pol^total` alone.
      call check(error, sum(terms%total), coul + disp + exrep + ct, thr=1.0e-14_dp, &
                 message="a pair total is not the sum of its four terms")
      if (allocated(error)) return

      ! Charge transfer switched off, since the original method leaves it out.
      terms = efp_pair_terms(frags, shifts, pairs, err, charge_transfer_on=.false.)
      call check(error,.not. err%has_error(), "the CT-free pair terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      call check(error, sum(terms%charge_transfer), 0.0_dp, thr=0.0_dp, &
                 message="charge transfer was switched off and still contributed")
      if (allocated(error)) return

      ! Every pair the other way round.
      back = reshape([2, 1, 3, 1, 3, 2], [2, 3])
      terms = efp_pair_terms(frags, shifts, pairs, err)
      swapped = efp_pair_terms(frags, shifts, back, err)
      call check(error,.not. err%has_error(), "the swapped pair terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return
      do k = 1, 3
         call check(error, swapped(k)%electrostatics, terms(k)%electrostatics, &
                    thr=1.0e-12_dp, message="electrostatics is not pair-symmetric")
         if (allocated(error)) return
         call check(error, swapped(k)%exchange_repulsion, terms(k)%exchange_repulsion, &
                    thr=1.0e-12_dp, message="exchange repulsion is not pair-symmetric")
         if (allocated(error)) return
         call check(error, swapped(k)%dispersion, terms(k)%dispersion, thr=1.0e-12_dp, &
                    message="dispersion is not pair-symmetric")
         if (allocated(error)) return
         call check(error, swapped(k)%charge_transfer, terms(k)%charge_transfer, &
                    thr=1.0e-12_dp, message="charge transfer is not pair-symmetric")
         if (allocated(error)) return
      end do
   end subroutine test_pair_terms

   subroutine compare(error, got, want, tol, what)
      !! Two rank-one blocks, by their largest disagreement
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: got(:), want(:)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: what

      call check(error, size(got), size(want), message=what//": different sizes")
      if (allocated(error)) return
      if (size(got) == 0) return
      call check(error, maxval(abs(got - want)), 0.0_dp, thr=tol, &
                 message=what//" differs between the converted and the read fragment")
   end subroutine compare

   subroutine compare2(error, got, want, tol, what)
      !! The same for a rank-two block
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: got(:, :), want(:, :)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: what

      call check(error, size(got, 1), size(want, 1), message=what//": different rows")
      if (allocated(error)) return
      call check(error, size(got, 2), size(want, 2), message=what//": different columns")
      if (allocated(error)) return
      if (size(got) == 0) return
      call check(error, maxval(abs(got - want)), 0.0_dp, thr=tol, &
                 message=what//" differs between the converted and the read fragment")
   end subroutine compare2

end module test_mqc_czt_efmo_pieces

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efmo_pieces, only: collect_mqc_czt_efmo_pieces_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efmo_pieces", collect_mqc_czt_efmo_pieces_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
