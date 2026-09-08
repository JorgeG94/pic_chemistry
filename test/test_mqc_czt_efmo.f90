!! The EFMO orchestrator, at the four limits where its total is known
module test_mqc_czt_efmo
   !! `run_efmo` assembles eq 6 of the EFMO paper (Sattasathuchana et al., JCTC
   !! 20, 2445 (2024)):
   !!
   !!     E = sum_I E_I^0
   !!       + sum_{R_IJ <= R_cut} (E_IJ^0 - E_I^0 - E_J^0 - E_IJ^pol)
   !!       + sum_{R_IJ >  R_cut} (E_IJ^Coul + E_IJ^disp + E_IJ^ExRep + E_IJ^CT)
   !!       + E_pol^total
   !!
   !! No reference for an EFMO energy exists yet -- Phase 3 gets those from
   !! GAMESS -- so what is checked here is not a number but the **identities the
   !! expression has by construction**, each of which fails if a term is
   !! dropped, double counted, or assembled with the wrong sign:
   !!
   !!   1. `R_cut` huge: every pair is quantum, so the expression collapses to
   !!      FMO2 in vacuo plus the *many-body* part of the induction.
   !!   2. `R_cut` zero: no pair is quantum, so it collapses to the in-vacuo
   !!      monomers plus the whole EFP-EFP interaction energy.
   !!   3. Two fragments, both quantum: their pair induction *is* the total
   !!      induction, the two cancel exactly, and every other term telescopes
   !!      away -- so the answer is the dimer's own RHF energy.
   !!   4. `R_cut` between the pair separations of a trimer: one quantum dimer
   !!      and two effective ones, against a total assembled here by hand from
   !!      the Phase 1 pieces.
   !!
   !! **Limit 1 is not "EFMO is nearly FMO2".** The induction is strongly
   !! non-additive -- three waters at four Angstrom carry 44 per cent of their
   !! total induction in terms no pair has -- so the difference in test one is a
   !! term of the method, of the same size as the pair energies. Only the
   !! identity is asserted, and the measured difference is reported by the run
   !! itself.
   !!
   !! Water in 6-31G, one trimer geometry for every case, and the reference
   !! potentials built once: MAKEFP is the whole cost of this file, and
   !! `run_efmo` builds its own set on each call.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_efmo, only: efmo_options_t, efmo_result_t, run_efmo
   use mqc_czt_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential
   use mqc_czt_efp_read, only: efp_fragment_t
   use mqc_czt_efp_convert, only: potential_to_fragment
   use mqc_czt_efp_energy, only: efp_energy_t, efp_interaction_energy, &
                                 efp_pair_energy_t, efp_pair_terms, &
                                 pair_polarization_energy
   use mqc_czt_efp_interaction, only: efp_system_t, build_efp_system, polarization_energy
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_atomic_guess, only: build_restricted_guess
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_scf_types, only: scf_numerics_t
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_mqc_czt_efmo_tests

   real(dp), parameter :: ANG = ANGSTROM_TO_BOHR

   character(len=*), parameter :: BASIS = "6-31g"
      !! s and p only, so the Cartesian form `make_efp_potential` forces and the
      !! spherical form `run_fmo2` builds are the same basis. Above d they are
      !! not, and test one would then compare two different models.

   real(dp), parameter :: TOL = 1.0e-9_dp
      !! Every identity below is exact, so what this has to clear is the SCF
      !! convergence the pieces are held to -- 1e-10 on the energy -- summed
      !! over the handful of SCFs a trimer needs.

   real(dp), parameter :: SPACING_12 = 4.0_dp
   real(dp), parameter :: SPACING_13 = 12.0_dp
      !! Where waters two and three sit along `x`, in Angstrom, from water one.
      !! Oxygen's Bondi radius is 1.52, so `R_IJ` is the separation over 3.04:
      !! 1.32 for the first pair, 3.95 and 2.63 for the other two. `R_cut = 2.0`
      !! therefore splits them one to two, which is what test four needs, and
      !! the two limits sweep past both ends of that range.

   ! Built once. Every test wants the same three potentials, and they are
   ! identical every time.
   type(efp_fragment_t), save :: cached_frag(3)
   real(dp), save :: cached_mono(3) = 0.0_dp
   real(dp), save :: cached_dimer_12 = 0.0_dp
   logical, save :: cached_ready = .false.

contains

   subroutine collect_mqc_czt_efmo_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("efmo_all_pairs_quantum_is_fmo2_plus_many_body_induction", &
                               test_all_quantum), &
                  new_unittest("efmo_no_pair_quantum_is_monomers_plus_efp", test_no_quantum), &
                  new_unittest("efmo_two_fragments_is_the_dimer_energy", test_two_fragments), &
                  new_unittest("efmo_trimer_split_one_quantum_two_effective", test_mixed) &
                  ]
   end subroutine collect_mqc_czt_efmo_tests

   subroutine water_geometry(z, symbols, coords)
      !! One water, in Bohr, in the `yz` plane
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

   subroutine water_chain(n, z, symbols, xyz, owner)
      !! `n` copies of that water along `x`, at the spacings above
      !!
      !! `n = 2` is exactly waters one and two of the trimer, at the same
      !! coordinates, so the two-fragment test and the trimer tests share a
      !! geometry and a dimer energy.
      integer, intent(in) :: n
      integer, intent(out) :: z(:), owner(:)
      character(len=2), intent(out) :: symbols(:)
      real(dp), intent(out) :: xyz(:, :)

      integer :: zw(3), k, f, at
      character(len=2) :: sw(3)
      real(dp) :: coords(3, 3), offset(3)

      call water_geometry(zw, sw, coords)
      at = 0
      do f = 1, n
         offset = 0.0_dp
         if (f == 2) offset(1) = SPACING_12*ANG
         if (f == 3) offset(1) = SPACING_13*ANG
         do k = 1, 3
            at = at + 1
            z(at) = zw(k)
            symbols(at) = sw(k)
            xyz(:, at) = coords(:, k) + offset
            owner(at) = f
         end do
      end do
   end subroutine water_chain

   subroutine efmo_settings(opts)
      !! The options every case runs with
      !!
      !! MAKEFP's own SCF defaults, which `run_efmo` uses for its dimers too:
      !! the reference energies computed here have to come from the same SCF, or
      !! an identity exact in exact arithmetic fails on convergence.
      type(efmo_options_t), intent(out) :: opts

      opts%basis = BASIS
      opts%scf_max_iter = 200
      opts%scf_energy_tol = 1.0e-10_dp
      opts%scf_density_tol = 1.0e-8_dp
      opts%scf_grad_tol = 1.0e-8_dp
      opts%scf%grad_tol = 1.0e-8_dp
   end subroutine efmo_settings

   subroutine build_reference(err)
      !! The three potentials, their `E_I^0`, and the 1-2 dimer's RHF energy
      type(error_t), intent(inout) :: err

      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)
      type(efp_potential_t) :: pot
      type(efmo_options_t) :: opts
      integer :: k
      integer, allocatable :: idx(:)

      if (cached_ready) return
      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)

      do k = 1, 3
         idx = atoms_of(owner, k)
         call make_efp_potential(z(idx), symbols(idx), xyz(:, idx), BASIS, "FRAG", pot, &
                                 err, charge=0, &
                                 energy_tol=opts%scf_energy_tol, &
                                 density_tol=opts%scf_density_tol, &
                                 grad_tol_in=opts%scf_grad_tol, scf_in=opts%scf, &
                                 max_iter_in=opts%scf_max_iter)
         if (err%has_error()) return
         cached_mono(k) = pot%scf_energy
         call potential_to_fragment(pot, cached_frag(k), err)
         call pot%destroy()
         if (err%has_error()) return
      end do

      call dimer_rhf(z(1:6), symbols(1:6), xyz(:, 1:6), 0, opts, cached_dimer_12, err)
      if (err%has_error()) return
      cached_ready = .true.
   end subroutine build_reference

   pure function atoms_of(owner, k) result(idx)
      !! The indices of fragment `k`'s atoms
      integer, intent(in) :: owner(:)
      integer, intent(in) :: k
      integer, allocatable :: idx(:)

      integer :: i, n

      n = count(owner == k)
      allocate (idx(n))
      n = 0
      do i = 1, size(owner)
         if (owner(i) == k) then
            n = n + 1
            idx(n) = i
         end if
      end do
   end function atoms_of

   subroutine dimer_rhf(z, symbols, xyz, charge, opts, energy, err)
      !! `E_IJ^0` computed here, the way `run_efmo` computes it
      !!
      !! Deliberately a second implementation rather than a call into the
      !! module: what test three and test four assert is that the orchestrator's
      !! dimer is *this* number, and reusing its own routine would assert
      !! nothing. Cartesian, because the monomer SCFs behind the potentials are.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: charge
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: guess_density(:, :)
      integer :: guess_kind

      energy = 0.0_dp
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err, force_cartesian=.true.)
      if (err%has_error()) return
      call build_restricted_guess(mol, "auto", guess_kind, guess_density, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, sum(z) - charge, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, .false., scf, err, &
                       guess=guess_kind, guess_density=guess_density, &
                       grad_tol=opts%scf_grad_tol, scf=opts%scf)
      call mol%destroy()
      if (err%has_error()) return
      energy = scf%energy
   end subroutine dimer_rhf

   function total_induction(frags, err) result(energy)
      !! `E_pol^total` over the fragments given
      type(efp_fragment_t), intent(in) :: frags(:)
      type(error_t), intent(inout) :: err
      real(dp) :: energy

      type(efp_system_t) :: system
      real(dp), allocatable :: shifts(:, :)

      allocate (shifts(3, size(frags)), source=0.0_dp)
      energy = 0.0_dp
      call build_efp_system(frags, shifts, system, err)
      if (err%has_error()) return
      energy = polarization_energy(system, frags, err)
      call system%destroy()
   end function total_induction

   function pair_induction_sum(frags, err) result(energy)
      !! `sum_{I<J} E_IJ^pol` over every pair of the fragments given
      type(efp_fragment_t), intent(in) :: frags(:)
      type(error_t), intent(inout) :: err
      real(dp) :: energy

      real(dp) :: zero(3)
      integer :: a, b

      zero = 0.0_dp
      energy = 0.0_dp
      do a = 1, size(frags) - 1
         do b = a + 1, size(frags)
            energy = energy + pair_polarization_energy(frags(a), frags(b), zero, zero, err)
            if (err%has_error()) return
         end do
      end do
   end function pair_induction_sum

   subroutine test_all_quantum(error)
      !! Limit one: `R_cut` huge, so EFMO is FMO2 in vacuo plus many-body induction
      !!
      !! With every pair quantum, eq 6 has no far sum and reads
      !!
      !!     sum_I E_I^0 + sum_IJ (E_IJ^0 - E_I^0 - E_J^0)  -  sum_IJ E_IJ^pol
      !!                                                    +  E_pol^total
      !!
      !! whose first line is exactly the many-body expansion at level two over
      !! in-vacuo fragments -- `run_fmo2` with the embedding off and the MBE
      !! assembly, which is a completely independent implementation of that sum.
      !! So the difference between the two totals must be the induction
      !! remainder and nothing else.
      !!
      !! **That remainder is not small.** It is the part of the induction no
      !! pair carries, and for this trimer it is a sizeable fraction of the
      !! total; the plan's earlier expectation that the two energies nearly
      !! agree was wrong. Only the identity is asserted.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(fmo_options_t) :: fmo_opts
      type(fmo_result_t) :: fmo_res
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9), e_pol_total, e_pol_pairs, expected

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 3, message="a huge cutoff left a pair effective")
      if (allocated(error)) return
      call check(error, res%n_efp_pairs, 0, message="a huge cutoff kept an EFP pair")
      if (allocated(error)) return

      ! FMO2 with no embedding and the MBE assembly: sum E_I + sum (E_IJ - E_I - E_J),
      ! every fragment solved in vacuo. The same first line as above, from other code.
      fmo_opts%basis = BASIS
      fmo_opts%esp = "none"
      fmo_opts%expansion = "mbe"
      fmo_opts%level = 2
      fmo_opts%scf_max_iter = 200
      fmo_opts%scf_energy_tol = 1.0e-10_dp
      fmo_opts%scf_density_tol = 1.0e-8_dp
      fmo_opts%scf%grad_tol = 1.0e-8_dp
      call run_fmo2(z, symbols, xyz, owner, fmo_opts, fmo_res, err)
      call check(error,.not. err%has_error(), "run_fmo2 failed: "//err%get_full_trace())
      if (allocated(error)) return

      e_pol_total = total_induction(cached_frag, err)
      e_pol_pairs = pair_induction_sum(cached_frag, err)
      call check(error,.not. err%has_error(), "the induction terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      expected = fmo_res%energy + (e_pol_total - e_pol_pairs)
      call check(error, res%energy, expected, thr=TOL, &
                 message="EFMO with every pair quantum is not FMO2 in vacuo plus the "// &
                 "many-body induction")
      if (allocated(error)) return

      ! And that the orchestrator's own two induction sums are the ones just
      ! computed, so a failure above says which half moved.
      call check(error, res%polarization_total, e_pol_total, thr=1.0e-12_dp, &
                 message="E_pol^total")
      if (allocated(error)) return
      call check(error, res%pair_polarization, e_pol_pairs, thr=1.0e-12_dp, &
                 message="sum E_IJ^pol")
      if (allocated(error)) return
      ! The remainder is a term of the method, not a residue: assert it is not
      ! negligible, so a version that silently dropped it would fail here.
      call check(error, abs(e_pol_total - e_pol_pairs) > 1.0e-6_dp, &
                 "the many-body induction vanished, so this geometry cannot tell a "// &
                 "dropped remainder from a correct one")
   end subroutine test_all_quantum

   subroutine test_no_quantum(error)
      !! Limit two: `R_cut` zero, so EFMO is in-vacuo monomers plus EFP-EFP
      !!
      !! With no quantum dimer there is no near sum at all, and eq 6 reads
      !! `sum_I E_I^0` plus the four pair terms over every pair plus
      !! `E_pol^total` -- which is exactly `sum_I E_I^0` plus what
      !! `efp_interaction_energy` returns for the same fragments, since that
      !! routine's five terms are those four and that induction.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(efp_energy_t) :: efp
      real(dp) :: shifts(3, 3), expected
      integer :: z(9), owner(9)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 0.0_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 0, message="a zero cutoff left a pair quantum")
      if (allocated(error)) return
      call check(error, res%n_efp_pairs, 3, message="a zero cutoff lost an EFP pair")
      if (allocated(error)) return
      call check(error, res%dimer_correction, 0.0_dp, thr=0.0_dp, &
                 message="there is no quantum dimer, so there is no dimer correction")
      if (allocated(error)) return
      call check(error, res%pair_polarization, 0.0_dp, thr=0.0_dp, &
                 message="there is no quantum dimer, so nothing subtracts pair induction")
      if (allocated(error)) return

      shifts = 0.0_dp
      efp = efp_interaction_energy(cached_frag, shifts, err)
      call check(error,.not. err%has_error(), "the EFP energy failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      expected = sum(cached_mono) + efp%total
      call check(error, res%energy, expected, thr=TOL, &
                 message="EFMO with no quantum dimer is not the in-vacuo monomers "// &
                 "plus the EFP-EFP interaction")
      if (allocated(error)) return
      call check(error, res%monomer_sum, sum(cached_mono), thr=TOL, &
                 message="the monomer sum is not the sum of the potentials' own SCFs")
   end subroutine test_no_quantum

   subroutine test_two_fragments(error)
      !! Limit three: two fragments, so EFMO is the dimer's own RHF energy
      !!
      !! `E_IJ^pol` on a two-fragment system *is* `E_pol^total` -- the same
      !! solver on the same system -- so the last term of eq 6 and the
      !! subtraction inside the near sum cancel exactly, and what is left is
      !! `E_I^0 + E_J^0 + (E_IJ^0 - E_I^0 - E_J^0) = E_IJ^0`.
      !!
      !! The sharpest of the four: it says the induction bookkeeping is right
      !! *and* that the monomer energies entering the difference are the ones
      !! the potentials were built from. Off by a wrong monomer energy, this
      !! fails by the size of an SCF energy rather than of an interaction.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      integer :: z(6), owner(6)
      character(len=2) :: symbols(6)
      real(dp) :: xyz(3, 6)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(2, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 1.0e6_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 1, message="the one pair is not quantum")
      if (allocated(error)) return
      call check(error, res%pair_polarization, res%polarization_total, thr=1.0e-12_dp, &
                 message="on two fragments the pair induction is not the total")
      if (allocated(error)) return
      call check(error, res%energy, cached_dimer_12, thr=TOL, &
                 message="EFMO on two fragments is not the dimer's in-vacuo RHF energy")
   end subroutine test_two_fragments

   subroutine test_mixed(error)
      !! Limit four: one quantum dimer and two effective ones, assembled by hand
      !!
      !! The trimer's separations are 1.32, 3.95 and 2.63 in contact units, so
      !! `R_cut = 2.0` puts the close pair in the quantum list and the other two
      !! in the effective one. This is the only case where both halves of eq 6
      !! are non-empty, and it is where a sign or a double count that cancels in
      !! the limits shows up.
      !!
      !! Every piece of the expected total comes from a Phase 1 routine called
      !! directly here: the monomer energies from the cached potentials, the
      !! dimer from an independent RHF, the pair induction and the four far
      !! terms from their own entries, and the total induction over all three.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(efmo_options_t) :: opts
      type(efmo_result_t) :: res
      type(efp_pair_energy_t), allocatable :: far(:)
      real(dp) :: shifts(3, 3), zero(3), e_pair_pol, e_pol_total, expected
      integer :: z(9), owner(9), far_pairs(2, 2)
      character(len=2) :: symbols(9)
      real(dp) :: xyz(3, 9)

      call build_reference(err)
      call check(error,.not. err%has_error(), "building the reference failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call water_chain(3, z, symbols, xyz, owner)
      call efmo_settings(opts)
      opts%rcut = 2.0_dp
      call run_efmo(z, symbols, xyz, owner, [0, 0, 0], opts, res, err)
      call check(error,.not. err%has_error(), "run_efmo failed: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, res%n_qm_pairs, 1, message="the split did not keep one QM dimer")
      if (allocated(error)) return
      call check(error, res%n_efp_pairs, 2, message="the split did not leave two EFP dimers")
      if (allocated(error)) return
      call check(error, res%pairs(1)%i == 1 .and. res%pairs(1)%j == 2, &
                 "the quantum dimer is not the close pair")
      if (allocated(error)) return
      ! The separations this geometry was built for, so a change to the radii or
      ! to the chain would fail here rather than silently reclassify a pair.
      call check(error, res%pairs(1)%r, SPACING_12/(2.0_dp*1.52_dp), thr=1.0e-10_dp, &
                 message="R_12 is not the oxygen separation over twice the oxygen radius")
      if (allocated(error)) return

      shifts = 0.0_dp
      zero = 0.0_dp
      far_pairs = reshape([1, 3, 2, 3], [2, 2])
      far = efp_pair_terms(cached_frag, shifts, far_pairs, err, charge_transfer_on=.true.)
      e_pair_pol = pair_polarization_energy(cached_frag(1), cached_frag(2), zero, zero, err)
      e_pol_total = total_induction(cached_frag, err)
      call check(error,.not. err%has_error(), "the reference terms failed: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      expected = sum(cached_mono) &
                 + (cached_dimer_12 - cached_mono(1) - cached_mono(2) - e_pair_pol) &
                 + sum(far%total) + e_pol_total
      call check(error, res%energy, expected, thr=TOL, &
                 message="the trimer total is not the sum assembled from the pieces")
      if (allocated(error)) return

      ! And each of the six sums separately, so a failure names its term.
      call check(error, res%monomer_sum, sum(cached_mono), thr=TOL, message="sum E_I^0")
      if (allocated(error)) return
      call check(error, res%dimer_correction, &
                 cached_dimer_12 - cached_mono(1) - cached_mono(2), thr=TOL, &
                 message="the QM dimer correction")
      if (allocated(error)) return
      call check(error, res%pair_polarization, e_pair_pol, thr=1.0e-12_dp, &
                 message="sum E_IJ^pol")
      if (allocated(error)) return
      call check(error, res%far_electrostatics, sum(far%electrostatics), thr=1.0e-12_dp, &
                 message="the far Coulomb sum")
      if (allocated(error)) return
      call check(error, res%far_dispersion, sum(far%dispersion), thr=1.0e-12_dp, &
                 message="the far dispersion sum")
      if (allocated(error)) return
      call check(error, res%far_exchange_repulsion, sum(far%exchange_repulsion), &
                 thr=1.0e-12_dp, message="the far exchange repulsion sum")
      if (allocated(error)) return
      call check(error, res%far_charge_transfer, sum(far%charge_transfer), thr=1.0e-12_dp, &
                 message="the far charge transfer sum")
      if (allocated(error)) return
      call check(error, res%polarization_total, e_pol_total, thr=1.0e-12_dp, &
                 message="E_pol^total")
   end subroutine test_mixed

end module test_mqc_czt_efmo

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_efmo, only: collect_mqc_czt_efmo_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_efmo", collect_mqc_czt_efmo_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
