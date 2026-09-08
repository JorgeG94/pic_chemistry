!! The effective fragment molecular orbital energy
module mqc_czt_efmo
   !! EFMO (Sattasathuchana, Xu, Bertoni, Kim, Leang, Pham, Gordon, JCTC 20,
   !! 2445 (2024), eq 6; Steinmann, Fedorov, Jensen, JPCA 114, 8705 (2010)):
   !!
   !!     E = sum_I E_I^0
   !!       + sum_{I<J, R_IJ <= R_cut} ( E_IJ^0 - E_I^0 - E_J^0 - E_IJ^pol )
   !!       + sum_{I<J, R_IJ >  R_cut} ( E_IJ^Coul + E_IJ^disp + E_IJ^ExRep + E_IJ^CT )
   !!       + E_pol^total
   !!
   !! **Nothing here is self-consistent across fragments.** `E_I^0` and `E_IJ^0`
   !! are *in vacuo* energies -- no embedding field, no monomer loop -- which is
   !! what separates EFMO from FMO and what lets a diffuse basis work: there is
   !! no neighbouring point charge for a diffuse function to fall onto. The
   !! coupling between fragments is carried entirely by the effective fragment
   !! potentials, one per fragment, each built by MAKEFP from the same SCF that
   !! produced `E_I^0`.
   !!
   !! **`E_IJ^pol` is subtracted from every near dimer, and it is not small.**
   !! A quantum dimer already contains the mutual induction of its two
   !! fragments, and `E_pol^total` -- the induction solved over every fragment
   !! at once -- contains it too, so one copy has to go. What survives,
   !! `E_pol^total - sum_IJ E_IJ^pol`, is the many-body part of the induction:
   !! measured on three waters at four Angstrom it is 44 per cent of the total,
   !! because the energy is quadratic in the field and the square of a sum keeps
   !! cross terms no pair has. It is a term of the method, not a residue.
   !!
   !! **Where each number comes from**, all of it Phase 1 work:
   !!
   !! * `E_I^0` is `pot%scf_energy`, the SCF `make_efp_potential` runs anyway.
   !! * the potential becomes an `efp_fragment_t` through `potential_to_fragment`
   !!   rather than through a written `.efp` file.
   !! * the QM/EFP split is `efmo_split_pairs` on the vdW-scaled `R_IJ` of eq 2.
   !! * the far pairs go through `efp_pair_terms`, which builds a two-fragment
   !!   system per pair with the charge-penetration screening on. **The
   !!   system-wide `electrostatic_energy` is deliberately never called here**:
   !!   it would put every fragment's points in one set, which is the same sum
   !!   only when *every* pair is a far pair, and there is no mask to exclude
   !!   the near ones.
   !! * `E_pol^total` is `polarization_energy` over every fragment at once,
   !!   which is the call `efp_interaction_energy` makes internally.
   !!
   !! **Cartesian throughout.** `make_efp_potential` forces 6d/10f, because a
   !! `.efp` is read by GAMESS, so `E_I^0` is a Cartesian energy. The dimer SCFs
   !! here are built the same way: with a spherical dimer the difference
   !! `E_IJ^0 - E_I^0 - E_J^0` would be taken between two different models and
   !! would not be an interaction energy at all. Below d functions the two forms
   !! coincide and the choice is invisible.
   !!
   !! **One rank, closed shell, whole molecules.** Distributing the monomers and
   !! dimers is Phase 4 and covalent fragments are Phase 5; a partition that
   !! cuts a bond is refused here rather than capped, since a cap's multipoles
   !! would act on the partner across the cut.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_calculation_defaults, only: DEFAULT_VDW_SCALE, DEFAULT_DYNAMIC_TOL, &
                                       DEFAULT_DYNAMIC_MAXITER, DEFAULT_RESPONSE_BATCH, &
                                       EFP_RESPONSE_AUTO
   use mqc_scf_types, only: scf_numerics_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_atomic_guess, only: build_restricted_guess
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential
   use mqc_czt_efp_read, only: efp_fragment_t
   use mqc_czt_efp_convert, only: potential_to_fragment
   use mqc_czt_efp_energy, only: efp_pair_energy_t, efp_pair_terms, &
                                 pair_polarization_energy
   use mqc_czt_efp_interaction, only: efp_system_t, build_efp_system, polarization_energy
   use mqc_czt_efmo_pairs, only: efmo_split_pairs, efmo_pair_distance
   implicit none
   private

   public :: efmo_options_t
   public :: efmo_pair_t
   public :: efmo_result_t
   public :: run_efmo

   type :: efmo_options_t
      !! What to run, and how hard
      character(len=64) :: basis = "6-31g"
      real(dp) :: rcut = 2.0_dp
         !! `R_cut` of eq 2, **unitless**: each interatomic distance is divided
         !! by the two van der Waals radii, so 1 is contact. A pair at or inside
         !! it is a quantum dimer, a pair beyond it is four EFP terms. At or
         !! below zero every pair is effective, which is EFP with in-vacuo
         !! monomers; huge, every pair is quantum, which is FMO2 in vacuo plus
         !! the many-body induction. Both limits run, and both are tested.
      logical :: charge_transfer = .true.
         !! Include `E_IJ^CT` in the far pairs. GAMESS's EFMO has it; the 2012
         !! method left it out, so it is switchable rather than assumed.
      character(len=32) :: guess = "auto"
         !! Initial guess for every SCF here, monomer and dimer alike.
      character(len=64) :: aux_basis = ""
         !! Fit the MAKEFP response Hessian against this basis. Empty is exact.
      type(scf_numerics_t) :: scf
         !! How every SCF is driven -- the accelerator, DIIS subspace, level
         !! shift, linear-dependence threshold and incremental Fock switch. Its
         !! `max_iter`, `energy_tol` and `density_tol` are not read: the four
         !! fields below are, and they are passed positionally so they win.
      integer :: scf_max_iter = 200
      real(dp) :: scf_energy_tol = 1.0e-10_dp
      real(dp) :: scf_density_tol = 1.0e-8_dp
      real(dp) :: scf_grad_tol = 1.0e-8_dp
         !! MAKEFP's own defaults, deliberately tighter than a whole-system
         !! run's. The dimer SCF has to be converged to the same place the
         !! monomer one is, because their difference is the interaction energy
         !! and is four orders smaller than either.
      real(dp) :: vdw_scale = DEFAULT_VDW_SCALE
      logical :: quadrupole_blocks = .true.
      real(dp) :: dynamic_tolerance = DEFAULT_DYNAMIC_TOL
      integer :: dynamic_maxiter = DEFAULT_DYNAMIC_MAXITER
      logical :: allow_crap_response = .false.
      integer :: response = EFP_RESPONSE_AUTO
      integer :: response_batch = DEFAULT_RESPONSE_BATCH
         !! `keywords.efp`, forwarded whole to `make_efp_potential`. Nothing
         !! here is read by this module; it configures the stages of MAKEFP
         !! after the SCF.
      logical :: verbose = .false.
         !! Let MAKEFP report its own stages. The EFMO table is written at info
         !! level either way.
   end type efmo_options_t

   type :: efmo_pair_t
      !! One fragment pair, in whichever of the two lists it landed
      integer :: i = 0, j = 0
      real(dp) :: r = 0.0_dp
         !! `R_IJ`, the vdW-scaled closest approach. Unitless.
      logical :: qm = .false.
         !! True: a dimer SCF ran and `e_dimer`/`e_pair_pol` are filled. False:
         !! the four EFP terms below are.
      real(dp) :: e_dimer = 0.0_dp             !! `E_IJ^0`, in vacuo
      real(dp) :: e_pair_pol = 0.0_dp          !! `E_IJ^pol`, subtracted
      real(dp) :: electrostatics = 0.0_dp
      real(dp) :: dispersion = 0.0_dp
      real(dp) :: exchange_repulsion = 0.0_dp
      real(dp) :: charge_transfer = 0.0_dp
   end type efmo_pair_t

   type :: efmo_result_t
      !! The total, and the six sums it is made of
      !!
      !! `energy` is exactly
      !! `monomer_sum + dimer_correction - pair_polarization + far_electrostatics
      !! + far_dispersion + far_exchange_repulsion + far_charge_transfer
      !! + polarization_total`, with `pair_polarization` held positive and
      !! subtracted, since that is how eq 6 writes it.
      real(dp) :: energy = 0.0_dp
      real(dp) :: monomer_sum = 0.0_dp
         !! `sum_I E_I^0`
      real(dp) :: dimer_correction = 0.0_dp
         !! `sum (E_IJ^0 - E_I^0 - E_J^0)` over the quantum dimers
      real(dp) :: pair_polarization = 0.0_dp
         !! `sum E_IJ^pol` over the quantum dimers, **subtracted** from the total
      real(dp) :: far_electrostatics = 0.0_dp
      real(dp) :: far_dispersion = 0.0_dp
      real(dp) :: far_exchange_repulsion = 0.0_dp
      real(dp) :: far_charge_transfer = 0.0_dp
         !! The four EFP terms, summed over the effective dimers
      real(dp) :: polarization_total = 0.0_dp
         !! `E_pol^total`, induction over every fragment at once
      real(dp), allocatable :: monomer_energy(:)   !! `E_I^0`
      type(efmo_pair_t), allocatable :: pairs(:)
         !! Every pair, quantum ones first, each carrying its `R_IJ`
      integer :: n_qm_pairs = 0
      integer :: n_efp_pairs = 0
   end type efmo_result_t

contains

   subroutine run_efmo(atomic_numbers, symbols, coordinates, owner, fragment_charges, &
                       opts, res, error)
      !! One EFMO energy, from a system already partitioned into fragments
      !!
      !! `owner(i)` is the fragment of atom `i`, numbered from one with no gaps
      !! -- the same partition `run_fmo2` takes. Coordinates are Bohr.
      !! `fragment_charges(k)` is fragment `k`'s net charge; a charged fragment
      !! needs nothing extra, its monomer SCF and every dimer holding it just
      !! carry the charge.
      integer, intent(in) :: atomic_numbers(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, n_atoms), Bohr
      integer, intent(in) :: owner(:)
      integer, intent(in) :: fragment_charges(:)    !! (n_fragments)
      type(efmo_options_t), intent(in) :: opts
      type(efmo_result_t), intent(out) :: res
      type(error_t), intent(inout) :: error

      type(efp_fragment_t), allocatable :: frags(:)
      type(efp_pair_energy_t), allocatable :: far(:)
      real(dp), allocatable :: shifts(:, :)
      integer, allocatable :: qm_pairs(:, :), efp_pairs(:, :)
      integer, allocatable :: idx_i(:), idx_j(:)
      integer, allocatable :: count_of(:)
      integer :: n_atoms, n_frag, k, p

      n_atoms = size(atomic_numbers)
      if (size(owner) /= n_atoms .or. size(coordinates, 2) /= n_atoms &
          .or. size(symbols) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "efmo: the owner list, the symbols and the "// &
                        "coordinates must cover every atom")
         return
      end if
      if (n_atoms < 1) then
         call error%set(ERROR_VALIDATION, "efmo: there are no atoms to fragment")
         return
      end if
      if (minval(owner) < 1) then
         call error%set(ERROR_VALIDATION, "efmo: every atom must belong to a fragment "// &
                        "numbered from one")
         return
      end if
      n_frag = maxval(owner)
      if (size(fragment_charges) /= n_frag) then
         call error%set(ERROR_VALIDATION, "efmo: the system has "//to_char(n_frag)// &
                        " fragments but "//to_char(size(fragment_charges))//" charges")
         return
      end if

      call fragment_counts(owner, n_frag, count_of, error)
      if (error%has_error()) return

      allocate (res%monomer_energy(n_frag), source=0.0_dp)
      allocate (frags(n_frag), shifts(3, n_frag))
      ! Every potential is built at the geometry it is used at, so no fragment
      ! is placed and no rigid transform is looked for. That is why
      ! `place_fragment` never appears here.
      shifts = 0.0_dp

      call build_potentials(atomic_numbers, symbols, coordinates, owner, count_of, &
                            fragment_charges, opts, frags, res%monomer_energy, error)
      if (error%has_error()) return
      res%monomer_sum = sum(res%monomer_energy)

      call efmo_split_pairs(owner, atomic_numbers, coordinates, opts%rcut, &
                            qm_pairs, efp_pairs, error)
      if (error%has_error()) return
      res%n_qm_pairs = size(qm_pairs, 2)
      res%n_efp_pairs = size(efp_pairs, 2)
      allocate (res%pairs(res%n_qm_pairs + res%n_efp_pairs))

      call quantum_dimers(atomic_numbers, symbols, coordinates, owner, &
                          fragment_charges, qm_pairs, frags, shifts, opts, res, error)
      if (error%has_error()) return

      ! The far half. `efp_pair_terms` takes the pair list directly, so the near
      ! pairs contribute nothing here -- which is the point, since their
      ! electrostatics, exchange and dispersion are inside their dimer SCF.
      far = efp_pair_terms(frags, shifts, efp_pairs, error, &
                           charge_transfer_on=opts%charge_transfer)
      if (error%has_error()) return
      do k = 1, res%n_efp_pairs
         p = res%n_qm_pairs + k
         res%pairs(p)%i = efp_pairs(1, k)
         res%pairs(p)%j = efp_pairs(2, k)
         res%pairs(p)%qm = .false.
         idx_i = gather(owner, efp_pairs(1, k))
         idx_j = gather(owner, efp_pairs(2, k))
         call efmo_pair_distance(atomic_numbers(idx_i), coordinates(:, idx_i), &
                                 atomic_numbers(idx_j), coordinates(:, idx_j), &
                                 res%pairs(p)%r, error)
         if (error%has_error()) return
         res%pairs(p)%electrostatics = far(k)%electrostatics
         res%pairs(p)%dispersion = far(k)%dispersion
         res%pairs(p)%exchange_repulsion = far(k)%exchange_repulsion
         res%pairs(p)%charge_transfer = far(k)%charge_transfer
      end do
      res%far_electrostatics = sum(far%electrostatics)
      res%far_dispersion = sum(far%dispersion)
      res%far_exchange_repulsion = sum(far%exchange_repulsion)
      res%far_charge_transfer = sum(far%charge_transfer)

      call total_polarization(frags, shifts, res%polarization_total, error)
      if (error%has_error()) return

      res%energy = res%monomer_sum + res%dimer_correction - res%pair_polarization &
                   + res%far_electrostatics + res%far_dispersion &
                   + res%far_exchange_repulsion + res%far_charge_transfer &
                   + res%polarization_total

      call report(res, opts)
   end subroutine run_efmo

   subroutine fragment_counts(owner, n_frag, count_of, error)
      !! How many atoms each fragment holds
      !!
      !! The atoms of a fragment need not be contiguous in the deck, so there is
      !! no offset to keep and `gather` builds the index list instead. What is
      !! checked here is that the numbering has no gap, which would otherwise
      !! reach `make_efp_potential` as a molecule with no atoms.
      integer, intent(in) :: owner(:)
      integer, intent(in) :: n_frag
      integer, allocatable, intent(out) :: count_of(:)
      type(error_t), intent(inout) :: error

      integer :: i

      allocate (count_of(n_frag))
      count_of = 0
      do i = 1, size(owner)
         count_of(owner(i)) = count_of(owner(i)) + 1
      end do
      if (any(count_of == 0)) then
         call error%set(ERROR_VALIDATION, "efmo: fragment numbering has a gap -- some "// &
                        "fragment between one and "//to_char(n_frag)//" holds no atoms")
      end if
   end subroutine fragment_counts

   pure function gather(owner, k) result(idx)
      !! The system indices of fragment `k`'s atoms, in deck order
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
   end function gather

   subroutine build_potentials(z, symbols, xyz, owner, count_of, charges, opts, &
                               frags, monomer_energy, error)
      !! One MAKEFP per fragment, and `E_I^0` off the same SCF
      !!
      !! **The whole cost of an EFMO run is here.** A potential is an SCF, a
      !! localization and twelve frequency-dependent response solves, against
      !! one SCF for a dimer, so the monomer loop dominates and is what Phase 4
      !! distributes.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: owner(:), count_of(:), charges(:)
      type(efmo_options_t), intent(in) :: opts
      type(efp_fragment_t), intent(out) :: frags(:)
      real(dp), intent(out) :: monomer_energy(:)
      type(error_t), intent(inout) :: error

      type(efp_potential_t) :: pot
      integer, allocatable :: idx(:)
      character(len=:), allocatable :: aux
      integer :: k

      do k = 1, size(count_of)
         idx = gather(owner, k)
         call logger%verbose("  efmo: fragment "//to_char(k)//" of "// &
                             to_char(size(count_of))//", "//to_char(size(idx))//" atoms")
         ! One call whether or not an auxiliary basis was named: an absent
         ! optional passed on as an actual argument arrives absent.
         if (len_trim(opts%aux_basis) > 0) then
            aux = trim(opts%aux_basis)
            call make_efp_potential(z(idx), symbols(idx), xyz(:, idx), trim(opts%basis), &
                                    "FRAG"//to_char(k), pot, error, charge=charges(k), &
                                    verbose=opts%verbose, aux_basis=aux, &
                                    guess=trim(opts%guess), &
                                    energy_tol=opts%scf_energy_tol, &
                                    density_tol=opts%scf_density_tol, &
                                    grad_tol_in=opts%scf_grad_tol, scf_in=opts%scf, &
                                    max_iter_in=opts%scf_max_iter, &
                                    vdwscl=opts%vdw_scale, &
                                    quadrupole_blocks=opts%quadrupole_blocks, &
                                    dynamic_tol=opts%dynamic_tolerance, &
                                    dynamic_maxiter=opts%dynamic_maxiter, &
                                    response=opts%response, &
                                    allow_crap_response=opts%allow_crap_response, &
                                    response_batch=opts%response_batch)
         else
            call make_efp_potential(z(idx), symbols(idx), xyz(:, idx), trim(opts%basis), &
                                    "FRAG"//to_char(k), pot, error, charge=charges(k), &
                                    verbose=opts%verbose, guess=trim(opts%guess), &
                                    energy_tol=opts%scf_energy_tol, &
                                    density_tol=opts%scf_density_tol, &
                                    grad_tol_in=opts%scf_grad_tol, scf_in=opts%scf, &
                                    max_iter_in=opts%scf_max_iter, &
                                    vdwscl=opts%vdw_scale, &
                                    quadrupole_blocks=opts%quadrupole_blocks, &
                                    dynamic_tol=opts%dynamic_tolerance, &
                                    dynamic_maxiter=opts%dynamic_maxiter, &
                                    response=opts%response, &
                                    allow_crap_response=opts%allow_crap_response, &
                                    response_batch=opts%response_batch)
         end if
         if (error%has_error()) return

         ! `E_I^0` *is* that SCF. Running a second one here would be a second
         ! determinant nothing compares against.
         monomer_energy(k) = pot%scf_energy
         call potential_to_fragment(pot, frags(k), error)
         call pot%destroy()
         if (error%has_error()) return
      end do
   end subroutine build_potentials

   subroutine quantum_dimers(z, symbols, xyz, owner, charges, qm_pairs, &
                             frags, shifts, opts, res, error)
      !! Every near pair: one in-vacuo dimer SCF, and its pair induction
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: owner(:), charges(:)
      integer, intent(in) :: qm_pairs(:, :)
      type(efp_fragment_t), intent(in) :: frags(:)
      real(dp), intent(in) :: shifts(:, :)
      type(efmo_options_t), intent(in) :: opts
      type(efmo_result_t), intent(inout) :: res
      type(error_t), intent(inout) :: error

      integer, allocatable :: idx_i(:), idx_j(:), idx(:)
      real(dp) :: e_dimer, e_pol
      integer :: k, a, b

      do k = 1, size(qm_pairs, 2)
         a = qm_pairs(1, k)
         b = qm_pairs(2, k)
         idx_i = gather(owner, a)
         idx_j = gather(owner, b)
         idx = [idx_i, idx_j]

         res%pairs(k)%i = a
         res%pairs(k)%j = b
         res%pairs(k)%qm = .true.
         call efmo_pair_distance(z(idx_i), xyz(:, idx_i), z(idx_j), xyz(:, idx_j), &
                                 res%pairs(k)%r, error)
         if (error%has_error()) return

         call logger%verbose("  efmo: dimer "//to_char(a)//"-"//to_char(b)// &
                             ", "//to_char(size(idx))//" atoms")
         call dimer_energy(z(idx), symbols(idx), xyz(:, idx), charges(a) + charges(b), &
                           opts, e_dimer, error)
         if (error%has_error()) return

         ! The same induction solver on the two fragments alone -- the same
         ! screening, the same static field rank, the same tolerance as the
         ! total below. A pair solved any other way would leave a residue in
         ! `E_pol^total - sum E_IJ^pol` that looks like three-body induction.
         e_pol = pair_polarization_energy(frags(a), frags(b), shifts(:, a), &
                                          shifts(:, b), error)
         if (error%has_error()) return

         res%pairs(k)%e_dimer = e_dimer
         res%pairs(k)%e_pair_pol = e_pol
         res%dimer_correction = res%dimer_correction + e_dimer &
                                - res%monomer_energy(a) - res%monomer_energy(b)
         res%pair_polarization = res%pair_polarization + e_pol
      end do
   end subroutine quantum_dimers

   subroutine dimer_energy(z, symbols, xyz, charge, opts, energy, error)
      !! One dimer's restricted Hartree-Fock energy, in vacuo
      !!
      !! Cartesian, to match the monomer SCFs `make_efp_potential` ran; see the
      !! module header. No embedding of any kind: the neighbouring fragments'
      !! potentials are not felt by this SCF, and their interaction with the
      !! pair is carried by the EFP terms and the total induction instead.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: charge
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: guess_density(:, :)
      integer :: guess_kind, nelec

      energy = 0.0_dp
      nelec = sum(z) - charge
      if (nelec < 2 .or. mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "efmo: a dimer with "//to_char(nelec)// &
                        " electrons is not closed-shell. EFMO is restricted "// &
                        "Hartree-Fock for now, so the fragment charges have to "// &
                        "leave every fragment and every dimer with an even count.")
         return
      end if

      call build_czt_molecule(z, symbols, xyz, trim(opts%basis), mol, error, &
                              force_cartesian=.true.)
      if (error%has_error()) return
      call build_restricted_guess(mol, trim(opts%guess), guess_kind, guess_density, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      call run_czt_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, opts%verbose, scf, error, &
                       guess=guess_kind, guess_density=guess_density, &
                       grad_tol=opts%scf_grad_tol, scf=opts%scf)
      call mol%destroy()
      if (error%has_error()) return
      if (.not. scf%converged .and. .not. opts%scf%allow_crap_scf) then
         call error%set(ERROR_VALIDATION, "efmo: a dimer SCF did not converge, so the "// &
                        "pair correction it feeds is not trustworthy. Set "// &
                        "keywords.scf.allow_crap_scf to finish anyway.")
         return
      end if
      energy = scf%energy
   end subroutine dimer_energy

   subroutine total_polarization(frags, shifts, energy, error)
      !! `E_pol^total`: induction over every fragment at once
      !!
      !! The same call `efp_interaction_energy` makes internally, on the same
      !! system, so the pair terms subtracted from the near dimers cancel
      !! against exactly what is here.
      type(efp_fragment_t), intent(in) :: frags(:)
      real(dp), intent(in) :: shifts(:, :)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(efp_system_t) :: system

      energy = 0.0_dp
      if (size(frags) < 2) return
      call build_efp_system(frags, shifts, system, error)
      if (error%has_error()) return
      energy = polarization_energy(system, frags, error)
      call system%destroy()
   end subroutine total_polarization

   subroutine report(res, opts)
      !! The breakdown, at info level
      type(efmo_result_t), intent(in) :: res
      type(efmo_options_t), intent(in) :: opts

      character(len=160) :: line
      integer :: k

      call logger%info("============================================================")
      call logger%info("  EFMO, R_cut = "//to_char(opts%rcut)//" (unitless)")
      call logger%info("------------------------------------------------------------")
      call logger%info("  fragment           E_I^0 / Hartree")
      do k = 1, size(res%monomer_energy)
         write (line, "(A,I6,F26.10)") "  ", k, res%monomer_energy(k)
         call logger%info(trim(line))
      end do
      call logger%info("------------------------------------------------------------")
      call logger%info("  pair      R_IJ  class      contribution / Hartree")
      do k = 1, size(res%pairs)
         if (res%pairs(k)%qm) then
            write (line, "(A,I4,A,I4,F8.3,A,F22.10)") "  ", res%pairs(k)%i, "-", &
               res%pairs(k)%j, res%pairs(k)%r, "  QM  ", &
               res%pairs(k)%e_dimer - res%monomer_energy(res%pairs(k)%i) &
               - res%monomer_energy(res%pairs(k)%j) - res%pairs(k)%e_pair_pol
         else
            write (line, "(A,I4,A,I4,F8.3,A,F22.10)") "  ", res%pairs(k)%i, "-", &
               res%pairs(k)%j, res%pairs(k)%r, "  EFP ", &
               res%pairs(k)%electrostatics + res%pairs(k)%dispersion &
               + res%pairs(k)%exchange_repulsion + res%pairs(k)%charge_transfer
         end if
         call logger%info(trim(line))
      end do
      call logger%info("------------------------------------------------------------")
      call logger%info("  monomers            sum E_I^0        "//to_char(res%monomer_sum))
      call logger%info("  QM dimers           E_IJ - E_I - E_J "//to_char(res%dimer_correction))
      call logger%info("  QM dimers           - sum E_IJ^pol   "//to_char(-res%pair_polarization))
      call logger%info("  EFP dimers          Coulomb          "//to_char(res%far_electrostatics))
      call logger%info("  EFP dimers          dispersion       "//to_char(res%far_dispersion))
      call logger%info("  EFP dimers          exchange rep.    "//to_char(res%far_exchange_repulsion))
      call logger%info("  EFP dimers          charge transfer  "//to_char(res%far_charge_transfer))
      call logger%info("  all fragments       E_pol^total      "//to_char(res%polarization_total))
      call logger%info("------------------------------------------------------------")
      call logger%info("  QM dimers "//to_char(res%n_qm_pairs)//", EFP dimers "// &
                       to_char(res%n_efp_pairs))
      call logger%info("  EFMO total energy   "//to_char(res%energy)//" Hartree")
      call logger%info("============================================================")
   end subroutine report

end module mqc_czt_efmo
