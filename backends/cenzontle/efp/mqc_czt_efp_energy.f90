!! The EFP2 interaction energy, all five terms
module mqc_czt_efp_energy
   !! What the rest of the EFP code is for: given placed fragments, one number.
   !!
   !! The five terms live in two modules because they need different things.
   !! `mqc_czt_efp_interaction` works on a flattened point set and covers the terms that
   !! are geometry and stored tensors -- electrostatics, polarization, and the
   !! undamped `E6`. `mqc_czt_efp_pair` works on fragment *pairs* and covers everything
   !! needing integrals over two fragments' basis sets at once: exchange repulsion,
   !! charge transfer, and the damped dispersion, whose damping is an overlap between
   !! localized orbitals on different fragments. This module is where the two meet,
   !! and it exists so no caller has to know which term came from where.
   !!
   !! **Every term here is separately validated against GAMESS**, all for the
   !! same water dimer, which is why the result keeps the breakdown rather than
   !! returning a bare total: a regression then names the term it broke.
   !!
   !! **Dispersion is the damped sum `E6 + E7 + E8`**, which is what GAMESS
   !! totals when the potential carries the tensors for all three. The undamped
   !! `dispersion_energy_e6` in `mqc_czt_efp_interaction` is not used here.
   !!
   !! **Fragments arrive already turned.** `place_fragment` gives the rotation a deck
   !! implies and `mqc_czt_efp_rotate` applies it, so by the time a fragment reaches this
   !! module its own stored frame *is* the working frame and all that is left is the
   !! translation each term takes as an offset.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_efp_read, only: efp_fragment_t
   use mqc_czt_efp_interaction, only: efp_system_t, build_efp_system, &
                                      electrostatic_energy, polarization_energy
   use mqc_czt_efp_pair, only: exchange_repulsion, charge_transfer, &
                               dispersion_e6_damped, dispersion_e7_damped, &
                               dispersion_e8_damped
   use mqc_czt_efp_rotate, only: superpose
   implicit none
   private

   public :: efp_energy_t
   public :: efp_interaction_energy
   public :: place_fragment
   ! What EFMO needs of this module: one pair at a time, since its energy
   ! expression treats a pair as either quantum-mechanical or effective and
   ! never as part of one system-wide sum.
   public :: efp_pair_energy_t
   public :: efp_pair_terms
   public :: pair_polarization_energy

   ! How far a deck atom may sit from where the potential's own geometry puts it,
   ! after the rigid shift, before the placement is refused. A fragment is rigid,
   ! so this is not a fitting tolerance -- it is the width of the round trip
   ! through the file, which carries ten decimals.
   real(dp), parameter :: PLACEMENT_TOL = 1.0e-6_dp

   ! Every multipole rank the electrostatics carries: charges, dipoles,
   ! quadrupoles and octupoles.
   integer, parameter :: MAX_RANK = 3

   type :: efp_pair_energy_t
      !! The EFP terms of **one** fragment pair, polarization excluded
      !!
      !! What eq 6 of the EFMO paper sums over the far pairs: Coulomb,
      !! dispersion, exchange repulsion and charge transfer. Induction is absent
      !! by construction rather than by omission -- it is many-body and lives in
      !! the one `E_pol^total` over every fragment.
      integer :: i = 0, j = 0                  !! Which two fragments, as given
      real(dp) :: electrostatics = 0.0_dp
      real(dp) :: exchange_repulsion = 0.0_dp
      real(dp) :: dispersion = 0.0_dp     !! Damped `E6 + E7 + E8`
      real(dp) :: dispersion_e6 = 0.0_dp
      real(dp) :: dispersion_e7 = 0.0_dp
      real(dp) :: dispersion_e8 = 0.0_dp
      real(dp) :: charge_transfer = 0.0_dp
      real(dp) :: total = 0.0_dp          !! The four above, no polarization
   end type efp_pair_energy_t

   type :: efp_energy_t
      !! One interaction energy, kept broken down by term
      real(dp) :: electrostatics = 0.0_dp
      real(dp) :: polarization = 0.0_dp
      real(dp) :: exchange_repulsion = 0.0_dp
      real(dp) :: dispersion = 0.0_dp     !! Damped `E6 + E7 + E8`
      real(dp) :: dispersion_e6 = 0.0_dp
      real(dp) :: dispersion_e7 = 0.0_dp
      real(dp) :: dispersion_e8 = 0.0_dp
      real(dp) :: charge_transfer = 0.0_dp
      real(dp) :: total = 0.0_dp
   end type efp_energy_t

contains

   function efp_interaction_energy(fragments, translations, error) result(energy)
      !! The interaction energy of a set of placed fragments
      !!
      !! `translations(:, k)` places `fragments(k)`, in Bohr, as a rigid shift of the
      !! geometry the potential itself carries.
      !!
      !! Electrostatics and polarization are computed over the whole system at once --
      !! polarization *must* be, since the induced dipoles are solved together and are
      !! not a sum of pair terms. The other three are pairwise and are summed over
      !! `a < b`. Nothing here assumes two fragments.
      !!
      !! A term whose data the potential does not carry is left at zero rather
      !! than erroring -- a potential written without the dynamic polarizabilities
      !! has no dispersion to contribute. A term failing for any other reason
      !! propagates.
      type(efp_fragment_t), intent(in) :: fragments(:)
      real(dp), intent(in) :: translations(:, :)   !! (3, n_fragments), Bohr
      type(error_t), intent(inout) :: error
      type(efp_energy_t) :: energy

      type(efp_system_t) :: system
      integer :: n, a, b
      logical :: have_dynamic, have_lmo, have_ct

      n = size(fragments)
      if (size(translations, 1) /= 3 .or. size(translations, 2) /= n) then
         call error%set(ERROR_VALIDATION, "efp: one translation per fragment is "// &
                        "needed, as (3, n_fragments)")
         return
      end if
      if (n < 2) return      ! a single fragment interacts with nothing

      call build_efp_system(fragments, translations, system, error)
      if (error%has_error()) return

      energy%electrostatics = electrostatic_energy(system, MAX_RANK, screen=.true.)
      energy%polarization = polarization_energy(system, fragments, error)
      if (error%has_error()) return

      do a = 1, n
         do b = a + 1, n
            have_dynamic = fragments(a)%has_dynamic .and. fragments(b)%has_dynamic
            have_lmo = fragments(a)%has_lmo .and. fragments(b)%has_lmo
            have_ct = fragments(a)%has_ctvec .and. fragments(b)%has_ctvec &
                      .and. fragments(a)%has_ctfok .and. fragments(b)%has_ctfok

            if (have_lmo) then
               energy%exchange_repulsion = energy%exchange_repulsion &
                                           + exchange_repulsion(fragments(a), fragments(b), &
                                                                translations(:, a), translations(:, b), error)
               if (error%has_error()) return
            end if

            if (have_dynamic .and. have_lmo) then
               energy%dispersion_e6 = energy%dispersion_e6 &
                                      + dispersion_e6_damped(fragments(a), fragments(b), &
                                                             translations(:, a), translations(:, b), error)
               if (error%has_error()) return
               ! E7 and E8 need the higher tensor blocks, which a potential may omit
               ! -- and GAMESS treats their absence as switching the term off rather
               ! than as an error (`efinp.src:5559, 5577`), so this does too.
               if (fragments(a)%has_dipquad .and. fragments(b)%has_dipquad) then
                  energy%dispersion_e7 = energy%dispersion_e7 &
                                         + dispersion_e7_damped(fragments(a), fragments(b), &
                                                                translations(:, a), translations(:, b), error)
                  if (error%has_error()) return
               end if
               if (fragments(a)%has_quadquad .and. fragments(b)%has_quadquad) then
                  energy%dispersion_e8 = energy%dispersion_e8 &
                                         + dispersion_e8_damped(fragments(a), fragments(b), &
                                                                translations(:, a), translations(:, b), error)
                  if (error%has_error()) return
               end if
            end if

            if (have_ct) then
               energy%charge_transfer = energy%charge_transfer &
                                        + charge_transfer(fragments(a), fragments(b), &
                                                          translations(:, a), translations(:, b), error)
               if (error%has_error()) return
            end if
         end do
      end do

      energy%dispersion = energy%dispersion_e6 + energy%dispersion_e7 &
                          + energy%dispersion_e8
      energy%total = energy%electrostatics + energy%polarization &
                     + energy%exchange_repulsion + energy%dispersion &
                     + energy%charge_transfer

      call system%destroy()
   end function efp_interaction_energy

   function pair_polarization_energy(frag_i, frag_j, translation_i, translation_j, &
                                     error) result(energy)
      !! The induction energy of one *isolated* pair of fragments
      !!
      !! **`E_IJ^pol` of the EFMO energy**, and the reason it exists: every
      !! quantum-mechanical dimer of eq 6 already holds the mutual induction of
      !! its two fragments, which `E_pol^total` -- the induction solved over
      !! every fragment at once -- would then count a second time. So the pair
      !! term is subtracted from each near dimer, and for that subtraction to be
      !! clean it has to be *the same quantity* the total is built from.
      !!
      !! Which is why this is the same solver on a two-fragment system rather
      !! than a pair formula: same static field truncated at the quadrupole, the
      !! same screening `build_efp_system` puts on the points, and the same
      !! iteration to the same default tolerance, because no optional is passed
      !! here and none is passed by `efp_interaction_energy` either. A pair
      !! solved with a different convergence or a different field rank would
      !! leave a residue in the total that looks like three-body induction.
      !!
      !! Cheap: two fragments carry a couple of dozen polarizable points.
      type(efp_fragment_t), intent(in) :: frag_i, frag_j
      real(dp), intent(in) :: translation_i(3), translation_j(3)   !! Bohr
      type(error_t), intent(inout) :: error
      real(dp) :: energy

      type(efp_fragment_t) :: pair(2)
      type(efp_system_t) :: system
      real(dp) :: shifts(3, 2)

      energy = 0.0_dp
      pair(1) = frag_i
      pair(2) = frag_j
      shifts(:, 1) = translation_i
      shifts(:, 2) = translation_j

      call build_efp_system(pair, shifts, system, error)
      if (error%has_error()) return
      energy = polarization_energy(system, pair, error)
      call system%destroy()
   end function pair_polarization_energy

   function efp_pair_terms(fragments, translations, pairs, error, charge_transfer_on) &
      result(terms)
      !! The EFP terms of a given list of pairs, one result per pair
      !!
      !! **The far half of the EFMO energy.** A pair beyond `R_cut` contributes
      !! `E_IJ^Coul + E_IJ^disp + E_IJ^ExRep + E_IJ^CT` and nothing else;
      !! induction is deliberately absent, being carried whole by the
      !! system-wide `E_pol^total`. `pairs(:, k)` names the two fragments of
      !! pair `k` as indices into `fragments`.
      !!
      !! **No `place_fragment`.** That routine finds the rigid transform between
      !! a deck's atoms and a potential's own geometry, which a potential made on
      !! the fly for the geometry it is used at does not need: the fragment is
      !! already in its working frame and the translation is whatever the caller
      !! passes, commonly zero.
      !!
      !! **Electrostatics is decomposed by building a two-fragment system per
      !! pair**, because `electrostatic_energy` works on a flattened point set
      !! and takes no pair mask. That is exact rather than an approximation: the
      !! energy is a sum over point pairs on *different* fragments, and the
      !! charge-penetration screening is itself pairwise, so a system of two
      !! fragments reproduces their contribution to the full sum term by term.
      !! The one thing it does not reproduce is the *absence* of screening: the
      !! full system switches penetration off entirely when any fragment lacks a
      !! `SCREEN2` block, where a pair of screened fragments keeps it.
      !!
      !! A term whose data a potential does not carry is left at zero, as in
      !! `efp_interaction_energy`.
      type(efp_fragment_t), intent(in) :: fragments(:)
      real(dp), intent(in) :: translations(:, :)   !! (3, n_fragments), Bohr
      integer, intent(in) :: pairs(:, :)           !! (2, n_pairs), into `fragments`
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: charge_transfer_on
         !! Include `E_IJ^CT`. Default true, which is what GAMESS's EFMO does;
         !! the original method left it out, so it is switchable.
      type(efp_pair_energy_t), allocatable :: terms(:)

      type(efp_fragment_t) :: pair(2)
      type(efp_system_t) :: system
      real(dp) :: shifts(3, 2)
      integer :: n, k, a, b
      logical :: have_dynamic, have_lmo, have_ct, want_ct

      n = size(pairs, 2)
      allocate (terms(n))
      if (size(pairs, 1) /= 2) then
         call error%set(ERROR_VALIDATION, "efp: a pair list must be (2, n_pairs)")
         return
      end if
      if (size(translations, 1) /= 3 .or. size(translations, 2) /= size(fragments)) then
         call error%set(ERROR_VALIDATION, "efp: one translation per fragment is "// &
                        "needed, as (3, n_fragments)")
         return
      end if
      want_ct = .true.
      if (present(charge_transfer_on)) want_ct = charge_transfer_on

      do k = 1, n
         a = pairs(1, k)
         b = pairs(2, k)
         if (a < 1 .or. b < 1 .or. a > size(fragments) .or. b > size(fragments) &
             .or. a == b) then
            call error%set(ERROR_VALIDATION, "efp: a pair names a fragment that is "// &
                           "not in the list, or names one twice")
            return
         end if
         terms(k)%i = a
         terms(k)%j = b

         pair(1) = fragments(a)
         pair(2) = fragments(b)
         shifts(:, 1) = translations(:, a)
         shifts(:, 2) = translations(:, b)
         call build_efp_system(pair, shifts, system, error)
         if (error%has_error()) return
         terms(k)%electrostatics = electrostatic_energy(system, MAX_RANK, screen=.true.)
         call system%destroy()

         have_dynamic = fragments(a)%has_dynamic .and. fragments(b)%has_dynamic
         have_lmo = fragments(a)%has_lmo .and. fragments(b)%has_lmo
         have_ct = fragments(a)%has_ctvec .and. fragments(b)%has_ctvec &
                   .and. fragments(a)%has_ctfok .and. fragments(b)%has_ctfok

         if (have_lmo) then
            terms(k)%exchange_repulsion = exchange_repulsion(fragments(a), fragments(b), &
                                                             translations(:, a), &
                                                             translations(:, b), error)
            if (error%has_error()) return
         end if

         if (have_dynamic .and. have_lmo) then
            terms(k)%dispersion_e6 = dispersion_e6_damped(fragments(a), fragments(b), &
                                                          translations(:, a), &
                                                          translations(:, b), error)
            if (error%has_error()) return
            if (fragments(a)%has_dipquad .and. fragments(b)%has_dipquad) then
               terms(k)%dispersion_e7 = dispersion_e7_damped(fragments(a), fragments(b), &
                                                             translations(:, a), &
                                                             translations(:, b), error)
               if (error%has_error()) return
            end if
            if (fragments(a)%has_quadquad .and. fragments(b)%has_quadquad) then
               terms(k)%dispersion_e8 = dispersion_e8_damped(fragments(a), fragments(b), &
                                                             translations(:, a), &
                                                             translations(:, b), error)
               if (error%has_error()) return
            end if
         end if

         if (have_ct .and. want_ct) then
            terms(k)%charge_transfer = charge_transfer(fragments(a), fragments(b), &
                                                       translations(:, a), &
                                                       translations(:, b), error)
            if (error%has_error()) return
         end if

         terms(k)%dispersion = terms(k)%dispersion_e6 + terms(k)%dispersion_e7 &
                               + terms(k)%dispersion_e8
         terms(k)%total = terms(k)%electrostatics + terms(k)%exchange_repulsion &
                          + terms(k)%dispersion + terms(k)%charge_transfer
      end do
   end function efp_pair_terms

   subroutine place_fragment(frag, coords, rot, translation, error)
      !! Where a deck's atoms put a fragment: a rotation and a shift
      !!
      !! A potential carries the geometry it was made at; a deck names atoms in its
      !! own coordinates. Placing the fragment is finding the rigid transform between
      !! the two, `deck = rot . own + translation`, which `superpose` builds from the
      !! frame three atoms define.
      !!
      !! **The rotation is returned, not applied.** Turning the fragment is
      !! `rotate_fragment`, a separate step because it rewrites every tensor the
      !! potential carries.
      !!
      !! Every atom is then required to land where the potential says, which is a
      !! test rather than a fit: an EFP fragment is rigid, so a residual means the
      !! deck and the potential disagree -- the wrong species, atoms listed in a
      !! different order, or a geometry relaxed since the potential was made.
      type(efp_fragment_t), intent(in) :: frag
      real(dp), intent(in) :: coords(:, :)     !! (3, n_atoms) from the deck, Bohr
      real(dp), intent(out) :: rot(3, 3)
      real(dp), intent(out) :: translation(3)
      type(error_t), intent(inout) :: error

      real(dp) :: rmsd
      integer :: n

      n = size(coords, 2)
      if (size(coords, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "efp: fragment coordinates must be (3, n)")
         return
      end if
      if (n /= frag%n_atoms) then
         call error%set(ERROR_VALIDATION, "efp: this fragment has "//count_text(n)// &
                        " atoms in the deck but its potential describes "// &
                        count_text(frag%n_atoms))
         return
      end if

      call superpose(frag%points(:, 1:n), coords, rot, translation, rmsd, error)
      if (error%has_error()) return

      if (rmsd > PLACEMENT_TOL) then
         call error%set(ERROR_VALIDATION, "efp: this fragment is not a rigid "// &
                        "placement of its potential's own geometry -- the atoms miss "// &
                        "by "//real_text(rmsd)//" Bohr after the best rigid "// &
                        "transform. Check that the potential describes this species, "// &
                        "that the atoms are listed in the potential's own order, and "// &
                        "that the geometry is the one the potential was built at.")
         return
      end if
   end subroutine place_fragment

   pure function real_text(x) result(text)
      !! A number as text, for a message
      ! TODO(mqc): this and `count_text` below duplicate `pic_io`'s `to_char`,
      ! which every other module in this backend uses for the same job.
      real(dp), intent(in) :: x
      character(len=:), allocatable :: text

      character(len=24) :: buffer

      write (buffer, "(ES10.3)") x
      text = trim(adjustl(buffer))
   end function real_text

   pure function count_text(n) result(text)
      !! A small integer as text, for a message
      integer, intent(in) :: n
      character(len=:), allocatable :: text

      character(len=16) :: buffer

      write (buffer, "(I0)") n
      text = trim(buffer)
   end function count_text

end module mqc_czt_efp_energy
