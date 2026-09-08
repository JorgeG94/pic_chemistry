!! The separation EFMO measures, and the pair split it decides
module mqc_czt_efmo_pairs
   !! `R_IJ` of eq 2 of the EFMO paper (Sattasathuchana et al., JCTC 20, 2445
   !! (2024)) and the one decision it makes:
   !!
   !!     R_IJ = min_{i in I, j in J} |r_i - r_j| / (r_i^vdW + r_j^vdW)
   !!
   !! **It is unitless, not a distance.** Each interatomic distance is divided
   !! by the sum of the two van der Waals radii, so `R_IJ = 1` is contact and the
   !! default cutoff `R_cut = 2.0` is twice that -- a number that means the same
   !! thing for a water pair and for two aromatic rings, which an Angstrom
   !! threshold does not. Nothing here is comparable to
   !! `keywords.fragmentation.cutoffs`, which MBE uses and which is in Angstrom.
   !!
   !! What EFMO does with it: a pair at or inside the cutoff is computed as a
   !! quantum-mechanical dimer in vacuo, with its EFP pair induction subtracted;
   !! a pair beyond it gets the four EFP pair terms instead. Every pair is in
   !! exactly one of the two lists, so a pair counted in neither -- or in both --
   !! is a hole in the energy that no term reports.
   !!
   !! `vdw_scaled_distance` here is the kernel `mqc_czt_fmo`'s
   !! `unitless_distance` is written in terms of, so FMO's point-charge
   !! approximation cutoff and EFMO's dimer split measure separation with one
   !! formula rather than two that can drift apart.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_atomic_radii, only: vdw_radius_fmo
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   use pic_io, only: to_char
   implicit none
   private

   public :: vdw_scaled_distance
   public :: efmo_pair_distance
   public :: efmo_split_pairs

   type :: fragment_atoms_t
      !! One fragment's atoms, gathered out of the owner list
      integer, allocatable :: z(:)
      real(dp), allocatable :: xyz(:, :)   !! (3, n_atoms), Bohr
   end type fragment_atoms_t

contains

   pure function vdw_scaled_distance(z_a, xyz_a, z_b, xyz_b) result(r)
      !! One atom pair's separation in units of their van der Waals contact
      !!
      !! The radii are GAMESS's `$FMO VDWRAD` table and **not** Bondi's. FMO's
      !! `RESPPC` and `RESDIM`, and EFMO's `R_cut`, are all quoted in the
      !! literature against that table; with Bondi's hydrogen and oxygen the same
      !! water pair comes out 6 per cent further apart, which is enough to move a
      !! dimer across `R_cut = 2.0` and change which pairs are solved quantum
      !! mechanically. Every element has a radius here -- GAMESS substitutes 2.5
      !! Angstrom for the ones it does not name -- so the result is never
      !! negative and the caller has no undefined case to answer for.
      integer, intent(in) :: z_a, z_b
      real(dp), intent(in) :: xyz_a(3), xyz_b(3)   !! Bohr
      real(dp) :: r

      real(dp) :: scale

      scale = (vdw_radius_fmo(z_a) + vdw_radius_fmo(z_b))*ANGSTROM_TO_BOHR
      r = norm2(xyz_a - xyz_b)/scale
   end function vdw_scaled_distance

   subroutine efmo_pair_distance(z_i, xyz_i, z_j, xyz_j, r, error)
      !! `R_IJ`: the closest approach of two fragments, in contact units
      !!
      !! The minimum over every atom of `I` against every atom of `J`, each pair
      !! scaled by its own radii sum -- so the pair that decides it is the
      !! closest *relative to how big its atoms are*, which need not be the pair
      !! at the shortest distance.
      integer, intent(in) :: z_i(:), z_j(:)
      real(dp), intent(in) :: xyz_i(:, :), xyz_j(:, :)   !! (3, n_atoms), Bohr
      real(dp), intent(out) :: r
      type(error_t), intent(inout) :: error

      integer :: ia, ib
      real(dp) :: this

      r = huge(1.0_dp)
      if (size(z_i) < 1 .or. size(z_j) < 1) then
         call error%set(ERROR_VALIDATION, "efmo: a fragment with no atoms has no "// &
                        "separation from anything")
         return
      end if
      if (size(xyz_i, 1) /= 3 .or. size(xyz_j, 1) /= 3 .or. &
          size(xyz_i, 2) /= size(z_i) .or. size(xyz_j, 2) /= size(z_j)) then
         call error%set(ERROR_VALIDATION, "efmo: each fragment needs one (3, n_atoms) "// &
                        "coordinate array matching its atomic numbers")
         return
      end if

      do ia = 1, size(z_i)
         do ib = 1, size(z_j)
            this = vdw_scaled_distance(z_i(ia), xyz_i(:, ia), z_j(ib), xyz_j(:, ib))
            if (this < 0.0_dp) then
               call error%set(ERROR_VALIDATION, "efmo: no van der Waals radius for "// &
                              "element "//to_char(z_i(ia))//" or "//to_char(z_j(ib))// &
                              ", so the separation the dimer cutoff is measured in "// &
                              "is undefined")
               return
            end if
            r = min(r, this)
         end do
      end do
   end subroutine efmo_pair_distance

   subroutine efmo_split_pairs(owner, z, xyz, rcut, qm_pairs, efp_pairs, error)
      !! Every fragment pair, sorted into the quantum-mechanical and the EFP list
      !!
      !! `owner(i)` is the fragment atom `i` belongs to, numbered from one and
      !! contiguous -- the same `owner(atom)` partition `run_fmo2` takes. The two
      !! lists come back as `(2, n_pairs)` with `I < J`, and together they hold
      !! every unordered pair exactly once: a pair is quantum-mechanical when
      !! `R_IJ <= rcut` and effective otherwise, which is eq 6's split.
      !!
      !! `rcut` at or below zero puts every pair in the EFP list and `rcut` huge
      !! puts every pair in the QM one; both are limits worth being able to run,
      !! since the first is EFP with in-vacuo monomers and the second is FMO2 in
      !! vacuo plus the induction correction.
      integer, intent(in) :: owner(:)          !! (n_atoms), fragment of each atom
      integer, intent(in) :: z(:)              !! (n_atoms)
      real(dp), intent(in) :: xyz(:, :)        !! (3, n_atoms), Bohr
      real(dp), intent(in) :: rcut             !! Unitless, the `R_cut` of eq 2
      integer, allocatable, intent(out) :: qm_pairs(:, :)    !! (2, n_qm)
      integer, allocatable, intent(out) :: efp_pairs(:, :)   !! (2, n_efp)
      type(error_t), intent(inout) :: error

      integer :: n_atoms, n_frag, i, j, k, n_qm, n_efp
      integer, allocatable :: count_of(:), at(:)
      integer, allocatable :: list(:, :)
      logical, allocatable :: near(:)
      real(dp) :: r
      type(fragment_atoms_t), allocatable :: frag(:)

      n_atoms = size(owner)
      if (size(z) /= n_atoms .or. size(xyz, 1) /= 3 .or. size(xyz, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "efmo: the owner list, the atomic numbers "// &
                        "and the coordinates must describe the same atoms")
         return
      end if
      if (n_atoms < 1) then
         call error%set(ERROR_VALIDATION, "efmo: there are no atoms to fragment")
         return
      end if
      n_frag = maxval(owner)
      if (minval(owner) < 1) then
         call error%set(ERROR_VALIDATION, "efmo: every atom must belong to a fragment "// &
                        "numbered from one")
         return
      end if

      ! Gathered per fragment once, rather than scanned per pair: the pair loop
      ! is quadratic in the fragment count and would otherwise be quadratic in
      ! the atom count as well.
      allocate (count_of(n_frag), at(n_frag), frag(n_frag))
      count_of = 0
      do i = 1, n_atoms
         count_of(owner(i)) = count_of(owner(i)) + 1
      end do
      if (any(count_of == 0)) then
         call error%set(ERROR_VALIDATION, "efmo: fragment numbering has a gap -- some "// &
                        "fragment between one and "//to_char(n_frag)//" holds no atoms")
         return
      end if
      do k = 1, n_frag
         allocate (frag(k)%z(count_of(k)), frag(k)%xyz(3, count_of(k)))
      end do
      at = 0
      do i = 1, n_atoms
         k = owner(i)
         at(k) = at(k) + 1
         frag(k)%z(at(k)) = z(i)
         frag(k)%xyz(:, at(k)) = xyz(:, i)
      end do

      allocate (list(2, n_frag*(n_frag - 1)/2))
      allocate (near(n_frag*(n_frag - 1)/2))
      k = 0
      do i = 1, n_frag - 1
         do j = i + 1, n_frag
            call efmo_pair_distance(frag(i)%z, frag(i)%xyz, frag(j)%z, frag(j)%xyz, &
                                    r, error)
            if (error%has_error()) return
            k = k + 1
            list(:, k) = [i, j]
            ! At the cutoff exactly the pair is quantum-mechanical: eq 6 writes
            ! the near sum as `R_IJ <= R_cut`.
            near(k) = r <= rcut
         end do
      end do

      n_qm = count(near)
      n_efp = k - n_qm
      allocate (qm_pairs(2, n_qm), efp_pairs(2, n_efp))
      n_qm = 0
      n_efp = 0
      do i = 1, k
         if (near(i)) then
            n_qm = n_qm + 1
            qm_pairs(:, n_qm) = list(:, i)
         else
            n_efp = n_efp + 1
            efp_pairs(:, n_efp) = list(:, i)
         end if
      end do
   end subroutine efmo_split_pairs

end module mqc_czt_efmo_pairs
