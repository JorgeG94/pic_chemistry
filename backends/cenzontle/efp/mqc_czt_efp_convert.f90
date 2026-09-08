!! A fragment potential held in memory, as the fragment the energy terms take
module mqc_czt_efp_convert
   !! The missing link between MAKEFP and the EFP energy: `efp_potential_t`, what
   !! `make_efp_potential` builds, turned into `efp_fragment_t`, what every
   !! interaction term consumes.
   !!
   !! Today the only producer of an `efp_fragment_t` is `read_efp_potential`, so
   !! a potential computed in this process has to be written to a `.efp` file and
   !! read back to be used. **EFMO makes a potential per fragment on the fly**,
   !! so that round trip would be one file write and one parse per monomer of
   !! every run, and every number would arrive truncated to the ten decimals the
   !! format carries.
   !!
   !! **The contract is that this equals the round trip.** Same fields, same
   !! `has_*` flags, same orderings -- the localized orbitals and `CTVEC` stay in
   !! GAMESS's AO order, as the reader leaves them, because `mqc_czt_efp_pair`
   !! inverts that order itself against a molecule it builds. Anything here that
   !! disagreed with the reader would be a term computed from data in a layout
   !! the term does not expect, which is silent. `test_mqc_czt_efmo_pieces`
   !! writes, reads and compares field by field.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_efp_potential, only: efp_potential_t
   use mqc_czt_efp_read, only: efp_fragment_t, read_projection_basis
   implicit none
   private

   public :: potential_to_fragment

   integer, parameter :: N_DIPQUAD = 27
   integer, parameter :: N_QUADQUAD = 81
      !! Values a dipole-quadrupole and a quadrupole-quadrupole record carries,
      !! the same counts the reader expects of the file's own records.

   integer, parameter :: MAX_LINE = 160
      !! The width the potential carries its projection basis lines at, matching
      !! the reader's.

contains

   subroutine potential_to_fragment(pot, frag, error)
      !! One computed potential as a placed-fragment description
      !!
      !! Field for field what `write_efp_potential` followed by
      !! `read_efp_potential` produces, without the file: the same blocks in the
      !! same layouts, at full precision rather than the format's ten decimals.
      !!
      !! **The two flat tensor blocks are packed here, not copied.** The
      !! potential holds `dipquad` and `quadquad` shaped `(3,3,3,...)` and
      !! `(3,3,3,3,...)`; a fragment holds them flat, because which slot is which
      !! index tuple is a convention of the file rather than something the reader
      !! can know. The packing below is the writer's, so a change to one is a
      !! change to both.
      !!
      !! **The projection basis goes through the reader.** The potential carries
      !! it as the text lines the writer emits -- GAMESS's own columns, its `L`
      !! shells, its primitive normalization folded into every coefficient -- and
      !! undoing all of that is exactly `read_projection_basis`, so it is called
      !! on those lines rather than reimplemented.
      type(efp_potential_t), intent(in) :: pot
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: error

      character(len=MAX_LINE), allocatable :: basis_lines(:)
      integer :: i, k, f, a, b, c, e, slot

      if (pot%n_points < 1) then
         call error%set(ERROR_VALIDATION, "efp: this potential carries no "// &
                        "expansion points, so there is nothing to interact")
         return
      end if

      if (allocated(pot%name)) frag%name = trim(pot%name)
      frag%multiplicity = pot%multiplicity

      ! --- the points and their multipoles --------------------------------------
      frag%n_points = pot%n_points
      allocate (frag%labels(frag%n_points), frag%points(3, frag%n_points), &
                frag%mass(frag%n_points), frag%charge(frag%n_points))
      do i = 1, frag%n_points
         frag%labels(i) = trim(adjustl(pot%labels(i)))
      end do
      frag%points = pot%points
      frag%mass = pot%mass
      frag%charge = pot%charge
      ! An atom is a point carrying a nuclear charge; a bond midpoint carries
      ! none. The reader tells them apart the same way rather than by counting.
      frag%n_atoms = count(frag%charge > 0.0_dp)

      allocate (frag%q_elec(frag%n_points), frag%q_nuc(frag%n_points))
      frag%q_elec = pot%q_elec
      frag%q_nuc = pot%q_nuc
      frag%dipole = pot%dipole
      frag%quadrupole = pot%quadrupole
      frag%octopole = pot%octopole

      ! --- charge-penetration screening -----------------------------------------
      ! Both damping forms are optional in the format; a potential that carries
      ! them sets the flags the electrostatics reads.
      if (allocated(pot%screen2)) then
         frag%screen2 = pot%screen2
         frag%has_screen2 = .true.
      end if
      if (allocated(pot%screen)) then
         frag%screen = pot%screen
         frag%has_screen = .true.
      end if

      ! --- polarization ---------------------------------------------------------
      ! The static polarizabilities sit at the localized-orbital centroids, which
      ! is why the fragment counts them separately from the dynamic set even
      ! though the two are the same points here.
      if (allocated(pot%static_pol)) then
         frag%n_pol = pot%n_lmo
         frag%static_pol = pot%static_pol
         frag%pol_points = pot%centroids
         frag%has_static_pol = .true.
      end if

      ! --- dispersion -----------------------------------------------------------
      if (allocated(pot%dynamic_pol) .and. allocated(pot%frequencies)) then
         frag%n_lmo = pot%n_lmo
         frag%n_freq = size(pot%frequencies)
         frag%dyn_pol = pot%dynamic_pol
         frag%centroids = pot%centroids
         frag%frequencies = pot%frequencies
         frag%has_dynamic = .true.
      end if

      if (frag%has_dynamic .and. allocated(pot%dipquad)) then
         frag%n_dipquad = N_DIPQUAD
         allocate (frag%dipquad(N_DIPQUAD, frag%n_lmo, frag%n_freq))
         do f = 1, frag%n_freq
            do k = 1, frag%n_lmo
               do a = 1, 3
                  do b = 1, 3
                     do c = 1, 3
                        ! The *first* quadrupole index runs fastest, which is
                        ! the writer's slot order and transposed from how the
                        ! `DQSHIFT` source reads.
                        frag%dipquad((a - 1)*9 + (c - 1)*3 + b, k, f) = &
                           pot%dipquad(a, b, c, k, f)
                     end do
                  end do
               end do
            end do
         end do
         frag%has_dipquad = .true.
      end if

      if (frag%has_dynamic .and. allocated(pot%quadquad)) then
         frag%n_quadquad = N_QUADQUAD
         allocate (frag%quadquad(N_QUADQUAD, frag%n_lmo, frag%n_freq))
         do f = 1, frag%n_freq
            do k = 1, frag%n_lmo
               slot = 0
               do a = 1, 3
                  do b = 1, 3
                     do c = 1, 3
                        do e = 1, 3
                           ! Last index fastest, no transposition: every
                           ! `QQSHIFT` term is symmetric within each index pair.
                           slot = slot + 1
                           frag%quadquad(slot, k, f) = pot%quadquad(a, b, c, e, k, f)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         frag%has_quadquad = .true.
      end if

      ! --- exchange repulsion and charge transfer -------------------------------
      ! Left in GAMESS's AO order, as the reader leaves them: converting needs
      ! the shell layout of a built molecule, which `mqc_czt_efp_pair` has and
      ! this does not.
      if (allocated(pot%orbitals)) then
         frag%n_lmo_proj = pot%n_lmo
         frag%nao_proj = pot%nao
         frag%lmo_gamess = pot%orbitals
         frag%has_lmo = .true.
      end if

      if (allocated(pot%fock_lmo)) then
         ! The file carries the lower triangle and the reader unpacks it
         ! symmetric, so what the fragment holds is the symmetric part of this.
         frag%fock_lmo = pot%fock_lmo
         frag%has_fock = .true.
      end if

      if (allocated(pot%canonical)) then
         frag%n_occ_ct = pot%n_occ
         frag%n_mo_ct = pot%nao
         frag%ctvec_gamess = pot%canonical
         frag%has_ctvec = .true.
         if (allocated(pot%eps_occ)) then
            frag%eps_occ = pot%eps_occ(1:pot%n_occ)
            frag%has_ctfok = .true.
         end if
      end if

      ! --- the projection basis --------------------------------------------------
      if (allocated(pot%basis_lines)) then
         allocate (basis_lines(size(pot%basis_lines) + 2))
         basis_lines(1) = " PROJECTION BASIS SET"
         do i = 1, size(pot%basis_lines)
            basis_lines(i + 1) = pot%basis_lines(i)
         end do
         basis_lines(size(basis_lines)) = " STOP"
         call read_projection_basis(basis_lines, size(basis_lines), frag, error)
         deallocate (basis_lines)
         if (error%has_error()) return
      end if
   end subroutine potential_to_fragment

end module mqc_czt_efp_convert
