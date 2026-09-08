!! Read a `.efp` fragment potential back in
module mqc_czt_efp_read
   !! The other half of `mqc_czt_efp_potential`: that module computes a potential and
   !! writes it, this one reads one and hands back the parameters an interaction
   !! energy needs.
   !!
   !! Sections are found and skipped by name rather than by position, so adding
   !! one later is a matter of naming it.
   !!
   !! **The format is GAMESS's, so it is read the way GAMESS writes it** rather
   !! than the way a reader would prefer:
   !!
   !!   * a section is a header line, records, then `STOP`
   !!   * a record is a label followed by numbers, and a record too long for a line
   !!     is continued with a trailing `>`
   !!   * a header is distinguished from a record by starting in column two and
   !!     carrying no numbers of its own -- except that some do (`MULTIPLICITY 1`,
   !!     `PROJECTION WAVEFUNCTION 4 19`), which is why the known names are matched
   !!     rather than the shape guessed at
   !!
   !! Units are as the file carries them: Bohr for the points, atomic units for the
   !! multipoles. Nothing is converted on the way in.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_IO
   use mqc_czt_efp_potential, only: gamess_primitive_norm
   implicit none
   private

   public :: efp_fragment_t
   public :: read_efp_potential
   ! Exposed so a potential held in memory can be turned into a fragment
   ! without a file: the projection basis is the one block a converter cannot
   ! copy field for field, because the potential carries it as the text lines
   ! the writer emits and only this routine knows how to read them.
   public :: read_projection_basis
   ! How many components each stored multipole carries. Exposed because anything
   ! that unpacks one needs the same count, and two files agreeing by coincidence
   ! is how a packing convention drifts.
   public :: N_DIPOLE, N_QUADRUPOLE, N_OCTUPOLE

   integer, parameter :: N_DIPOLE = 3
   integer, parameter :: N_QUADRUPOLE = 6
   integer, parameter :: N_OCTUPOLE = 10
      !! Multipole components at one expansion point, in the file's own order.

   integer, parameter :: N_DYNAMIC_RECORD = 12
      !! A dynamic polarizability record: a centroid then a 3x3 tensor.

   ! Row and column of each of GAMESS's nine polarizability slots. The diagonal
   ! comes first and the off-diagonal triples are the transpose of what the
   ! labels suggest -- the same map `mqc_czt_efp_potential` writes with. Reading the
   ! nine as a row-major 3x3 gives a tensor whose trace is negative, which is
   ! how the mistake announces itself.
   !
   ! **`efinp.src` appears to contradict this, and does not.** Its reader
   ! (`efinp.src:8456-8464`) and its writer (`efinp.src:7552-7561`) agree that
   ! slot 4 is `DYNDD_LMO(1,2)`, the transpose of the map below. Both are true
   ! because the two arrays are indexed oppositely: GAMESS's is
   ! `alpha(field, dipole)` and ours is `alpha(dipole, field)`. The maps compose
   ! to the identity, so the *files* agree. Do not "fix" this to match
   ! `efinp.src`.
   !
   ! What does follow from the transpose is that anything mirroring a GAMESS
   ! expression over `DYNDD_LMO` must swap the two indices when it reads
   ! `dyn_pol`. `dipquad` and `quadquad` are stored flat in file order and need
   ! no swap.
   integer, parameter :: N_POL_SLOTS = 9
   integer, parameter :: POL_ROW(N_POL_SLOTS) = [1, 2, 3, 2, 3, 3, 1, 1, 2]
   integer, parameter :: POL_COL(N_POL_SLOTS) = [1, 2, 3, 1, 1, 2, 2, 3, 3]

   integer, parameter :: MAX_LINE = 160
      !! Longest line a potential is expected to carry, matching the writer's.

   type :: efp_fragment_t
      !! One fragment's parameters, as read
      character(len=:), allocatable :: name        !! From `$FRAGNAME`, without the `$`
      integer :: n_points = 0                      !! Atoms, then bond midpoints
      integer :: n_atoms = 0                       !! Points carrying a nuclear charge
      integer :: multiplicity = 1
      character(len=8), allocatable :: labels(:)
      real(dp), allocatable :: points(:, :)        !! (3, n_points), Bohr
      real(dp), allocatable :: mass(:)             !! amu, zero at a midpoint
      real(dp), allocatable :: charge(:)           !! Z, zero at a midpoint
      real(dp), allocatable :: q_elec(:)           !! Electronic monopole
      real(dp), allocatable :: q_nuc(:)            !! Nuclear monopole
      real(dp), allocatable :: dipole(:, :)        !! (3, n_points)
      real(dp), allocatable :: quadrupole(:, :)    !! (6, n_points)
      real(dp), allocatable :: octopole(:, :)      !! (10, n_points)
      real(dp), allocatable :: screen(:)           !! Gaussian damping exponent per point
      real(dp), allocatable :: screen2(:)          !! Exponential damping exponent per point
      logical :: has_screen = .false.
      logical :: has_screen2 = .false.

      integer :: n_lmo = 0
      integer :: n_freq = 0
      real(dp), allocatable :: dyn_pol(:, :, :, :)
         !! Dynamic polarizabilities, for dispersion. `(3, 3, n_lmo, n_freq)`.
      real(dp), allocatable :: centroids(:, :)     !! (3, n_lmo), Bohr
      real(dp), allocatable :: frequencies(:)      !! Imaginary, a.u.
      logical :: has_dynamic = .false.

      integer :: n_pol = 0
      real(dp), allocatable :: static_pol(:, :, :)  !! (3, 3, n_pol)
         !! Static polarizabilities at the same centroids, for polarization.
      real(dp), allocatable :: pol_points(:, :)     !! (3, n_pol), Bohr
      logical :: has_static_pol = .false.

      integer :: n_dipquad = 0                    !! Values per record, expected 27
      integer :: n_quadquad = 0                   !! Values per record, expected 81
      real(dp), allocatable :: dipquad(:, :, :)   !! (n_dipquad, n_lmo, n_freq)
      real(dp), allocatable :: quadquad(:, :, :)  !! (n_quadquad, n_lmo, n_freq)
         !! The two higher dispersion tensor sets, **stored flat**: which of the
         !! 27 or 81 slots is which index tuple is not established, and shaping
         !! them would bake in a guess.
      logical :: has_dipquad = .false.
      logical :: has_quadquad = .false.

      integer :: n_shells = 0
      integer, allocatable :: shell_atom(:), shell_l(:), shell_first(:), shell_nprim(:)
      real(dp), allocatable :: prim_expo(:), prim_coef(:)
         !! The projection basis, with GAMESS's primitive normalization divided
         !! back out, so these are the raw contraction coefficients a basis
         !! object wants. An `L` shell arrives as two shells over shared
         !! exponents.
      logical :: has_basis = .false.

      integer :: n_lmo_proj = 0
      integer :: nao_proj = 0
      real(dp), allocatable :: lmo_gamess(:, :)   !! (nao_proj, n_lmo_proj)
         !! The localized orbitals, in **GAMESS's** AO order, not libcint's.
         !! Converting them needs the shell layout of a built molecule, so it is
         !! left to whoever builds one; see `mqc_czt_efp_pair`.
      logical :: has_lmo = .false.

      real(dp), allocatable :: fock_lmo(:, :)     !! (n_lmo_proj, n_lmo_proj)
         !! The Fock operator over the localized orbitals, unpacked to a full
         !! symmetric matrix. The file stores its lower triangle row by row and
         !! carries no labels.
      logical :: has_fock = .false.

      integer :: n_occ_ct = 0
      integer :: n_mo_ct = 0
      real(dp), allocatable :: ctvec_gamess(:, :)  !! (nao_proj, n_mo_ct)
         !! `CTVEC`, which charge transfer needs. A wider orbital set than the
         !! projection wavefunction's -- occupied *and* virtual, since charge
         !! transfer moves density into the latter -- and in GAMESS's AO order,
         !! so it needs the same inversion.
      real(dp), allocatable :: eps_occ(:)          !! (n_occ_ct), from CTFOK
      logical :: has_ctvec = .false.
      logical :: has_ctfok = .false.
   contains
      procedure :: destroy => fragment_destroy
      procedure :: net_charge => fragment_net_charge
   end type efp_fragment_t

contains

   subroutine fragment_destroy(self)
      class(efp_fragment_t), intent(inout) :: self

      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%labels)) deallocate (self%labels)
      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%mass)) deallocate (self%mass)
      if (allocated(self%charge)) deallocate (self%charge)
      if (allocated(self%q_elec)) deallocate (self%q_elec)
      if (allocated(self%q_nuc)) deallocate (self%q_nuc)
      if (allocated(self%dipole)) deallocate (self%dipole)
      if (allocated(self%quadrupole)) deallocate (self%quadrupole)
      if (allocated(self%octopole)) deallocate (self%octopole)
      if (allocated(self%screen)) deallocate (self%screen)
      if (allocated(self%screen2)) deallocate (self%screen2)
      if (allocated(self%dyn_pol)) deallocate (self%dyn_pol)
      if (allocated(self%centroids)) deallocate (self%centroids)
      if (allocated(self%frequencies)) deallocate (self%frequencies)
      self%n_lmo = 0
      self%n_freq = 0
      self%has_dynamic = .false.
      if (allocated(self%static_pol)) deallocate (self%static_pol)
      if (allocated(self%pol_points)) deallocate (self%pol_points)
      self%n_pol = 0
      self%has_static_pol = .false.
      if (allocated(self%dipquad)) deallocate (self%dipquad)
      if (allocated(self%quadquad)) deallocate (self%quadquad)
      self%n_dipquad = 0
      self%n_quadquad = 0
      self%has_dipquad = .false.
      self%has_quadquad = .false.
      if (allocated(self%shell_atom)) deallocate (self%shell_atom)
      if (allocated(self%shell_l)) deallocate (self%shell_l)
      if (allocated(self%shell_first)) deallocate (self%shell_first)
      if (allocated(self%shell_nprim)) deallocate (self%shell_nprim)
      if (allocated(self%prim_expo)) deallocate (self%prim_expo)
      if (allocated(self%prim_coef)) deallocate (self%prim_coef)
      self%n_shells = 0
      self%has_basis = .false.
      if (allocated(self%lmo_gamess)) deallocate (self%lmo_gamess)
      self%n_lmo_proj = 0
      self%nao_proj = 0
      self%has_lmo = .false.
      if (allocated(self%fock_lmo)) deallocate (self%fock_lmo)
      self%has_fock = .false.
      if (allocated(self%ctvec_gamess)) deallocate (self%ctvec_gamess)
      if (allocated(self%eps_occ)) deallocate (self%eps_occ)
      self%n_occ_ct = 0
      self%n_mo_ct = 0
      self%has_ctvec = .false.
      self%has_ctfok = .false.
      self%n_points = 0
      self%n_atoms = 0
      self%has_screen = .false.
      self%has_screen2 = .false.
   end subroutine fragment_destroy

   pure function fragment_net_charge(self) result(q)
      !! The fragment's net charge, nuclear plus electronic
      !!
      !! The cheapest check that a potential was read correctly: the monopoles
      !! must sum to the charge the fragment was computed for.
      class(efp_fragment_t), intent(in) :: self
      real(dp) :: q

      q = 0.0_dp
      if (allocated(self%q_elec)) q = q + sum(self%q_elec)
      if (allocated(self%q_nuc)) q = q + sum(self%q_nuc)
   end function fragment_net_charge

   subroutine read_efp_potential(path, frag, error)
      !! Read the electrostatic parameters out of a `.efp` file
      character(len=*), intent(in) :: path
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: error

      character(len=MAX_LINE), allocatable :: lines(:)
      character(len=MAX_LINE), allocatable :: labels(:)
      real(dp), allocatable :: values(:, :)
      integer :: n_lines, i, n_rec
      character(len=:), allocatable :: section

      call slurp(path, lines, n_lines, error)
      if (error%has_error()) return

      call fragment_name(lines, n_lines, frag%name)

      ! COORDINATES first: it fixes how many points every other section must carry,
      ! so a section that disagrees is caught rather than silently truncated.
      call section_records(lines, n_lines, "COORDINATES", 5, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec == 0) then
         call error%set(ERROR_VALIDATION, "efp: '"//trim(path)// &
                        "' carries no COORDINATES section")
         return
      end if
      frag%n_points = n_rec
      allocate (frag%labels(n_rec), frag%points(3, n_rec), frag%mass(n_rec), &
                frag%charge(n_rec))
      do i = 1, n_rec
         frag%labels(i) = trim(labels(i))
         frag%points(:, i) = values(1:3, i)
         frag%mass(i) = values(4, i)
         frag%charge(i) = values(5, i)
      end do
      frag%n_atoms = count(frag%charge > 0.0_dp)

      call one_section(lines, n_lines, "MONOPOLES", 2, frag%n_points, values, error)
      if (error%has_error()) return
      allocate (frag%q_elec(frag%n_points), frag%q_nuc(frag%n_points))
      frag%q_elec = values(1, :)
      frag%q_nuc = values(2, :)

      call one_section(lines, n_lines, "DIPOLES", N_DIPOLE, frag%n_points, values, error)
      if (error%has_error()) return
      allocate (frag%dipole(N_DIPOLE, frag%n_points))
      frag%dipole = values

      call one_section(lines, n_lines, "QUADRUPOLES", N_QUADRUPOLE, frag%n_points, &
                       values, error)
      if (error%has_error()) return
      allocate (frag%quadrupole(N_QUADRUPOLE, frag%n_points))
      frag%quadrupole = values

      call one_section(lines, n_lines, "OCTUPOLES", N_OCTUPOLE, frag%n_points, &
                       values, error)
      if (error%has_error()) return
      allocate (frag%octopole(N_OCTUPOLE, frag%n_points))
      frag%octopole = values

      ! The two screening sections are optional; without them there is no
      ! penetration term. Each record is `1.0 alpha`, the leading one being a
      ! weight GAMESS always writes as unity.
      call section_records(lines, n_lines, "SCREEN2", 2, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec > 0) then
         call expect_points(n_rec, frag%n_points, "SCREEN2", error)
         if (error%has_error()) return
         allocate (frag%screen2(frag%n_points))
         frag%screen2 = values(2, :)
         frag%has_screen2 = .true.
      end if

      call section_records(lines, n_lines, "SCREEN", 2, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec > 0) then
         call expect_points(n_rec, frag%n_points, "SCREEN", error)
         if (error%has_error()) return
         allocate (frag%screen(frag%n_points))
         frag%screen = values(2, :)
         frag%has_screen = .true.
      end if

      call read_static_pol(lines, n_lines, frag, error)
      if (error%has_error()) return

      call read_dynamic(lines, n_lines, frag, error)
      if (error%has_error()) return

      call read_tensor_block(lines, n_lines, "DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS", &
                             27, frag%dipquad, frag%n_dipquad, frag%has_dipquad, &
                             frag%n_lmo, frag%n_freq, error)
      if (error%has_error()) return
      call read_tensor_block(lines, n_lines, "LMOQQPOL DYNAMIC POLARIZABLE POINTS", &
                             81, frag%quadquad, frag%n_quadquad, frag%has_quadquad, &
                             frag%n_lmo, frag%n_freq, error)
      if (error%has_error()) return

      call read_projection_wavefunction(lines, n_lines, frag, error)
      if (error%has_error()) return

      call read_fock_lmo(lines, n_lines, frag, error)
      if (error%has_error()) return

      call read_ctvec(lines, n_lines, frag, error)
      if (error%has_error()) return

      call read_projection_basis(lines, n_lines, frag, error)
      if (error%has_error()) return

      call multiplicity(lines, n_lines, frag%multiplicity)
      deallocate (lines)
      if (allocated(labels)) deallocate (labels)
      if (allocated(values)) deallocate (values)
   end subroutine read_efp_potential

   subroutine read_ctvec(lines, n_lines, frag, error)
      !! `CTVEC` and `CTFOK`: the orbitals and orbital energies charge transfer needs
      !!
      !! `CTVEC`'s header carries two counts, occupied and total -- `CTVEC 5 19`
      !! is five occupied among nineteen. **Not the same set as `PROJECTION
      !! WAVEFUNCTION`**, which is valence-occupied only; which of the whole
      !! canonical set or occupied-plus-valence-virtuals arrived depends on
      !! GAMESS's `CTVVO` and is what the header count says.
      !!
      !! Rows are chunked by column, as in `PROJECTION WAVEFUNCTION`, and the
      !! coefficients are in GAMESS's AO order.
      !!
      !! `CTFOK` follows as a bare unlabelled list of the occupied orbital
      !! energies, and is read here because it is meaningless without the count
      !! `CTVEC` declares.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer, parameter :: WIDTH = 15
      integer :: at, i, stat, col, filled, pos, n_occ, n_mo
      character(len=MAX_LINE) :: text
      real(dp) :: value

      at = 0
      do i = 1, n_lines
         if (index(adjustl(lines(i)), "CTVEC") == 1) then
            at = i
            exit
         end if
      end do
      if (at == 0) return
      text = adjustl(lines(at))
      read (text(6:), *, iostat=stat) n_occ, n_mo
      if (stat /= 0 .or. n_occ < 1 .or. n_mo < n_occ) then
         call error%set(ERROR_VALIDATION, "efp: the CTVEC header does not carry an "// &
                        "occupied and a total orbital count")
         return
      end if
      if (frag%nao_proj == 0) return

      frag%n_occ_ct = n_occ
      frag%n_mo_ct = n_mo
      allocate (frag%ctvec_gamess(frag%nao_proj, n_mo))
      frag%ctvec_gamess = 0.0_dp

      col = 0
      filled = 0
      do i = at + 1, n_lines
         text = lines(i)
         if (len_trim(text) == 0) cycle
         if (trim(adjustl(text)) == "STOP") exit
         read (text(1:2), *, iostat=stat) col
         if (stat /= 0 .or. col < 1 .or. col > n_mo) exit
         if (filled >= frag%nao_proj) filled = 0
         pos = 6
         do while (pos + WIDTH - 1 <= len_trim(text))
            read (text(pos:pos + WIDTH - 1), *, iostat=stat) value
            if (stat /= 0) exit
            filled = filled + 1
            if (filled > frag%nao_proj) then
               call error%set(ERROR_VALIDATION, "efp: CTVEC carries more coefficients "// &
                              "than there are basis functions")
               return
            end if
            frag%ctvec_gamess(filled, col) = value
            pos = pos + WIDTH
         end do
      end do
      frag%has_ctvec = .true.

      ! CTFOK: the occupied orbital energies, unlabelled.
      at = 0
      do i = 1, n_lines
         if (index(adjustl(lines(i)), "CTFOK") == 1) then
            at = i
            exit
         end if
      end do
      if (at == 0) return
      allocate (frag%eps_occ(n_occ))
      frag%eps_occ = 0.0_dp
      filled = 0
      do i = at + 1, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         rest_of_ctfok: block
            character(len=:), allocatable :: rest
            integer :: cut
            rest = lines(i)
            cut = index(rest, ">")
            if (cut > 0) rest = rest(1:cut - 1)
            rest = adjustl(rest)
            do
               if (len_trim(rest) == 0) exit
               call next_number(rest, value, stat)
               if (stat /= 0) exit
               filled = filled + 1
               if (filled > n_occ) exit
               frag%eps_occ(filled) = value
            end do
         end block rest_of_ctfok
         if (filled >= n_occ) exit
      end do
      frag%has_ctfok = filled == n_occ
   end subroutine read_ctvec

   subroutine read_fock_lmo(lines, n_lines, frag, error)
      !! `FOCK MATRIX ELEMENTS`: the Fock operator over the localized orbitals
      !!
      !! Exchange repulsion needs it, and it is the one section with no labels at
      !! all -- just the lower triangle, row by row, over continuation lines.
      !! Its size fixes the orbital count independently, so a mismatch against
      !! the projection wavefunction is caught rather than reshaped around.
      !!
      !! Unpacked to a full symmetric matrix here, so consumers index `(i, j)`
      !! without thinking about packing.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer :: start, finish, i, k, n_want, filled, row, col, stat
      real(dp), allocatable :: packed(:)
      character(len=:), allocatable :: rest
      character(len=MAX_LINE) :: text
      real(dp) :: value

      if (frag%n_lmo_proj == 0) return
      start = 0
      do i = 1, n_lines
         if (index(adjustl(lines(i)), "FOCK MATRIX ELEMENTS") == 1) then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return
      finish = start - 1
      do i = start, n_lines
         text = adjustl(lines(i))
         ! No labels here, so the section ends at `STOP` or at the next header -- a
         ! line that carries a letter where a number belongs.
         if (trim(text) == "STOP") exit
         if (len_trim(text) > 0) then
            if (index("0123456789-+.", text(1:1)) == 0) exit
         end if
         finish = i
      end do

      n_want = frag%n_lmo_proj*(frag%n_lmo_proj + 1)/2
      allocate (packed(n_want))
      packed = 0.0_dp
      filled = 0
      do i = start, finish
         rest = lines(i)
         k = index(rest, ">")
         if (k > 0) rest = rest(1:k - 1)
         rest = adjustl(rest)
         do
            if (len_trim(rest) == 0) exit
            call next_number(rest, value, stat)
            if (stat /= 0) exit
            filled = filled + 1
            if (filled > n_want) then
               call error%set(ERROR_VALIDATION, "efp: FOCK MATRIX ELEMENTS carries "// &
                              "more values than the orbital count allows")
               return
            end if
            packed(filled) = value
         end do
      end do
      if (filled /= n_want) then
         call error%set(ERROR_VALIDATION, "efp: FOCK MATRIX ELEMENTS is not a lower "// &
                        "triangle over the projection orbitals")
         return
      end if

      allocate (frag%fock_lmo(frag%n_lmo_proj, frag%n_lmo_proj))
      k = 0
      do row = 1, frag%n_lmo_proj
         do col = 1, row
            k = k + 1
            frag%fock_lmo(row, col) = packed(k)
            frag%fock_lmo(col, row) = packed(k)
         end do
      end do
      deallocate (packed)
      frag%has_fock = .true.
   end subroutine read_fock_lmo

   subroutine read_projection_wavefunction(lines, n_lines, frag, error)
      !! `PROJECTION WAVEFUNCTION`: the localized orbitals exchange repulsion needs
      !!
      !! The header carries the two counts -- `PROJECTION WAVEFUNCTION 4 19` is four
      !! localized orbitals over nineteen basis functions -- so nothing has to be
      !! inferred from the record count. Each row is an orbital index, a chunk index,
      !! then five coefficients in fixed columns, continuing over as many lines as the
      !! orbital needs.
      !!
      !! **Left in GAMESS's AO order.** Converting to libcint's needs the shell
      !! layout of a molecule, which this module does not build, and the d and f
      !! ordering maps in `mqc_czt_efp_potential` have to be applied in *reverse*.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer, parameter :: WIDTH = 15
      integer :: at, i, stat, lmo, filled, pos, n_lmo, nao
      character(len=MAX_LINE) :: text
      real(dp) :: value

      at = 0
      do i = 1, n_lines
         if (index(adjustl(lines(i)), "PROJECTION WAVEFUNCTION") == 1) then
            at = i
            exit
         end if
      end do
      if (at == 0) return

      text = adjustl(lines(at))
      read (text(24:), *, iostat=stat) n_lmo, nao
      if (stat /= 0 .or. n_lmo < 1 .or. nao < 1) then
         call error%set(ERROR_VALIDATION, "efp: the PROJECTION WAVEFUNCTION header "// &
                        "does not carry an orbital and a function count")
         return
      end if

      frag%n_lmo_proj = n_lmo
      frag%nao_proj = nao
      allocate (frag%lmo_gamess(nao, n_lmo))
      frag%lmo_gamess = 0.0_dp

      lmo = 0
      filled = 0
      do i = at + 1, n_lines
         text = lines(i)
         if (len_trim(text) == 0) cycle
         if (trim(adjustl(text)) == "STOP") exit
         ! A row belongs to the orbital its first token names, so a new orbital is
         ! recognised from the data rather than by counting lines per orbital -- the
         ! last row of each is short.
         read (text(1:2), *, iostat=stat) lmo
         if (stat /= 0 .or. lmo < 1 .or. lmo > n_lmo) exit
         if (lmo == 1 .and. filled >= nao) exit
         if (filled >= nao) filled = 0
         pos = 6
         do while (pos + WIDTH - 1 <= len_trim(text))
            read (text(pos:pos + WIDTH - 1), *, iostat=stat) value
            if (stat /= 0) exit
            filled = filled + 1
            if (filled > nao) then
               call error%set(ERROR_VALIDATION, "efp: PROJECTION WAVEFUNCTION carries "// &
                              "more coefficients than its header declares")
               return
            end if
            frag%lmo_gamess(filled, lmo) = value
            pos = pos + WIDTH
         end do
      end do
      frag%has_lmo = .true.
   end subroutine read_projection_wavefunction

   subroutine read_projection_basis(lines, n_lines, frag, error)
      !! `PROJECTION BASIS SET`, with GAMESS's normalization taken back out
      !!
      !! Needed to build a molecule spanning two fragments, which the
      !! inter-fragment overlaps behind exchange repulsion, charge transfer and
      !! the dispersion damping are computed from.
      !!
      !! **The printed coefficients are not the contraction coefficients.**
      !! GAMESS folds the primitive normalization in, so each is divided by
      !! `gamess_primitive_norm` at its own exponent and angular momentum --
      !! exact rather than fitted, since the writer multiplies by that same
      !! function.
      !!
      !! An `L` shell is stored as two shells sharing exponents, and the angular
      !! momentum used for the division is **per coefficient column** -- zero for the
      !! first, one for the second -- not one value for the shell.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer, parameter :: MAX_SHELLS = 512, MAX_PRIMS = 4096
      integer :: start, finish, i, atom, nprim, k, ncol, col, stat, l
      integer :: first_of(2)
      integer :: n_sh, n_pr
      integer :: sh_atom(MAX_SHELLS), sh_l(MAX_SHELLS), sh_first(MAX_SHELLS)
      integer :: sh_n(MAX_SHELLS)
      real(dp) :: expo(MAX_PRIMS), coef(MAX_PRIMS)
      real(dp) :: values(3)
      character(len=MAX_LINE) :: text
      character(len=:), allocatable :: rest
      character(len=1) :: letter

      start = 0
      do i = 1, n_lines
         if (index(adjustl(lines(i)), "PROJECTION BASIS SET") == 1) then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return
      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      atom = 0
      n_sh = 0
      n_pr = 0
      i = start
      do while (i <= finish)
         text = lines(i)
         if (len_trim(text) == 0) then
            i = i + 1
            cycle
         end if
         if (text(1:1) /= " ") then
            ! An atom header: label, centre, valence charge.
            atom = atom + 1
            i = i + 1
            cycle
         end if
         rest = adjustl(text)
         letter = rest(1:1)
         if (index("SPDFGL", letter) == 0) then
            i = i + 1
            cycle
         end if
         read (rest(2:), *, iostat=stat) nprim
         if (stat /= 0) then
            call error%set(ERROR_VALIDATION, "efp: a projection basis shell has no "// &
                           "primitive count")
            return
         end if
         ncol = 1
         if (letter == "L") ncol = 2
         if (n_sh + ncol > MAX_SHELLS .or. n_pr + ncol*nprim > MAX_PRIMS) then
            call error%set(ERROR_VALIDATION, "efp: the projection basis is larger "// &
                           "than this reader allows")
            return
         end if
         do col = 1, ncol
            n_sh = n_sh + 1
            sh_atom(n_sh) = atom
            sh_l(n_sh) = shell_l_of(letter, col)
            sh_first(n_sh) = n_pr + (col - 1)*nprim + 1
            sh_n(n_sh) = nprim
            first_of(col) = sh_first(n_sh)
         end do
         i = i + 1
         do k = 1, nprim
            if (i > finish) then
               call error%set(ERROR_VALIDATION, "efp: a projection basis shell is "// &
                              "short of its primitives")
               return
            end if
            ! index, exponent, then one coefficient per column
            rest = strip_tokens(lines(i), 1)
            do col = 1, ncol + 1
               call next_number(rest, values(col), stat)
               if (stat /= 0) then
                  call error%set(ERROR_VALIDATION, "efp: a projection basis "// &
                                 "primitive is short of its values")
                  return
               end if
            end do
            do col = 1, ncol
               l = shell_l_of(letter, col)
               expo(first_of(col) + k - 1) = values(1)
               coef(first_of(col) + k - 1) = values(1 + col) &
                                             /gamess_primitive_norm(l, values(1))
            end do
            i = i + 1
         end do
         n_pr = n_pr + ncol*nprim
      end do
      if (n_sh == 0) return

      frag%n_shells = n_sh
      allocate (frag%shell_atom(n_sh), frag%shell_l(n_sh), frag%shell_first(n_sh), &
                frag%shell_nprim(n_sh), frag%prim_expo(n_pr), frag%prim_coef(n_pr))
      frag%shell_atom = sh_atom(1:n_sh)
      frag%shell_l = sh_l(1:n_sh)
      frag%shell_first = sh_first(1:n_sh)
      frag%shell_nprim = sh_n(1:n_sh)
      frag%prim_expo = expo(1:n_pr)
      frag%prim_coef = coef(1:n_pr)
      frag%has_basis = .true.
   end subroutine read_projection_basis

   pure function shell_l_of(letter, col) result(l)
      !! The angular momentum of one coefficient column of a shell
      character(len=1), intent(in) :: letter
      integer, intent(in) :: col
      integer :: l

      select case (letter)
      case ("S")
         l = 0
      case ("P")
         l = 1
      case ("D")
         l = 2
      case ("F")
         l = 3
      case ("G")
         l = 4
      case ("L")
         ! A shared-exponent sp pair: the first column is the s, the second the p.
         l = col - 1
      case default
         l = -1
      end select
   end function shell_l_of

   subroutine read_tensor_block(lines, n_lines, name, per_record, store, n_values, &
                                present_flag, n_lmo, n_freq, error)
      !! One of the higher dispersion tensor sections, read flat
      !!
      !! `DIPOLE-QUADRUPOLE` carries 27 numbers per record and `LMOQQPOL` 81, both
      !! laid out as a block per frequency with one record per localized orbital --
      !! the same shape as the dipole dynamic section, two-token labels and all, so
      !! the same joiner serves.
      !!
      !! **Stored flat, deliberately.** Which of the 27 or 81 slots is which
      !! index tuple is not established, and reshaping now would bake in a
      !! guess.
      !!
      !! The orbital and frequency counts come from the dipole section rather
      !! than being counted again, and a section that disagrees with them is an
      !! error -- they describe the same localized orbitals at the same
      !! frequencies.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=*), intent(in) :: name
      integer, intent(in) :: per_record
      real(dp), allocatable, intent(out) :: store(:, :, :)
      integer, intent(out) :: n_values
      logical, intent(out) :: present_flag
      integer, intent(in) :: n_lmo, n_freq
      type(error_t), intent(inout) :: error

      integer :: start, finish, i, k, n_rec, lmo, freq, stat
      character(len=:), allocatable :: joined, rest
      character(len=MAX_LINE) :: text
      real(dp) :: value
      logical :: new_block

      n_values = 0
      present_flag = .false.
      if (n_lmo == 0 .or. n_freq == 0) return

      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (index(text, trim(name)) == 1) then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      n_rec = 0
      i = start
      do while (i <= finish)
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) > 0) n_rec = n_rec + 1
      end do
      if (n_rec == 0) return
      if (n_rec /= n_lmo*n_freq) then
         call error%set(ERROR_VALIDATION, "efp: "//trim(name)//" does not carry one "// &
                        "record per orbital per frequency")
         return
      end if

      n_values = per_record
      allocate (store(per_record, n_lmo, n_freq))
      store = 0.0_dp
      i = start
      freq = 0
      lmo = 0
      do while (i <= finish)
         new_block = index(lines(i), "FOR W=") > 0
         if (new_block) then
            freq = freq + 1
            lmo = 0
         end if
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         lmo = lmo + 1
         ! Two label tokens then the centroid, which the dipole section already
         ! carries, so only the tensor values are kept here.
         rest = strip_tokens(joined, 5)
         do k = 1, per_record
            call next_number(rest, value, stat)
            if (stat /= 0) then
               call error%set(ERROR_VALIDATION, "efp: a record of "//trim(name)// &
                              " is short of its values")
               return
            end if
            store(k, lmo, freq) = value
         end do
      end do
      present_flag = .true.
   end subroutine read_tensor_block

   subroutine read_static_pol(lines, n_lines, frag, error)
      !! `POLARIZABLE POINTS`: the static polarizability at each orbital centroid
      !!
      !! The same record shape as the dynamic section -- a label, a centroid,
      !! nine tensor components in GAMESS's slot order -- with one block rather
      !! than twelve, and **one** label token rather than two: `CT1` here
      !! against `CT  1` there. Taking the wrong number of tokens off the front
      !! shifts every value by one, so the count is per section rather than
      !! shared.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer :: start, finish, i, k, n_rec, at, stat
      character(len=:), allocatable :: joined, rest
      character(len=MAX_LINE) :: text
      real(dp) :: values(N_DYNAMIC_RECORD)

      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         ! Matched exactly: `DYNAMIC POLARIZABLE POINTS` and
         ! `LMOQQPOL DYNAMIC POLARIZABLE POINTS` both contain this name.
         if (trim(text) == "POLARIZABLE POINTS") then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      n_rec = 0
      i = start
      do while (i <= finish)
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) > 0) n_rec = n_rec + 1
      end do
      if (n_rec == 0) return

      frag%n_pol = n_rec
      allocate (frag%static_pol(3, 3, n_rec), frag%pol_points(3, n_rec))
      frag%static_pol = 0.0_dp
      i = start
      at = 0
      do while (i <= finish)
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         at = at + 1
         rest = strip_tokens(joined, 1)
         do k = 1, N_DYNAMIC_RECORD
            call next_number(rest, values(k), stat)
            if (stat /= 0) then
               call error%set(ERROR_VALIDATION, "efp: a polarizable point is short "// &
                              "of its twelve numbers")
               return
            end if
         end do
         frag%pol_points(:, at) = values(1:3)
         do k = 1, N_POL_SLOTS
            frag%static_pol(POL_ROW(k), POL_COL(k), at) = values(3 + k)
         end do
      end do
      frag%has_static_pol = .true.
   end subroutine read_static_pol

   subroutine read_dynamic(lines, n_lines, frag, error)
      !! `DYNAMIC POLARIZABLE POINTS`: a 3x3 tensor per localized orbital per
      !! frequency, which is what dispersion is built from
      !!
      !! The section's shape is a block per frequency, each block one record per
      !! localized orbital: a label, the orbital's centroid, then nine tensor
      !! components over continuation lines. Only the *first* record of a block
      !! carries the frequency, tagged on the end of its label line as
      !! `-- FOR W= 0.002792I A.U.`, so a new block is recognised by that tag
      !! appearing rather than by counting records.
      !!
      !! Absent is not an error: a potential written without the dynamic response
      !! is still usable for electrostatics.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer :: start, finish, i, k, n_rec, n_blocks, lmo, freq
      integer :: per_block
      character(len=:), allocatable :: joined, rest
      character(len=MAX_LINE) :: text
      real(dp) :: values(N_DYNAMIC_RECORD)
      integer :: stat
      logical :: new_block

      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (index(text, "DYNAMIC POLARIZABLE POINTS") == 1) then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      ! First pass: how many records, and how many of them before the frequency tag
      ! reappears -- that count is the number of localized orbitals.
      n_rec = 0
      n_blocks = 0
      per_block = 0
      i = start
      do while (i <= finish)
         new_block = index(lines(i), "FOR W=") > 0
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         if (new_block) n_blocks = n_blocks + 1
         n_rec = n_rec + 1
         if (n_blocks == 1) per_block = per_block + 1
      end do
      if (n_rec == 0 .or. n_blocks == 0) return
      if (per_block*n_blocks /= n_rec) then
         call error%set(ERROR_VALIDATION, "efp: the dynamic polarizability section "// &
                        "does not hold the same number of orbitals at every frequency")
         return
      end if

      frag%n_lmo = per_block
      frag%n_freq = n_blocks
      allocate (frag%dyn_pol(3, 3, per_block, n_blocks))
      allocate (frag%centroids(3, per_block), frag%frequencies(n_blocks))
      frag%frequencies = 0.0_dp

      i = start
      freq = 0
      lmo = 0
      do while (i <= finish)
         new_block = index(lines(i), "FOR W=") > 0
         if (new_block) then
            freq = freq + 1
            lmo = 0
            call frequency_of(lines(i), frag%frequencies(freq))
         end if
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         lmo = lmo + 1
         ! The label is two tokens, `CT` and its index, so three coordinates and
         ! nine tensor components follow from the third.
         rest = strip_tokens(joined, 2)
         do k = 1, N_DYNAMIC_RECORD
            call next_number(rest, values(k), stat)
            if (stat /= 0) then
               call error%set(ERROR_VALIDATION, "efp: a dynamic polarizability "// &
                              "record is short of its twelve numbers")
               return
            end if
         end do
         if (freq == 1) frag%centroids(:, lmo) = values(1:3)
         do k = 1, N_POL_SLOTS
            frag%dyn_pol(POL_ROW(k), POL_COL(k), lmo, freq) = values(3 + k)
         end do
      end do
      frag%has_dynamic = .true.
   end subroutine read_dynamic

   subroutine join_dynamic(lines, i, finish, joined)
      !! One dynamic-polarizability record: its label line plus its tensor lines
      !!
      !! **This section cannot use `join_record`.** There, a record is continued
      !! only where a line ends in `>`, and here the label line does not -- the
      !! tensor starts on the line *after* it, and the `>` marks continue within
      !! the tensor.
      !!
      !! The frequency tag is dropped: everything from `--` onwards is a comment
      !! on the label line, not data.
      character(len=*), intent(in) :: lines(:)
      integer, intent(inout) :: i
      integer, intent(in) :: finish
      character(len=:), allocatable, intent(out) :: joined

      character(len=MAX_LINE) :: text
      integer :: cut
      logical :: continued

      joined = ""
      if (i > finish) return
      text = lines(i)
      i = i + 1
      cut = index(text, "--")
      if (cut > 0) then
         joined = trim(text(1:cut - 1))
      else
         joined = trim(text)
      end if

      do
         if (i > finish) exit
         text = lines(i)
         ! A new label line means this record is done; do not consume it.
         if (index(adjustl(text), "CT") == 1) exit
         i = i + 1
         cut = index(text, ">")
         continued = cut > 0
         if (continued) then
            joined = joined//" "//trim(text(1:cut - 1))
         else
            joined = joined//" "//trim(text)
         end if
         if (.not. continued) exit
      end do
      joined = adjustl(joined)
   end subroutine join_dynamic

   function strip_tokens(text, n) result(rest)
      !! `text` with its first `n` whitespace-separated tokens removed
      character(len=*), intent(in) :: text
      integer, intent(in) :: n
      character(len=:), allocatable :: rest

      integer :: k, space

      rest = adjustl(text)
      do k = 1, n
         space = index(trim(rest), " ")
         if (space == 0) then
            rest = ""
            return
         end if
         rest = adjustl(rest(space:))
      end do
   end function strip_tokens

   subroutine frequency_of(line, w)
      !! The imaginary frequency out of `-- FOR W= 0.002792I A.U.`
      character(len=*), intent(in) :: line
      real(dp), intent(out) :: w

      integer :: at, stop_at, stat

      w = 0.0_dp
      at = index(line, "FOR W=")
      if (at == 0) return
      stop_at = index(line(at:), "I")
      if (stop_at == 0) return
      read (line(at + 6:at + stop_at - 2), *, iostat=stat) w
      if (stat /= 0) w = 0.0_dp
   end subroutine frequency_of

   subroutine one_section(lines, n_lines, name, n_values, n_expected, values, error)
      !! A section that must be present and must carry one record per point
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=*), intent(in) :: name
      integer, intent(in) :: n_values, n_expected
      real(dp), allocatable, intent(out) :: values(:, :)
      type(error_t), intent(inout) :: error

      character(len=MAX_LINE), allocatable :: labels(:)
      integer :: n_rec

      call section_records(lines, n_lines, name, n_values, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec == 0) then
         call error%set(ERROR_VALIDATION, "efp: no "//name//" section")
         return
      end if
      call expect_points(n_rec, n_expected, name, error)
      if (allocated(labels)) deallocate (labels)
   end subroutine one_section

   subroutine expect_points(got, want, name, error)
      integer, intent(in) :: got, want
      character(len=*), intent(in) :: name
      type(error_t), intent(inout) :: error

      character(len=32) :: a, b

      if (got /= want) then
         write (a, "(I0)") got
         write (b, "(I0)") want
         call error%set(ERROR_VALIDATION, "efp: "//name//" carries "//trim(a)// &
                        " records but COORDINATES has "//trim(b)//" points")
      end if
   end subroutine expect_points

   subroutine section_records(lines, n_lines, name, n_values, labels, values, n_rec, error)
      !! Every record of one named section, continuations joined
      !!
      !! Returns `n_rec = 0` rather than an error when the section is absent, so an
      !! optional section costs the caller one branch.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=*), intent(in) :: name
      integer, intent(in) :: n_values
      character(len=MAX_LINE), allocatable, intent(out) :: labels(:)
      real(dp), allocatable, intent(out) :: values(:, :)
      integer, intent(out) :: n_rec
      type(error_t), intent(inout) :: error

      integer :: i, start, finish, k, count_rec
      character(len=:), allocatable :: joined
      character(len=MAX_LINE) :: text

      n_rec = 0
      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (trim(text) == trim(name) .or. index(trim(text), trim(name)//" ") == 1) then
            ! A header, not a record: records carry a label in column one after
            ! adjustment and this matched the section name outright.
            start = i + 1
            exit
         end if
      end do
      if (start == 0) then
         allocate (labels(0), values(n_values, 0))
         return
      end if

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      ! First pass: count records, a record being a line that does not continue the
      ! one before it.
      count_rec = 0
      i = start
      do while (i <= finish)
         call join_record(lines, i, finish, joined)
         if (len_trim(joined) > 0) count_rec = count_rec + 1
      end do
      allocate (labels(count_rec), values(n_values, count_rec))
      values = 0.0_dp

      i = start
      k = 0
      do while (i <= finish)
         call join_record(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         k = k + 1
         call split_record(joined, n_values, labels(k), values(:, k), name, error)
         if (error%has_error()) return
      end do
      n_rec = k
   end subroutine section_records

   subroutine join_record(lines, i, finish, joined)
      !! One logical record starting at `i`, following `>` continuations
      !!
      !! `i` is advanced past everything consumed, so the caller's loop needs no
      !! separate bookkeeping for how many physical lines a record took.
      character(len=*), intent(in) :: lines(:)
      integer, intent(inout) :: i
      integer, intent(in) :: finish
      character(len=:), allocatable, intent(out) :: joined

      character(len=MAX_LINE) :: text
      integer :: cut

      joined = ""
      do
         if (i > finish) exit
         text = lines(i)
         i = i + 1
         cut = index(text, ">")
         if (cut > 0) then
            joined = joined//" "//trim(text(1:cut - 1))
         else
            joined = joined//" "//trim(text)
            exit
         end if
      end do
      joined = adjustl(joined)
   end subroutine join_record

   subroutine split_record(record, n_values, label, values, name, error)
      !! A record's label and its numbers
      character(len=*), intent(in) :: record
      integer, intent(in) :: n_values
      character(len=*), intent(out) :: label
      real(dp), intent(out) :: values(:)
      character(len=*), intent(in) :: name
      type(error_t), intent(inout) :: error

      integer :: space, stat, k
      character(len=:), allocatable :: rest
      character(len=32) :: a

      space = index(trim(record), " ")
      if (space == 0) then
         call error%set(ERROR_VALIDATION, "efp: "//name//": record '"// &
                        trim(record)//"' carries no values")
         return
      end if
      label = record(1:space - 1)
      rest = adjustl(record(space:))
      do k = 1, n_values
         call next_number(rest, values(k), stat)
         if (stat /= 0) then
            write (a, "(I0)") n_values
            call error%set(ERROR_VALIDATION, "efp: "//name//": record '"// &
                           trim(label)//"' needs "//trim(a)//" values")
            return
         end if
      end do
   end subroutine split_record

   subroutine next_number(text, value, stat)
      !! Read one number off the front of `text` and consume it
      character(len=:), allocatable, intent(inout) :: text
      real(dp), intent(out) :: value
      integer, intent(out) :: stat

      integer :: space

      value = 0.0_dp
      text = adjustl(text)
      if (len_trim(text) == 0) then
         stat = 1
         return
      end if
      space = index(trim(text), " ")
      if (space == 0) then
         read (text, *, iostat=stat) value
         text = ""
      else
         read (text(1:space - 1), *, iostat=stat) value
         text = adjustl(text(space:))
      end if
   end subroutine next_number

   subroutine fragment_name(lines, n_lines, name)
      !! `$FRAGNAME` without its `$`, or a placeholder if the file omits it
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=:), allocatable, intent(out) :: name

      integer :: i
      character(len=MAX_LINE) :: text

      name = "FRAGMENT"
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (text(1:1) == "$") then
            name = trim(text(2:))
            return
         end if
      end do
   end subroutine fragment_name

   subroutine multiplicity(lines, n_lines, mult)
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(out) :: mult

      integer :: i, stat
      character(len=MAX_LINE) :: text

      mult = 1
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (index(text, "MULTIPLICITY") == 1) then
            read (text(13:), *, iostat=stat) mult
            if (stat /= 0) mult = 1
            return
         end if
      end do
   end subroutine multiplicity

   subroutine slurp(path, lines, n_lines, error)
      !! The whole file as lines
      !!
      !! Read twice rather than grown, so the count is known before allocating.
      character(len=*), intent(in) :: path
      character(len=MAX_LINE), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: n_lines
      type(error_t), intent(inout) :: error

      integer :: unit, stat, i
      character(len=MAX_LINE) :: text

      n_lines = 0
      open (newunit=unit, file=path, status="old", action="read", iostat=stat)
      if (stat /= 0) then
         call error%set(ERROR_IO, "efp: cannot open '"//trim(path)//"'")
         return
      end if
      do
         read (unit, "(A)", iostat=stat) text
         if (stat /= 0) exit
         n_lines = n_lines + 1
      end do
      rewind (unit)
      allocate (lines(max(n_lines, 1)))
      lines = ""
      do i = 1, n_lines
         read (unit, "(A)", iostat=stat) lines(i)
         if (stat /= 0) exit
      end do
      close (unit)
   end subroutine slurp

end module mqc_czt_efp_read
