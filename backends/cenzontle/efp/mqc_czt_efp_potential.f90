!! A complete effective fragment potential, computed and written here
module mqc_czt_efp_potential
   !! `RUNTYP=MAKEFP`: a geometry and a basis name in, a `.efp` file out.
   !!
   !! The assembly: the SCF, the localization, the multipoles, the static and
   !! dynamic polarizabilities, the projection data and both screening fits,
   !! written as one file. Each block has a `validation/check_*` program of its
   !! own that compares it against GAMESS's printed numbers.
   !!
   !! **All seventeen sections GAMESS's reader recognises are written**:
   !! electrostatics to octupole with charge-penetration screening, polarization,
   !! exchange repulsion, dispersion through `E6`, `E7` and `E8`, and charge
   !! transfer. Seventeen and not eighteen because `CTFOK` is a *subsection* of
   !! `CTVEC` rather than a section -- GAMESS accepts it only directly behind one
   !! and aborts on a standalone one, so the two share a `STOP`.
   !!
   !! **The formats target GAMESS's reader, not byte-identity with its writer.**
   !! What is checked is that GAMESS accepts the file and agrees with the energies
   !! it computes from it; `tools/efp_validation/dimer_energy.py` asks that.
   !!
   !! **`LMOQQPOL` is the one block validated by its energy rather than its
   !! values.** Its per-orbital tensors differ from GAMESS's written ones, in the
   !! part antisymmetric under exchanging the two index pairs, which cancels on
   !! summing: the response summed over orbitals reproduces GAMESS's own
   !! *molecular* quadrupole-quadrupole polarizability more closely than its own
   !! written per-orbital block does.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_mass
   use mqc_physical_constants, only: PI
   use mqc_calculation_defaults, only: DEFAULT_VDW_SCALE
   use mqc_cgto, only: molecular_basis_type
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, shell_dim
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_atomic_guess, only: build_restricted_guess, guess_display_name
   use mqc_czt_localize, only: boys_localize
   use mqc_czt_dma, only: dma_result_t, distributed_multipoles
   use mqc_czt_cphf, only: response_hessian_t, distributed_polarizability, &
                           distributed_dynamic_polarizability, &
                           distributed_dynamic_cross, &
                           casimir_polder_frequencies, N_CASIMIR_POLDER
   use mqc_czt_multipole, only: multipole_matrices
   use mqc_czt_screening, only: fit_screening, screening_target_t, &
                                SCREEN_EXPONENTIAL, SCREEN_GAUSSIAN
   use pic_timer, only: timer_type
   use libcint_fortran, only: LIBCINT_ANG_OF
   use pic_logger, only: logger => global_logger
   use mqc_scf_types, only: scf_numerics_t, print_scf_config
   use mqc_diis, only: parse_accelerator_name, ACCEL_DIIS
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: efp_potential_t
   public :: make_efp_potential
   public :: write_efp_potential
   ! Exposed so a reader can invert it: the printed contraction coefficients
   ! carry this factor, and recovering the raw ones is an exact division by it.
   public :: gamess_primitive_norm
   public :: from_gamess_ao_order
   ! Its inverse, for anything that changes an orbital and must hand the result
   ! back to a fragment still storing GAMESS's order.
   public :: to_gamess_ao_order
   public :: frozen_core

   integer, parameter :: MAX_LINE = 160
      !! Longest line any section emits, with room to spare.

   real(dp), parameter :: MAKEFP_DENSITY_TOL = 1.0e-8_dp
      !! What a fragment potential is fitted at, and this module's default for
      !! both the density threshold and the commutator gate derived from it.
      !! Named because the warning below compares against the same number.

   integer, parameter :: N_CART_PAIR = 9
      !! Flattened extents of the Cartesian tensors the polarizability blocks
      !! carry: a pair of directions, a triple, and a quadruple. Nine is also the
      !! number of slots GAMESS writes a polarizability in.
   integer, parameter :: N_CART_TRIPLE = 27
   integer, parameter :: N_CART_QUAD = 81

   integer, parameter :: N_CART_D = 6
      !! Components in a Cartesian d and f shell.
   integer, parameter :: N_CART_F = 10

   integer, parameter :: POL_ROW(N_CART_PAIR) = [1, 2, 3, 2, 3, 3, 1, 1, 2]
      !! Row and column of each of GAMESS's nine polarizability slots. **The
      !! off-diagonal triples are the transpose of what its labels suggest**,
      !! measured in `validation/check_distributed_polarizability.py` -- the one
      !! convention here that a symmetric test tensor would not have caught.
      !!
      !! Also the transpose of what `efinp.src:7552-7561` writes, which is not a
      !! conflict: GAMESS indexes the tensor `(field, dipole)` where this code
      !! indexes it `(dipole, field)`, so both put the same number in slot 4. See
      !! the note on `POL_ROW` in `mqc_czt_efp_read`, which carries the measurement.
   integer, parameter :: POL_COL(N_CART_PAIR) = [1, 2, 3, 1, 1, 2, 2, 3, 3]

   integer, parameter :: QXX = 1, QXY = 2, QXZ = 3, QYX = 4, QYY = 5
   integer, parameter :: QYZ = 6, QZX = 7, QZY = 8, QZZ = 9
      !! libcint's full-Cartesian quadrupole slots, which run xx,xy,xz,yx,...,zz.

   integer, parameter :: D_FROM_LIBCINT(N_CART_D) = [1, 4, 6, 2, 3, 5]
   real(dp), parameter :: D_NORMALIZATION = 1.585330892_dp
      !! libcint's index for each of GAMESS's six Cartesian d slots, and the
      !! normalization between the two codes' d functions. Both established in
      !! `validation/check_projection.py` against GAMESS's own coefficients.

   integer, parameter :: F_FROM_LIBCINT(N_CART_F) = [1, 7, 10, 2, 3, 4, 8, 6, 9, 5]
   real(dp), parameter :: F_NORMALIZATION = 1.339849174_dp
      !! The same for the ten Cartesian f slots, read off GAMESS's own
      !! coefficients in `validation/check_projection` and solved for there.
      !!
      !! **The molecule has to be in a frame with no zero coordinate for this to
      !! be solvable at all**: planar water puts an exact zero in every function
      !! with an odd power of y, and a slot that is zero on both sides admits any
      !! scale factor, so a plausible-looking map can send several GAMESS slots
      !! to one of ours.

   integer, parameter :: F_CLASS(N_CART_F) = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3]
      !! Which of the three normalization classes each GAMESS f slot belongs to.
      !! The measured scales come out as `F_NORMALIZATION` divided by one of 1,
      !! sqrt(5) and sqrt(15), exact to eight figures.

   type :: efp_potential_t
      !! Every parameter a `.efp` carries that we can compute
      character(len=:), allocatable :: name        !! `$FRAGNAME`, without the `$`
      character(len=:), allocatable :: basis_name
      integer :: n_points = 0     !! Expansion points: atoms then bond midpoints
      integer :: n_atoms = 0
      integer :: nao = 0
      integer :: n_occ = 0        !! Including the core, which `CTFOK` needs
      integer :: n_lmo = 0        !! Valence localized orbitals
      integer :: multiplicity = 1
      real(dp) :: scf_energy = 0.0_dp
         !! The RHF total energy of the monomer this potential was made from,
         !! nuclear repulsion included. **EFMO's `E_I^0`**: the fragment sum of
         !! eq 6 is these numbers, so the SCF behind a potential is not run a
         !! second time to get them.
      logical :: quadrupole_blocks = .true.
         !! Whether `dipquad` and `quadquad` are computed and written
      real(dp) :: vdwscl = DEFAULT_VDW_SCALE
         !! The screening grid's van der Waals scale. `fit_screening` is handed
         !! this rather than reading the default for itself, so the grid and the
         !! written header agree.
      character(len=8), allocatable :: labels(:)      !! `A01O`, `BO21`, ...
      real(dp), allocatable :: points(:, :)           !! (3, n_points), Bohr
      real(dp), allocatable :: mass(:)                !! amu, zero at a midpoint
      real(dp), allocatable :: charge(:)              !! Z, zero at a midpoint
      real(dp), allocatable :: q_elec(:), q_nuc(:)    !! MONOPOLES
      real(dp), allocatable :: dipole(:, :)           !! (3, n_points)
      real(dp), allocatable :: quadrupole(:, :)       !! (6, n_points)
      real(dp), allocatable :: octopole(:, :)         !! (10, n_points)
      real(dp), allocatable :: centroids(:, :)        !! (3, n_lmo)
      real(dp), allocatable :: static_pol(:, :, :)    !! (3, 3, n_lmo)
      real(dp), allocatable :: dynamic_pol(:, :, :, :)  !! (3, 3, n_lmo, n_freq)
      real(dp), allocatable :: dipquad_pre(:, :, :, :, :)
         !! The dipole-quadrupole tensor *before* the translation, which is what
         !! `QQSHIFT` takes as an input -- so it has to be kept, not just the
         !! shifted form the file carries.
      real(dp), allocatable :: quadquad(:, :, :, :, :, :)
         !! `(3, 3, 3, 3, n_lmo, n_freq)`, after the write-time translation.
      real(dp), allocatable :: dipquad(:, :, :, :, :)
         !! `(3, 3, 3, n_lmo, n_freq)` as `A'(a,b,c)`, **after** the write-time
         !! translation to each centroid, which is the form the file carries.
      real(dp), allocatable :: frequencies(:)         !! Imaginary, a.u.
      real(dp), allocatable :: fock_lmo(:, :)         !! (n_lmo, n_lmo)
      real(dp), allocatable :: orbitals(:, :)         !! LMOs in GAMESS's AO order
      real(dp), allocatable :: canonical(:, :)
         !! All the canonical MOs, in GAMESS's AO order. `CTVEC` in its
         !! canonical-orbital form is exactly this matrix.
      real(dp), allocatable :: eps_occ(:)             !! CTFOK
      real(dp), allocatable :: screen2(:)             !! Exponential alpha per point
      real(dp), allocatable :: screen(:)              !! Gaussian alpha per point
      character(len=MAX_LINE), allocatable :: basis_lines(:)
   contains
      procedure :: destroy => potential_destroy
   end type efp_potential_t

contains

   subroutine potential_destroy(self)
      class(efp_potential_t), intent(inout) :: self

      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%basis_name)) deallocate (self%basis_name)
      if (allocated(self%labels)) deallocate (self%labels)
      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%mass)) deallocate (self%mass)
      if (allocated(self%charge)) deallocate (self%charge)
      if (allocated(self%q_elec)) deallocate (self%q_elec)
      if (allocated(self%q_nuc)) deallocate (self%q_nuc)
      if (allocated(self%dipole)) deallocate (self%dipole)
      if (allocated(self%quadrupole)) deallocate (self%quadrupole)
      if (allocated(self%octopole)) deallocate (self%octopole)
      if (allocated(self%centroids)) deallocate (self%centroids)
      if (allocated(self%static_pol)) deallocate (self%static_pol)
      if (allocated(self%dynamic_pol)) deallocate (self%dynamic_pol)
      if (allocated(self%dipquad)) deallocate (self%dipquad)
      if (allocated(self%quadquad)) deallocate (self%quadquad)
      if (allocated(self%dipquad_pre)) deallocate (self%dipquad_pre)
      if (allocated(self%frequencies)) deallocate (self%frequencies)
      if (allocated(self%fock_lmo)) deallocate (self%fock_lmo)
      if (allocated(self%orbitals)) deallocate (self%orbitals)
      if (allocated(self%canonical)) deallocate (self%canonical)
      if (allocated(self%eps_occ)) deallocate (self%eps_occ)
      if (allocated(self%screen2)) deallocate (self%screen2)
      if (allocated(self%screen)) deallocate (self%screen)
      if (allocated(self%basis_lines)) deallocate (self%basis_lines)
      self%scf_energy = 0.0_dp
      self%n_points = 0
      self%n_atoms = 0
      self%nao = 0
      self%n_occ = 0
      self%n_lmo = 0
   end subroutine potential_destroy

   subroutine make_efp_potential(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, pot, error, charge, n_core, &
                                 vdwscl, verbose, aux_basis, guess, &
                                 energy_tol, density_tol, grad_tol_in, &
                                 scf_in, max_iter_in, dynamic_tol, &
                                 dynamic_maxiter, response, allow_crap_response, &
                                 response_batch, quadrupole_blocks, scf_out)
      !! The whole pipeline: SCF, localization, and every parameter block
      !!
      !! The order is forced by what depends on what: the SCF gives the density
      !! the multipoles are taken from and the orbitals everything else needs,
      !! localization gives the centres the polarizabilities and the exchange
      !! repulsion data sit on, and the screening is fitted last because what it
      !! fits is the error the *damped multipole* potential makes.
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)    !! (3, natm), Bohr
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: name         !! Fragment name, e.g. `WATER`
      type(efp_potential_t), intent(out) :: pot
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge
         !! Net charge. Ignoring it makes the electron count wrong, and odd for a
         !! cation, so a closed-shell reference is then refused outright.
      integer, intent(in), optional :: n_core
         !! Orbitals excluded from the localized set. Default is the standard
         !! frozen core, which is what MAKEFP uses: its polarizable points and
         !! its exchange-repulsion orbitals are valence only.
      logical, intent(in), optional :: quadrupole_blocks
         !! `keywords.efp.dispersion`: the dipole-quadrupole and
         !! quadrupole-quadrupole dynamic blocks too, or the dipole-dipole one
         !! alone. Absent is all three.
      real(dp), intent(in), optional :: vdwscl
         !! `keywords.efp.vdw_scale`. Where the innermost layer of the screening
         !! grid sits, as a fraction of a van der Waals radius, and what the
         !! written `SCREEN`/`SCREEN2` headers report having used. One number for
         !! both.
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: guess
         !! Initial-guess name from the deck (`keywords.guess.type`). Default is
         !! "auto", which resolves to SAD on this backend -- the same guess the
         !! Energy driver uses.
      character(len=*), intent(in), optional :: aux_basis
         !! Fit the dynamic response rather than building its Hessian exactly.
         !! That build is `n_ov` Fock builds and is most of what a potential
         !! costs; fitted, it is two matrix products. **The approximation is real
         !! and is measured by `validation/check_df_hessian`**, so it is asked for
         !! rather than inferred.
         !!
         !! The auxiliary basis must match the orbital basis in angular form,
         !! which `build_df_mo_block` checks: libcint builds all three centres of
         !! a fitting integral in one form.
      real(dp), intent(in), optional :: energy_tol, density_tol
         !! SCF convergence thresholds. Present only when a deck named
         !! `keywords.scf.tolerance` / `keywords.scf.density_tolerance`; absent,
         !! the tight defaults at the SCF call stand.
      type(scf_numerics_t), intent(in), optional :: scf_in
         !! The deck's `keywords.scf`, for every setting this routine has no
         !! opinion about: the level shift, the accelerator, the DIIS subspace,
         !! the linear-dependence threshold and incremental Fock building. Absent
         !! leaves the SCF's own defaults.
      integer, intent(in), optional :: max_iter_in
         !! `keywords.scf.maxiter`, present only when the deck named it. The 200
         !! below is this routine's own, deliberately larger than the shared
         !! default of 100 because a fragment potential is converged tightly.
      real(dp), intent(in), optional :: grad_tol_in
         !! `keywords.scf.gradient_tolerance`, present only when the deck named
         !! it. Wins outright over the rules at the SCF call below.
      real(dp), intent(in), optional :: dynamic_tol
         !! `keywords.efp.dynamic_tolerance`. What the frequency-dependent
         !! response solve converges its residual to. **Alone among the EFP keys
         !! this one moves the numbers in the file rather than the route to
         !! them**: the dynamic polarizabilities, and every dispersion energy
         !! taken from them, are only as converged as this says.
      logical, intent(in), optional :: allow_crap_response
         !! `keywords.efp.allow_crap_response`. Accept whatever the response
         !! solve reached. The potential is wrong; see `efp_config_t`.
      integer, intent(in), optional :: response_batch
         !! `keywords.efp.response_batch`. Tuning only; the answer is unchanged.
      integer, intent(in), optional :: dynamic_maxiter
         !! `keywords.efp.dynamic_maxiter`. Iterations that solve gets before it
         !! reports a failure to converge.
      integer, intent(in), optional :: response
         !! `keywords.efp.response`, as one of the `EFP_RESPONSE_*` codes: build the
         !! response operator, never build it, or let the size rule decide. Passed
         !! through to `dynamic_polarizability`, which is where the choice is made.

      type(rhf_result_t), intent(out), optional :: scf_out
         !! The converged monomer SCF, handed back rather than discarded.
         !!
         !! **This is what EFMO's `E_I^0` is.** A fragment's in-vacuo energy and
         !! its potential come from the same reference determinant, so an EFMO
         !! run that asked for both separately would run every monomer's SCF
         !! twice; `pot%scf_energy` carries the total alone and this carries the
         !! orbitals a correlated `E_I^0` would continue from.

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(dma_result_t) :: dma
      type(response_hessian_t) :: shared_hessian
      real(dp), allocatable :: u_static(:, :, :)
      real(dp), allocatable :: b_ao(:, :)
      type(screening_target_t) :: screen_target
      real(dp), allocatable :: loc(:, :), ovl(:, :), sc(:, :), w(:, :), scaled(:, :)
      real(dp), allocatable :: alpha(:)
      real(dp) :: rms_exp, rms_gauss
      integer :: natm, core, i, j, k, n_valence, n_electrons
      integer :: guess_kind
      real(dp) :: e_tol, d_tol, g_tol
      character(len=16) :: tol_text, ref_text
      integer :: n_iter, scf_diis, accel_kind
      real(dp) :: shift, lindep
      logical :: incr, accel_ok
      type(scf_numerics_t) :: echo
      real(dp), allocatable :: guess_total(:, :)
      character(len=:), allocatable :: guess_name
      type(timer_type) :: stage
      logical :: talk

      character(len=MAX_LINE_LENGTH) :: line

      talk = .false.
      if (present(verbose)) talk = verbose
      if (talk) call stage%start()
      natm = size(atomic_numbers)
      n_electrons = sum(atomic_numbers)
      if (present(charge)) n_electrons = n_electrons - charge
      pot%name = trim(name)
      pot%basis_name = trim(basis_name)
      pot%n_atoms = natm
      if (present(vdwscl)) pot%vdwscl = vdwscl
      if (present(quadrupole_blocks)) pot%quadrupole_blocks = quadrupole_blocks

      if (size(coordinates, 1) /= 3 .or. size(coordinates, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "makefp: coordinates must be (3, natm)")
         return
      end if

      ! **Cartesian, because a `.efp` is a GAMESS file.** The Basis Set Exchange
      ! is not consistent about Pople sets -- 6-31G* declares its d Cartesian
      ! while 6-311++G(3df,3pd) declares its d and f spherical -- and
      ! `build_czt_molecule` follows the declaration unless told otherwise.
      ! Here there is nothing to decide: the potential is read back by GAMESS,
      ! whose ISPHER default is -1, and the AO ordering map in
      ! `from_gamess_ao_order` is a Cartesian map. It is also what lets an f
      ! basis work at all.
      !
      ! **The accelerator name is validated before the molecule, the basis and
      ! the guess.** With the default "auto" guess this routine builds a SAD
      ! density from free-atom SCFs first, so validating afterwards meant paying
      ! for every one of them before saying the word was wrong.
      if (present(scf_in)) then
         call parse_accelerator_name(scf_in%accelerator, accel_kind, accel_ok)
         if (.not. accel_ok) then
            call error%set(ERROR_VALIDATION, "keywords.scf.accelerator '"// &
                           trim(scf_in%accelerator)//"' is not one of diis, adiis, ediis")
            return
         end if
      end if

      call build_czt_molecule(atomic_numbers, element_symbols, coordinates, &
                              basis_name, mol, error, force_cartesian=.true.)
      if (error%has_error()) return
      pot%nao = mol%nao

      if (talk) then
         write (line, "(A,A,A,I0,A)") "  basis ", trim(basis_name), ", ", mol%nao, " functions"
         call logger%info(trim(line))
      end if

      if (present(aux_basis)) then
         ! Read the fitting set in whatever angular form the orbital basis is in.
         ! libcint builds all three centres of a fitting integral in one form, so
         ! the two have to agree -- and the writer needs the orbital basis
         ! Cartesian, while every fitting set on hand is declared spherical.
         !
         ! Legitimate because an auxiliary basis is a fitting space, not a
         ! wavefunction. It does make the space redundant -- the Cartesian d
         ! shells carry an s contaminant that duplicates the aux s functions --
         ! so the metric is more nearly singular, and what that costs is measured
         ! in `validation/check_df_hessian`.
         call build_czt_molecule(atomic_numbers, element_symbols, coordinates, &
                                 aux_basis, aux, error, &
                                 force_cartesian=mol%cartesian)
         if (error%has_error()) then
            call mol%destroy()
            return
         end if
         if (talk) then
            write (line, "(A,A,A,I0,A)") "  density fitting: ", trim(aux_basis), ", ", aux%nao, &
               " functions (the SCF and the dynamic response)"
            call logger%info(trim(line))
         end if
      else if (talk) then
         call logger%info("  density fitting: off (the SCF and the dynamic response are exact)")
      end if

      ! Checked here, not where the ordering map needs it. The map runs after the
      ! SCF, the localization and every response solve, so a basis this cannot emit
      ! would otherwise be refused several minutes in.
      call check_angular_form(mol, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! Honour the deck's guess, the same as the Energy path. `guess_total` is
      ! left unallocated for core/gwh and the SCF reads it only for the atomic
      ! guesses, so passing it unconditionally is safe.
      guess_name = "auto"
      if (present(guess)) guess_name = guess
      call build_restricted_guess(mol, guess_name, guess_kind, guess_total, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) call logger%info("  initial guess: "//guess_display_name(guess_kind))

      ! 1e-10 energy / 1e-8 density *when the deck does not say otherwise*: the
      ! density is what the multipoles and the response are taken from, and 1e-8
      ! is tight enough that the potential does not move in any digit it reports.
      ! The energy keeps the 100:1 ratio the rest of the code uses, so the
      ! density is the binding criterion. Name the key and it is honoured; leave
      ! it alone and the tight pair stands.
      !
      ! `aux` present means density-fit the SCF too, not just the response.
      e_tol = 1.0e-10_dp
      d_tol = MAKEFP_DENSITY_TOL
      if (present(energy_tol)) e_tol = energy_tol
      if (present(density_tol)) d_tol = density_tol
      ! **The commutator threshold is stated, not derived.** Everything a
      ! fragment potential carries is fitted to the SCF *density* -- the
      ! multipoles, the polarizabilities, the screening -- and the density's
      ! error goes as the commutator where the energy's goes as its square, so
      ! `sqrt(e_tol)` would be 1e-5 and the multipoles drift off their GAMESS
      ! references. Tightening `e_tol` instead would need 1e-16, below what a
      ! molecular energy resolves.
      !
      ! A convergence measure the user names is the measure that decides:
      !
      !   * `gradient_tolerance` named -- it is the commutator threshold.
      !   * `density_tolerance` named -- the commutator follows it.
      !   * only `tolerance` named -- derive it as every other SCF here does,
      !     `sqrt(e_tol)`, so that loosening the energy loosens the run.
      !   * nothing named -- the tight pair stands.
      !
      ! The third case is not overridden. Nor is it warned about here: the
      ! bound is where the SCF is *allowed* to stop, and with DIIS the energy
      ! criterion binds first and the commutator lands orders below it. What
      ! the density was actually left at is checked after the SCF, below.
      g_tol = d_tol
      if (present(grad_tol_in)) then
         g_tol = grad_tol_in
      else if (.not. present(density_tol) .and. present(energy_tol)) then
         g_tol = sqrt(e_tol)
      end if
      ! Everything the deck said about how an SCF runs, forwarded.
      n_iter = 200
      if (present(max_iter_in)) n_iter = max_iter_in
      scf_diis = 8
      accel_kind = ACCEL_DIIS
      shift = 0.0_dp
      lindep = 0.0_dp
      incr = .true.
      if (present(scf_in)) then
         scf_diis = scf_in%diis_size
         if (.not. scf_in%use_diis) scf_diis = 0
         shift = scf_in%level_shift
         lindep = scf_in%linear_dependence
         incr = scf_in%incremental_fock
         ! Parsed again rather than re-derived: the spelling was already
         ! validated at the top of this routine, so `accel_ok` cannot be false
         ! here.
         call parse_accelerator_name(scf_in%accelerator, accel_kind, accel_ok)
      end if
      ! Echoed before the SCF runs, so a deck that set something and saw no
      ! effect can be checked against what arrived.
      !
      ! Assembled whether or not it is printed, because this is what the SCF is
      ! *given* and not merely what gets reported -- a configuration echoed but
      ! not passed would be the failure this echo exists to catch.
      echo = scf_numerics_t()
      if (present(scf_in)) echo = scf_in
      echo%max_iter = n_iter
      echo%energy_tol = e_tol
      echo%density_tol = d_tol
      echo%grad_tol = g_tol
      echo%diis_size = scf_diis
      echo%level_shift = shift
      echo%linear_dependence = lindep
      echo%incremental_fock = incr
      ! The resolved guess, not the deck spelling: echoing "auto" would name the
      ! request rather than the run.
      echo%guess = guess_display_name(guess_kind)
      if (talk) call print_scf_config(echo, "MAKEFP SCF")
      if (present(aux_basis)) then
         call run_czt_rhf(mol, n_electrons, n_iter, e_tol, d_tol, &
                          talk, scf, error, guess=guess_kind, guess_density=guess_total, &
                          aux=aux, grad_tol=g_tol, diis_vectors=scf_diis, &
                          level_shift=shift, linear_dependence=lindep, &
                          accelerator=accel_kind, incremental_fock=incr, &
                          scf=echo, b_ao_out=b_ao)
      else
         call run_czt_rhf(mol, n_electrons, n_iter, e_tol, d_tol, &
                          talk, scf, error, guess=guess_kind, guess_density=guess_total, &
                          grad_tol=g_tol, diis_vectors=scf_diis, &
                          level_shift=shift, linear_dependence=lindep, &
                          accelerator=accel_kind, incremental_fock=incr, &
                          scf=echo)
      end if
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "makefp: the SCF did not converge")
         call mol%destroy()
         return
      end if
      ! Judged on the commutator the SCF reached, not the one it was allowed
      ! to stop at: the multipoles, the polarizabilities and the screening
      ! are all fitted to this density, and its error goes as the commutator.
      if (scf%commutator > MAKEFP_DENSITY_TOL) then
         write (tol_text, "(es9.2)") scf%commutator
         write (ref_text, "(es9.2)") MAKEFP_DENSITY_TOL
         call logger%warning("  MAKEFP: the SCF stopped at a commutator of "// &
                             trim(adjustl(tol_text))//", above the "// &
                             trim(adjustl(ref_text))//" a "// &
                             "fragment potential is fitted at, so the multipoles "// &
                             "and polarizabilities may differ from a reference in "// &
                             "their last digits. Tighten keywords.scf.tolerance, "// &
                             "or name gradient_tolerance, to converge it further.")
      end if
      pot%n_occ = scf%n_occupied
      pot%scf_energy = scf%energy
      if (present(scf_out)) scf_out = scf
      if (talk) call report(stage, "SCF", talk)
      if (talk) then
         write (line, "(A,F18.10)") "  RHF energy ", scf%energy
         call logger%info(trim(line))
      end if

      core = frozen_core(atomic_numbers)
      if (present(n_core)) core = n_core
      n_valence = pot%n_occ - core
      if (n_valence < 1) then
         call error%set(ERROR_VALIDATION, "makefp: no valence orbitals to localize")
         call mol%destroy()
         return
      end if
      pot%n_lmo = n_valence

      allocate (pot%eps_occ(pot%n_occ))
      pot%eps_occ = scf%orbital_energies(1:pot%n_occ)

      ! --- localization: the centres everything below is expressed on ------------
      call boys_localize(mol, scf%orbitals(:, core + 1:pot%n_occ), n_valence, &
                         loc, pot%centroids, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) then
         write (line, "(A,I0,A,I0,A)") "  localized ", n_valence, " of ", pot%n_occ, " occupied orbitals"
         call logger%info(trim(line))
      end if
      if (talk) call report(stage, "localization", talk)

      ! --- electrostatics -------------------------------------------------------
      call distributed_multipoles(mol, scf%density, atomic_numbers, dma, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      pot%n_points = size(dma%labels)
      allocate (pot%labels(pot%n_points), pot%points(3, pot%n_points))
      pot%labels = dma%labels
      pot%points = dma%points
      allocate (pot%q_elec(pot%n_points), pot%q_nuc(pot%n_points))
      pot%q_elec = dma%electronic
      pot%q_nuc = dma%nuclear
      allocate (pot%dipole(3, pot%n_points), pot%quadrupole(6, pot%n_points))
      allocate (pot%octopole(10, pot%n_points))
      pot%dipole = dma%dipole
      pot%quadrupole = dma%quadrupole
      pot%octopole = dma%octopole

      ! Mass and charge sit on atoms only; a bond midpoint carries neither, which
      ! is how GAMESS's reader tells the two kinds of point apart.
      allocate (pot%mass(pot%n_points), pot%charge(pot%n_points))
      pot%mass = 0.0_dp
      pot%charge = 0.0_dp
      do i = 1, natm
         pot%mass(i) = element_mass(atomic_numbers(i))
         pot%charge(i) = real(atomic_numbers(i), dp)
      end do

      ! --- polarization ---------------------------------------------------------
      if (talk) call report(stage, "distributed multipoles", talk)

      allocate (pot%frequencies(N_CASIMIR_POLDER))
      pot%frequencies = casimir_polder_frequencies()
      ! One Hessian for all three dynamic blocks: it depends on the reference
      ! alone, so rebuilding it per block would be three identical builds. Only
      ! this call can build it -- the two blocks after are handed the built one --
      ! so this is the one place the auxiliary basis has to reach, and an optional
      ! dummy cannot be passed conditionally from a local, hence the branch.
      !
      ! Ahead of the static block, which is solved inside this call: the static
      ! response is the zero-frequency member of the family these blocks solve,
      ! so it rides along as one more frequency and comes back in `u_static`,
      ! on either route. Matrix free, that is the difference between one
      ! batched solve and two, and the second one used to run at the solver's
      ! own defaults rather than at the deck's.
      if (present(aux_basis)) then
         call dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, pot, &
                                      shared_hessian, error, progress=talk, aux=aux, &
                                      max_iter=dynamic_maxiter, tol=dynamic_tol, &
                                      route=response, &
                                      allow_unconverged=allow_crap_response, &
                                      batch=response_batch, static_response=u_static, &
                                      b_ao=b_ao)
      else
         call dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, pot, &
                                      shared_hessian, error, progress=talk, &
                                      max_iter=dynamic_maxiter, tol=dynamic_tol, &
                                      route=response, &
                                      allow_unconverged=allow_crap_response, &
                                      batch=response_batch, static_response=u_static)
      end if
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) then
         write (line, "(A,I0,A)") "  polarizabilities at ", N_CASIMIR_POLDER, " imaginary frequencies"
         call logger%info(trim(line))
         write (line, "(A)") "  dipole-quadrupole dispersion tensors"
         call logger%info(trim(line))
      end if
      if (talk) call report(stage, "all three dynamic blocks, with the Hessian build", talk)

      ! The static block, from the response the dynamic blocks already solved
      ! for; only the localization and the contraction are left to do here.
      call distributed_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                      pot%n_occ, pot%static_pol, pot%centroids, &
                                      error, n_core=core, response=u_static)
      deallocate (u_static)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      if (talk) call report(stage, "static polarizability", talk)

      ! --- exchange repulsion: the LMO Fock matrix and the orbitals themselves --
      ! F in the LMO basis is W^T diag(eps) W with W = C_occ^T S C_loc, so no AO
      ! Fock matrix is needed and nothing here depends on basis function ordering.
      call mol%overlap(ovl)
      allocate (sc(mol%nao, n_valence), w(pot%n_occ, n_valence))
      call pic_gemm(ovl, loc, sc)
      call pic_gemm(scf%orbitals(:, 1:pot%n_occ), sc, w, transa="T")
      allocate (scaled(pot%n_occ, n_valence), pot%fock_lmo(n_valence, n_valence))
      do j = 1, n_valence
         do k = 1, pot%n_occ
            scaled(k, j) = scf%orbital_energies(k)*w(k, j)
         end do
      end do
      call pic_gemm(w, scaled, pot%fock_lmo, transa="T")

      call to_gamess_ao_order(mol, loc, pot%orbitals, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! The charge-transfer basis, in the canonical form: `$MAKEFP CTVVO=.FALSE.`
      ! writes the whole canonical MO matrix under the header `CTVEC NA NUM`,
      ! where GAMESS's default path writes `NOCC` occupied orbitals plus
      ! quasi-atomic valence virtuals from `VVOS`. The canonical form needs no
      ! extra machinery and is what GAMESS recommends when the valence virtuals
      ! cannot be formed.
      call to_gamess_ao_order(mol, scf%orbitals, pot%canonical, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      call projection_basis_lines(atomic_numbers, pot%labels, pot%points, &
                                  element_symbols, basis_name, pot%basis_lines, &
                                  error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      ! --- charge penetration screening ----------------------------------------
      ! Last, because it is fitted to the error the multipoles above make.
      !
      ! One grid and one quantum potential for both damping forms: they are fitted
      ! to the same target and differ only in the damping term, so the first fit
      ! hands its target to the second. Building the target is nearly the whole
      ! cost -- tens of thousands of grid points, each needing an integral over
      ! every shell pair.
      call fit_screening(mol, scf%density, dma, atomic_numbers, SCREEN_EXPONENTIAL, &
                         alpha, error, target=screen_target, residual=rms_exp, &
                         vdw_scale=pot%vdwscl)
      if (error%has_error()) then
         call screen_target%destroy()
         call mol%destroy()
         return
      end if
      allocate (pot%screen2(pot%n_points))
      pot%screen2 = alpha
      deallocate (alpha)
      call fit_screening(mol, scf%density, dma, atomic_numbers, SCREEN_GAUSSIAN, &
                         alpha, error, target=screen_target, residual=rms_gauss, &
                         vdw_scale=pot%vdwscl)
      call screen_target%destroy()
      if (error%has_error()) then
         call mol%destroy()
         return
      end if
      allocate (pot%screen(pot%n_points))
      pot%screen = alpha
      deallocate (alpha)
      if (talk) then
         write (line, "(A,F0.4,A,F0.4,A)") "  screening fitted: exponential misses by ", &
            rms_exp, " kcal/mol, Gaussian by ", rms_gauss, " kcal/mol"
         call logger%info(trim(line))
      end if
      if (talk) call report(stage, "charge-penetration screening", talk)

      call shared_hessian%destroy()
      call mol%destroy()
      deallocate (loc, ovl, sc, w, scaled)
   end subroutine make_efp_potential

   subroutine dipole_quadrupole_block(mol, scf, coordinates, atomic_numbers, core, &
                                      pot, hessian, error, progress, aux, &
                                      max_iter, tol, route, allow_unconverged, batch, &
                                      static_response, b_ao)
      !! `DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS`, ready to write
      !!
      !! Three conventions here were established by
      !! `validation/check_dipquad_sumrule`, which pins them by structure rather
      !! than by fitting. **Getting any one of them wrong leaves a tensor that
      !! passes every internal check and disagrees with GAMESS.**
      !!
      !!   * **The quadrupole measures and the dipole drives**, both expanded about
      !!     the centre of mass, and the quadrupole is the traceless Buckingham
      !!     form. Per orbital that is not the same as the reverse assignment,
      !!     because the projector onto the localized set does not commute with the
      !!     response operator; summed over orbitals they agree.
      !!   * **The write-time translation to each centroid**, `DQSHIFT`, whose
      !!     `delta_bc` term takes the dipole-dipole tensor **transposed**:
      !!     `alpha(a,d)`, not `alpha(d,a)` as the rest of the formula reads.
      !!   * **The dipole-dipole tensor in the shift is the dynamic one at the same
      !!     frequency**, not the static one.
      type(czt_molecule_t), intent(in) :: mol
      type(rhf_result_t), intent(in) :: scf
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: core
      type(efp_potential_t), intent(inout) :: pot
      type(response_hessian_t), intent(inout) :: hessian
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: progress
      type(czt_molecule_t), intent(in), optional :: aux
         !! Fit the Hessian rather than build it exactly; passed straight through.
      integer, intent(in), optional :: max_iter
         !! This and the four below are the `keywords.efp` settings the response
         !! solve reads, passed straight through. Absent here means absent there,
         !! which leaves the solver on its own defaults.
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: route
      logical, intent(in), optional :: allow_unconverged
      integer, intent(in), optional :: batch
      real(dp), allocatable, intent(out), optional :: static_response(:, :, :)
         !! `(n_vir, n_occ, 3)`: the static dipole response over the whole
         !! occupied space, in `cphf_solve`'s convention, solved as the
         !! zero-frequency member of the same batch. What
         !! `distributed_polarizability` takes as `response`.
      real(dp), allocatable, intent(inout), optional :: b_ao(:, :)
         !! The fitted SCF's AO tensor, for the response to transform rather
         !! than rebuild. Consumed by the solve.

      real(dp), allocatable :: dip(:, :, :), quad(:, :, :), buck(:, :, :)
      real(dp), allocatable :: with_zero(:), u_all(:, :, :)
      real(dp), allocatable :: both(:, :, :), all_blocks(:, :, :, :)
      real(dp), allocatable :: drives(:, :, :), solved(:, :, :, :)
      integer :: drive_of(N_CART_PAIR)
      real(dp), allocatable :: raw(:, :, :, :), centroids(:, :), qq(:, :, :, :)
      real(dp) :: com(3), r(3), alpha(3, 3)
      real(dp) :: mass_total, isotropic
      integer :: i, a, b, c, d, k, f, n_freq, n_both

      com = 0.0_dp
      mass_total = 0.0_dp
      do i = 1, size(atomic_numbers)
         com = com + element_mass(atomic_numbers(i))*coordinates(:, i)
         mass_total = mass_total + element_mass(atomic_numbers(i))
      end do
      com = com/mass_total

      call multipole_matrices(mol, com, 1, dip, error)
      if (error%has_error()) return
      call multipole_matrices(mol, com, 2, quad, error)
      if (error%has_error()) return

      ! The traceless Buckingham quadrupole, as GAMESS builds it, kept as all nine
      ! Cartesian slots so that no expansion of six unique values into nine has to
      ! be guessed at.
      allocate (buck(mol%nao, mol%nao, 9))
      buck(:, :, QXX) = 0.5_dp*(2.0_dp*quad(:, :, QXX) - quad(:, :, QYY) &
                                - quad(:, :, QZZ))
      buck(:, :, QYY) = 0.5_dp*(2.0_dp*quad(:, :, QYY) - quad(:, :, QXX) &
                                - quad(:, :, QZZ))
      buck(:, :, QZZ) = 0.5_dp*(2.0_dp*quad(:, :, QZZ) - quad(:, :, QXX) &
                                - quad(:, :, QYY))
      buck(:, :, QXY) = 1.5_dp*quad(:, :, QXY)
      buck(:, :, QXZ) = 1.5_dp*quad(:, :, QXZ)
      buck(:, :, QYZ) = 1.5_dp*quad(:, :, QYZ)
      buck(:, :, QYX) = buck(:, :, QXY)
      buck(:, :, QZX) = buck(:, :, QXZ)
      buck(:, :, QZY) = buck(:, :, QYZ)

      ! All three dynamic blocks from one solve: they are three contractions of
      ! one response, not three responses. The dipole-dipole block measures and
      ! drives with `dip`, the mixed one measures with `buck` and drives with
      ! `dip`, and the quadrupole one uses `buck` on both sides -- so the driving
      ! operators are `dip` and `buck` together, twelve of them, and every block
      ! is a slice of the twelve-by-twelve result.
      ! Dipole only: three operators measure and three drive, and the two
      ! quadrupole blocks are never formed.
      if (pot%quadrupole_blocks) then
         n_both = 3 + size(buck, 3)
      else
         n_both = 3
      end if
      allocate (both(mol%nao, mol%nao, n_both))
      both(:, :, 1:3) = dip
      if (pot%quadrupole_blocks) both(:, :, 4:n_both) = buck

      ! Twelve operators measure, eight drive. The nine quadrupole slots hold
      ! only five independent operators -- three are their own transposes and
      ! the traceless form fixes `zz` as `-(xx + yy)` -- and a response solve is
      ! linear in its driving operator, so the responses to the other four are
      ! sums of responses already in hand. Solving them was a third of every
      ! matrix-free response; `drive_of` says which solved column each of the
      ! nine slots reads, with `zz` read as the negative of two.
      allocate (drives(mol%nao, mol%nao, merge(8, 3, pot%quadrupole_blocks)))
      drives(:, :, 1:3) = dip
      if (pot%quadrupole_blocks) then
         drives(:, :, 4) = buck(:, :, QXX)
         drives(:, :, 5) = buck(:, :, QYY)
         drives(:, :, 6) = buck(:, :, QXY)
         drives(:, :, 7) = buck(:, :, QXZ)
         drives(:, :, 8) = buck(:, :, QYZ)
         drive_of = 0
         drive_of(QXX) = 4
         drive_of(QYY) = 5
         drive_of(QXY) = 6
         drive_of(QYX) = 6
         drive_of(QXZ) = 7
         drive_of(QZX) = 7
         drive_of(QYZ) = 8
         drive_of(QZY) = 8
      end if

      ! The static response is the zero-frequency member of the same family,
      ! so it is solved here as a thirteenth frequency, appended after the
      ! Casimir-Polder twelve so that every slice below keeps its indexing,
      ! and handed out through `static_response` rather than solved again.
      n_freq = size(pot%frequencies)
      allocate (with_zero(n_freq + 1))
      with_zero(1:n_freq) = pot%frequencies
      with_zero(n_freq + 1) = 0.0_dp
      ! One call rather than one per combination of present arguments: `aux` here
      ! is an optional dummy, not a local, and an absent one passed on as an actual
      ! argument arrives absent at the other end.
      call distributed_dynamic_cross(mol, scf%orbitals, scf%orbital_energies, &
                                     pot%n_occ, with_zero, both, drives, &
                                     solved, centroids, error, n_core=core, &
                                     hessian=hessian, progress=progress, aux=aux, &
                                     max_iter=max_iter, tol=tol, route=route, &
                                     allow_unconverged=allow_unconverged, batch=batch, &
                                     static_response=u_all, b_ao=b_ao)
      if (error%has_error()) return
      deallocate (both, drives, with_zero)

      ! The twelve-by-twelve table every slice below reads, the four dependent
      ! drive columns filled by linearity.
      allocate (all_blocks(n_both, n_both, size(solved, 3), size(solved, 4)))
      all_blocks(:, 1:3, :, :) = solved(:, 1:3, :, :)
      if (pot%quadrupole_blocks) then
         do d = 1, 9
            if (drive_of(d) > 0) then
               all_blocks(:, 3 + d, :, :) = solved(:, drive_of(d), :, :)
            else
               all_blocks(:, 3 + d, :, :) = -solved(:, drive_of(QXX), :, :) &
                                            - solved(:, drive_of(QYY), :, :)
            end if
         end do
      end if
      deallocate (solved)
      if (present(static_response)) then
         allocate (static_response, source=u_all(:, :, 1:3))
      end if
      deallocate (u_all)

      if (allocated(pot%centroids)) deallocate (pot%centroids)
      allocate (pot%centroids, source=centroids)
      allocate (pot%dynamic_pol(3, 3, pot%n_lmo, n_freq))
      pot%dynamic_pol = all_blocks(1:3, 1:3, :, 1:n_freq)
      if (.not. pot%quadrupole_blocks) then
         deallocate (all_blocks)
         return
      end if

      allocate (raw(size(buck, 3), 3, pot%n_lmo, n_freq))
      raw = all_blocks(4:n_both, 1:3, :, 1:n_freq)
      allocate (pot%dipquad(3, 3, 3, pot%n_lmo, n_freq))
      allocate (pot%dipquad_pre(3, 3, 3, pot%n_lmo, n_freq))
      do f = 1, n_freq
         do k = 1, pot%n_lmo
            r = pot%centroids(:, k) - com
            alpha = pot%dynamic_pol(:, :, k, f)
            do a = 1, 3
               ! The delta_bc term's transpose: alpha(a,d), not alpha(d,a).
               isotropic = 0.0_dp
               do d = 1, 3
                  isotropic = isotropic + r(d)*alpha(a, d)
               end do
               do b = 1, 3
                  do c = 1, 3
                     ! raw is (quadrupole slot, dipole, orbital, frequency); the
                     ! nine quadrupole slots run with the second index fastest.
                     pot%dipquad_pre(a, b, c, k, f) = raw((b - 1)*3 + c, a, k, f)
                     pot%dipquad(a, b, c, k, f) = raw((b - 1)*3 + c, a, k, f) &
                                                  - 1.5_dp*(r(b)*alpha(c, a) &
                                                            + r(c)*alpha(a, b))
                     if (b == c) then
                        pot%dipquad(a, b, c, k, f) = pot%dipquad(a, b, c, k, f) &
                                                     + isotropic
                     end if
                  end do
               end do
            end do
         end do
      end do

      ! --- the quadrupole-quadrupole block ------------------------------------
      ! Same operator on both sides, and the factor is 1/3 rather than 1:
      ! `LQQPOL` carries a factor of a third that the dipole-quadrupole one does
      ! not. Confirmed against GAMESS's own molecular `QUAD-QUAD POLARIZABILITY`,
      ! which `$MAKEFP MOLPOL=.TRUE.` writes with no translation applied.
      allocate (qq(size(buck, 3), size(buck, 3), pot%n_lmo, n_freq))
      qq = all_blocks(4:n_both, 4:n_both, :, 1:n_freq)
      deallocate (all_blocks)

      allocate (pot%quadquad(3, 3, 3, 3, pot%n_lmo, n_freq))
      do f = 1, n_freq
         do k = 1, pot%n_lmo
            r = pot%centroids(:, k) - com
            alpha = pot%dynamic_pol(:, :, k, f)
            call qq_shift(qq(:, :, k, f)/3.0_dp, pot%dipquad_pre(:, :, :, k, f), &
                          alpha, r, pot%quadquad(:, :, :, :, k, f))
         end do
      end do

      deallocate (dip, quad, buck, raw, qq, centroids)
   end subroutine dipole_quadrupole_block

   subroutine qq_shift(qq, dq, alpha, r, shifted)
      !! `QQSHIFT`, the quadrupole-quadrupole translation
      !!
      !! Transcribed term for term. It mixes in both the dipole-dipole and the
      !! *pre-shift* dipole-quadrupole tensors, which is why the two blocks are
      !! built together.
      real(dp), intent(in) :: qq(N_CART_PAIR, N_CART_PAIR)       !! pre-shift, already scaled by 1/3
      real(dp), intent(in) :: dq(3, 3, 3)    !! pre-shift dipole-quadrupole
      real(dp), intent(in) :: alpha(3, 3)
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: shifted(3, 3, 3, 3)

      real(dp) :: a1, a2, rralph, adelt1, adelt2, rrad1, rrad2, rrdalph
      integer :: a, b, c, e, i, j

      do a = 1, 3
         do b = 1, 3
            do c = 1, 3
               do e = 1, 3
                  a1 = r(a)*dq(b, c, e) + r(b)*dq(a, c, e)
                  a2 = r(c)*dq(e, a, b) + r(e)*dq(c, a, b)
                  rralph = r(a)*r(c)*alpha(b, e) + r(a)*r(e)*alpha(b, c) &
                           + r(b)*r(c)*alpha(a, e) + r(b)*r(e)*alpha(a, c)
                  adelt1 = 0.0_dp
                  adelt2 = 0.0_dp
                  rrad1 = 0.0_dp
                  rrad2 = 0.0_dp
                  rrdalph = 0.0_dp
                  do i = 1, 3
                     if (a == b) adelt1 = adelt1 + r(i)*dq(i, c, e)
                     if (c == e) adelt2 = adelt2 + r(i)*dq(i, a, b)
                     if (a == b) rrad1 = rrad1 + r(c)*r(i)*alpha(i, e) &
                                         + r(e)*r(i)*alpha(i, c)
                     if (c == e) rrad2 = rrad2 + r(a)*r(i)*alpha(b, i) &
                                         + r(b)*r(i)*alpha(a, i)
                     if (a == b .and. c == e) then
                        do j = 1, 3
                           rrdalph = rrdalph + r(i)*r(j)*alpha(i, j)
                        end do
                     end if
                  end do
                  shifted(a, b, c, e) = qq((a - 1)*3 + b, (c - 1)*3 + e) &
                                        - 0.5_dp*(a1 + a2 + rrad1 + rrad2) &
                                        + (adelt1 + adelt2 + rrdalph)/3.0_dp &
                                        + 0.75_dp*rralph
               end do
            end do
         end do
      end do
   end subroutine qq_shift

   pure function frozen_core(atomic_numbers) result(n)
      !! The standard frozen core, which is the set MAKEFP excludes
      ! TODO(mqc): disagrees with `core_orbital_count` in `mqc_czt_bridge`
      ! above xenon -- that one adds the 4d 5s 5p shells for Z > 54 where this
      ! stays at 18 -- and both are documented as the standard frozen core.
      integer, intent(in) :: atomic_numbers(:)
      integer :: n
      integer :: i, z

      n = 0
      do i = 1, size(atomic_numbers)
         z = atomic_numbers(i)
         if (z > 2 .and. z <= 10) then
            n = n + 1
         else if (z > 10 .and. z <= 18) then
            n = n + 5
         else if (z > 18 .and. z <= 36) then
            n = n + 9
         else if (z > 36) then
            n = n + 18
         end if
      end do
   end function frozen_core

   subroutine report(stage, what, talk)
      !! Seconds for the stage just finished, then restart the clock
      !!
      !! The dynamic response dominates a mid-sized fragment, but the split moves
      !! with the fragment: the localization is quadratic in the occupied count
      !! and the multipoles work over primitive pairs, so neither can be assumed
      !! small for something larger.
      type(timer_type), intent(inout) :: stage
      character(len=*), intent(in) :: what
      logical, intent(in) :: talk

      character(len=MAX_LINE_LENGTH) :: line

      if (.not. talk) return
      call stage%stop()
      write (line, "(A,F9.1,A,A)") "      ", stage%get_elapsed_time(), " s  ", what
      call logger%info(trim(line))
      flush (6)
      call stage%start()
   end subroutine report

   subroutine check_angular_form(mol, error)
      !! Refuse a basis whose angular form this cannot write
      !!
      !! **The ordering map is Cartesian s, p, d and f, and both codes have to
      !! agree on that.** GAMESS defaults ISPHER = -1, Cartesian, which is what a
      !! Pople set wants.
      !!
      !! A Dunning or def2 set is not simply unsupported: the Basis Set Exchange
      !! declares those spherical and GAMESS will not run them Cartesian, so both
      !! sides would agree on spherical. What is missing is the spherical ordering
      !! map and `ISPHER=1` in the deck.
      ! TODO(mqc): the message below still says only s, p and d are mapped, while
      ! `to_gamess_ao_order` has mapped Cartesian f since `F_FROM_LIBCINT` was
      ! added and refuses only g and higher.
      type(czt_molecule_t), intent(in) :: mol
      type(error_t), intent(inout) :: error

      integer :: ish
      logical :: has_high_l

      ! Only an issue if there is actually a shell it applies to: s and p are the
      ! same either way, so 6-31G* on hydrogen alone is not a spherical basis in
      ! any sense that matters.
      has_high_l = .false.
      do ish = 1, mol%nbas
         if (mol%bas(LIBCINT_ANG_OF, ish) >= 2) has_high_l = .true.
      end do
      if (.not. mol%cartesian .and. has_high_l) then
         call error%set(ERROR_VALIDATION, &
                        "makefp: this basis is spherical and only Cartesian s, p "// &
                        "and d are mapped to GAMESS's ordering. GAMESS reads a "// &
                        "spherical potential with ISPHER=1, so this needs a "// &
                        "spherical ordering map rather than a Cartesian basis.")
      end if
   end subroutine check_angular_form

   subroutine from_gamess_ao_order(mol, mapped, coefficients, error)
      !! The inverse of `to_gamess_ao_order`
      !!
      !! A potential stores its orbitals in GAMESS's AO order, so anything reading
      !! one back for use with our own integrals has to undo the permutation and
      !! the normalization.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: mapped(:, :)
      real(dp), allocatable, intent(out) :: coefficients(:, :)
      type(error_t), intent(inout) :: error

      integer :: ish, l, off, dim, slot
      real(dp) :: scale

      allocate (coefficients(size(mapped, 1), size(mapped, 2)))
      coefficients = mapped
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         off = mol%shell_offset(ish)
         dim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         if (l == 2 .and. dim == 6) then
            do slot = 1, 6
               scale = D_NORMALIZATION
               if (slot > 3) scale = D_NORMALIZATION/sqrt(3.0_dp)
               coefficients(off + D_FROM_LIBCINT(slot), :) = mapped(off + slot, :)/scale
            end do
         else if (l == 3 .and. dim == 10) then
            do slot = 1, 10
               select case (F_CLASS(slot))
               case (1)
                  scale = F_NORMALIZATION
               case (2)
                  scale = F_NORMALIZATION/sqrt(5.0_dp)
               case default
                  scale = F_NORMALIZATION/sqrt(15.0_dp)
               end select
               coefficients(off + F_FROM_LIBCINT(slot), :) = mapped(off + slot, :)/scale
            end do
         else if (l >= 2) then
            call error%set(ERROR_VALIDATION, "efp: this basis has g functions or "// &
                           "higher, whose ordering against GAMESS is not mapped")
            return
         end if
      end do
   end subroutine from_gamess_ao_order

   subroutine to_gamess_ao_order(mol, coefficients, mapped, error)
      !! Orbital coefficients in the AO order and normalization GAMESS reads
      !!
      !! Only the Cartesian d and f shells move. Our s and p already agree,
      !! including the interleaving inside a shared-exponent `L` shell.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coefficients(:, :)
      real(dp), allocatable, intent(out) :: mapped(:, :)
      type(error_t), intent(inout) :: error

      integer :: ish, l, off, dim, slot
      real(dp) :: scale

      allocate (mapped(size(coefficients, 1), size(coefficients, 2)))
      mapped = coefficients
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         off = mol%shell_offset(ish)
         dim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         if (l == 2 .and. dim == 6) then
            do slot = 1, 6
               scale = D_NORMALIZATION
               if (slot > 3) scale = D_NORMALIZATION/sqrt(3.0_dp)
               mapped(off + slot, :) = coefficients(off + D_FROM_LIBCINT(slot), :)*scale
            end do
         else if (l == 3 .and. dim == 10) then
            do slot = 1, 10
               select case (F_CLASS(slot))
               case (1)
                  scale = F_NORMALIZATION
               case (2)
                  scale = F_NORMALIZATION/sqrt(5.0_dp)
               case default
                  scale = F_NORMALIZATION/sqrt(15.0_dp)
               end select
               mapped(off + slot, :) = coefficients(off + F_FROM_LIBCINT(slot), :)*scale
            end do
         else if (l >= 2) then
            ! g and up need their own permutation and normalizations against
            ! GAMESS's ordering, derived the way the d and f ones were: read off
            ! its own printed coefficients for a basis that has them.
            call error%set(ERROR_VALIDATION, &
                           "makefp: this basis has g functions or higher, and only "// &
                           "Cartesian s, p, d and f are mapped to GAMESS's ordering "// &
                           "so far")
            return
         end if
      end do
   end subroutine to_gamess_ao_order

   subroutine projection_basis_lines(atomic_numbers, labels, points, symbols, &
                                     basis_name, lines, error)
      !! `PROJECTION BASIS SET`, in GAMESS's columns and its normalization
      !!
      !! Two conventions are GAMESS's rather than ours. Shells are named the way
      !! it names them, including `L` for a shared-exponent sp pair -- our own
      !! reader emits one shell per coefficient column, so an `L` has to be
      !! recognised by finding an s and a p over identical exponents. And the
      !! printed coefficient has the primitive normalization folded in, which is
      !! the factor `gamess_primitive_norm` supplies.
      !!
      !! The primitive counter runs across the whole file rather than restarting
      !! per atom.
      integer, intent(in) :: atomic_numbers(:)
      character(len=8), intent(in) :: labels(:)
      real(dp), intent(in) :: points(:, :)
      character(len=*), intent(in) :: symbols(:)
      character(len=*), intent(in) :: basis_name
      character(len=MAX_LINE), allocatable, intent(out) :: lines(:)
      type(error_t), intent(inout) :: error

      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path
      character(len=MAX_LINE), allocatable :: buffer(:)
      character(len=MAX_LINE) :: text
      integer :: n, natm, iatom, ish, nsh, k, primitive, valence, ncol, col
      integer :: l_of(2)
      logical :: is_l
      real(dp) :: expo, coefficient

      natm = size(atomic_numbers)
      call find_basis_file(basis_name, path, error)
      if (error%has_error()) return
      call build_molecular_basis_json(path, symbols, basis, error)
      if (error%has_error()) return

      allocate (buffer(4096))
      n = 0
      primitive = 0
      do iatom = 1, natm
         valence = atomic_numbers(iatom) - frozen_core([atomic_numbers(iatom)])*2
         ! Ten wide here, where COORDINATES uses eight -- what MAKEFP writes, and
         ! the two sections genuinely differ. Written as A8 then two blanks rather
         ! than A10, because A10 on an eight-character string pads on the left and
         ! would indent the label instead of the number.
         write (text, "(A8,A2,3F15.10,F7.1)") labels(iatom), "  ", &
            points(:, iatom), real(valence, dp)
         n = n + 1
         buffer(n) = text

         nsh = basis%elements(iatom)%nshells
         ish = 1
         do while (ish <= nsh)
            ! An sp pair arrives as two consecutive shells over the same
            ! exponents; GAMESS wants them as one "L" with two coefficient
            ! columns.
            is_l = .false.
            if (ish < nsh) then
               if (basis%elements(iatom)%shells(ish)%ang_mom == 0 .and. &
                   basis%elements(iatom)%shells(ish + 1)%ang_mom == 1) then
                  if (size(basis%elements(iatom)%shells(ish)%exponents) == &
                      size(basis%elements(iatom)%shells(ish + 1)%exponents)) then
                     is_l = maxval(abs(basis%elements(iatom)%shells(ish)%exponents - &
                                       basis%elements(iatom)%shells(ish + 1)%exponents)) &
                            < 1.0e-12_dp
                  end if
               end if
            end if

            if (is_l) then
               ncol = 2
               l_of = [0, 1]
               write (text, "(A,I11)") "   L", &
                  size(basis%elements(iatom)%shells(ish)%exponents)
            else
               ncol = 1
               l_of(1) = basis%elements(iatom)%shells(ish)%ang_mom
               write (text, "(A,A,I11)") "   ", &
                  shell_letter(l_of(1)), size(basis%elements(iatom)%shells(ish)%exponents)
            end if
            n = n + 1
            buffer(n) = text

            do k = 1, size(basis%elements(iatom)%shells(ish)%exponents)
               primitive = primitive + 1
               expo = basis%elements(iatom)%shells(ish)%exponents(k)
               write (text, "(I6,F21.10)") primitive, expo
               do col = 1, ncol
                  coefficient = basis%elements(iatom)%shells(ish + col - 1) &
                                %coefficients(k)
                  write (text(len_trim(text) + 1:), "(F15.8)") &
                     coefficient*gamess_primitive_norm(l_of(col), expo)
               end do
               n = n + 1
               buffer(n) = text
            end do
            ish = ish + ncol
         end do
         n = n + 1
         buffer(n) = "  "
      end do

      allocate (lines(n))
      lines = buffer(1:n)
      deallocate (buffer)
      call basis%destroy()
   end subroutine projection_basis_lines

   pure function shell_letter(l) result(c)
      integer, intent(in) :: l
      character(len=1) :: c
      character(len=*), parameter :: LETTERS = "SPDFGH"

      if (l >= 0 .and. l < len(LETTERS)) then
         c = LETTERS(l + 1:l + 1)
      else
         c = "?"
      end if
   end function shell_letter

   pure function gamess_primitive_norm(l, exponent) result(factor)
      !! The normalization GAMESS folds into a printed contraction coefficient
      integer, intent(in) :: l
      real(dp), intent(in) :: exponent
      real(dp) :: factor
      real(dp) :: double_factorial
      integer :: n

      double_factorial = 1.0_dp
      n = 2*l - 1
      do while (n > 1)
         double_factorial = double_factorial*real(n, dp)
         n = n - 2
      end do
      factor = (2.0_dp*exponent/PI)**0.75_dp*(4.0_dp*exponent)**(0.5_dp*real(l, dp)) &
               /sqrt(double_factorial)
   end function gamess_primitive_norm

   subroutine write_efp_potential(pot, path, error, omitted)
      !! The potential as a `.efp` file
      !!
      !! `omitted` comes back with the sections a complete MAKEFP would have
      !! written and this did not, so a caller can say so rather than let their
      !! absence pass unremarked.
      type(efp_potential_t), intent(in) :: pot
      character(len=*), intent(in) :: path
      type(error_t), intent(inout) :: error
      character(len=:), allocatable, intent(out), optional :: omitted

      integer :: unit, i, k, f, stat, a, b, c, e
      character(len=8) :: label
      real(dp) :: tensor(N_CART_PAIR)
      real(dp) :: wide(N_CART_TRIPLE)
      real(dp) :: broad(N_CART_QUAD)

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call error%set(ERROR_VALIDATION, "makefp: cannot write "//trim(path))
         return
      end if

      write (unit, "(A)") &
         "          RUNTYP=MAKEFP EFFECTIVE FRAGMENT POTENTIAL DATA FOLLOWS..."
      write (unit, "(A)") "          "//pot%name//" GENERATED BY METALQUICHA"
      write (unit, "(A)") " $"//pot%name
      write (unit, "(A)") "EFP DATA FOR "//pot%name//" SCFTYP=RHF     "// &
         "... GENERATED WITH BASIS SET="//pot%basis_name

      write (unit, "(A)") " COORDINATES (BOHR)"
      do i = 1, pot%n_points
         write (unit, "(A8,3F15.10,F12.7,F5.1)") pot%labels(i), pot%points(:, i), &
            pot%mass(i), pot%charge(i)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " MONOPOLES "
      do i = 1, pot%n_points
         write (unit, "(A8,F15.10,F10.5)") pot%labels(i), pot%q_elec(i), pot%q_nuc(i)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " DIPOLES "
      do i = 1, pot%n_points
         call write_record(unit, pot%labels(i), pot%dipole(:, i), 16, 10, 3)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " QUADRUPOLES "
      do i = 1, pot%n_points
         call write_record(unit, pot%labels(i), pot%quadrupole(:, i), 16, 10, 4)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " OCTUPOLES  "
      do i = 1, pot%n_points
         call write_record(unit, pot%labels(i), pot%octopole(:, i), 17, 9, 4)
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A)") " POLARIZABLE POINTS"
      do k = 1, pot%n_lmo
         write (label, "(A,I0)") "CT", k
         do i = 1, 9
            tensor(i) = pot%static_pol(POL_ROW(i), POL_COL(i), k)
         end do
         call write_tensor_point(unit, label, pot%centroids(:, k), tensor)
      end do
      write (unit, "(A)") " STOP"

      ! The dynamic blocks run frequency-outermost, each stamped with its own
      ! frequency, and the label carries a gap -- "CT  2" where every other
      ! section writes "CT2".
      write (unit, "(A)") " DYNAMIC POLARIZABLE POINTS"
      do f = 1, size(pot%frequencies)
         do k = 1, pot%n_lmo
            write (label, "(A,I3)") "CT", k
            do i = 1, 9
               tensor(i) = pot%dynamic_pol(POL_ROW(i), POL_COL(i), k, f)
            end do
            ! Only the first point of a block carries the frequency, which is how
            ! MAKEFP writes it: the stamp opens a block rather than labelling a
            ! point. GAMESS's reader tolerates one per point, but a parser keying
            ! on the stamp -- ours included -- would then see every point as a new
            ! block.
            if (k == 1) then
               call write_tensor_point(unit, label, pot%centroids(:, k), tensor, &
                                       frequency=pot%frequencies(f))
            else
               call write_tensor_point(unit, label, pot%centroids(:, k), tensor)
            end if
         end do
      end do
      write (unit, "(A)") " STOP"

      if (allocated(pot%dipquad)) then
         ! The dipole-quadrupole block, 27 values a point. The slot order is
         ! `(a-1)*9 + (c-1)*3 + b` -- the *first* quadrupole index runs fastest, which
         ! is transposed from how the `DQSHIFT` source reads, and was pinned by
         ! requiring the pre-shift tensor's symmetry in `bc` to come back.
         write (unit, "(A)") " DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS"
         do f = 1, size(pot%frequencies)
            do k = 1, pot%n_lmo
               write (label, "(A,I3)") "CT", k
               do a = 1, 3
                  do b = 1, 3
                     do c = 1, 3
                        wide((a - 1)*9 + (c - 1)*3 + b) = pot%dipquad(a, b, c, k, f)
                     end do
                  end do
               end do
               if (k == 1) then
                  write (unit, "(A,3F15.10,A,F9.6,A)") trim(label), &
                     pot%centroids(:, k), " -- FOR W=", pot%frequencies(f), "I A.U."
               else
                  write (unit, "(A,3F15.10)") trim(label), pot%centroids(:, k)
               end if
               call write_values(unit, wide, 16, 10, 4)
            end do
         end do
         write (unit, "(A)") " STOP"
      end if

      if (allocated(pot%quadquad)) then
         ! The quadrupole-quadrupole block, 81 values a point, written with the last
         ! index fastest, the last of the four varying first. No
         ! transposition here, unlike the dipole-quadrupole slots: every `QQSHIFT`
         ! term is symmetric within each index pair, so the written values are too.
         write (unit, "(A)") " LMOQQPOL DYNAMIC POLARIZABLE POINTS"
         do f = 1, size(pot%frequencies)
            do k = 1, pot%n_lmo
               write (label, "(A,I3)") "CT", k
               i = 0
               do a = 1, 3
                  do b = 1, 3
                     do c = 1, 3
                        do e = 1, 3
                           i = i + 1
                           broad(i) = pot%quadquad(a, b, c, e, k, f)
                        end do
                     end do
                  end do
               end do
               if (k == 1) then
                  write (unit, "(A,3F15.10,A,F9.6,A)") trim(label), &
                     pot%centroids(:, k), " -- FOR W=", pot%frequencies(f), "I A.U."
               else
                  write (unit, "(A,3F15.10)") trim(label), pot%centroids(:, k)
               end if
               call write_values(unit, broad, 16, 10, 4)
            end do
         end do
         write (unit, "(A)") " STOP"
      end if

      write (unit, "(A)") " PROJECTION BASIS SET"
      do i = 1, size(pot%basis_lines)
         write (unit, "(A)") trim(pot%basis_lines(i))
      end do
      write (unit, "(A)") " STOP"

      write (unit, "(A,I5)") " MULTIPLICITY", pot%multiplicity
      write (unit, "(A)") " STOP"

      write (unit, "(A,2I7)") " PROJECTION WAVEFUNCTION", pot%n_lmo, pot%nao
      call write_wavefunction(unit, pot%orbitals)

      write (unit, "(A)") " FOCK MATRIX ELEMENTS"
      call write_lower_triangle(unit, pot%fock_lmo)

      write (unit, "(A)") " LMO CENTROIDS"
      do k = 1, pot%n_lmo
         write (label, "(A,I0)") "CT", k
         write (unit, "(A3,3F15.10)") adjustl(label), pot%centroids(:, k)
      end do
      write (unit, "(A)") " STOP"

      ! Charge transfer, in the canonical-orbital form: the header carries the
      ! occupied count and the number of vectors, then the whole MO matrix in the
      ! same five-to-a-line layout the projection wavefunction uses. `CTFOK` is a
      ! *subsection* of this one, not a section of its own -- GAMESS's reader looks
      ! for it only directly behind a `CTVEC` block and aborts on a standalone one --
      ! so the two are written together and share one `STOP`.
      write (unit, "(A,I8,A,I8)") " CTVEC   ", pot%n_occ, "  ", pot%nao
      call write_wavefunction(unit, pot%canonical)
      write (unit, "(A)") " CTFOK   "
      call write_values(unit, pot%eps_occ, 16, 10, 4)
      write (unit, "(A)") " STOP"

      ! Beta is frozen at one, as MAKEFP freezes it: its ICFIX flag fixes the
      ! prefactor and fits the exponent alone.
      write (unit, "(A,F8.3,A)") "SCREEN2      (FROM VDWSCL=", pot%vdwscl, ")"
      do i = 1, pot%n_points
         write (unit, "(1X,A8,2F14.9)") pot%labels(i), 1.0_dp, pot%screen2(i)
      end do
      write (unit, "(A)") "STOP"

      write (unit, "(A,F8.3,A)") "SCREEN       (FROM VDWSCL=", pot%vdwscl, ")"
      do i = 1, pot%n_points
         write (unit, "(1X,A8,2F14.9)") pot%labels(i), 1.0_dp, pot%screen(i)
      end do
      write (unit, "(A)") "STOP"

      write (unit, "(A)") " $END"
      close (unit)

      if (present(omitted)) then
         omitted = "nothing"
      end if
   end subroutine write_efp_potential

   subroutine write_record(unit, label, values, width, decimals, per_line)
      !! A labelled record, continued with `>` and indented under its label
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: width, decimals, per_line

      character(len=MAX_LINE) :: line
      character(len=16) :: value_format
      integer :: i, first, last, at

      write (value_format, "(A,I0,A,I0,A)") "(F", width, ".", decimals, ")"
      first = 1
      do while (first <= size(values))
         last = min(first + per_line - 1, size(values))
         line = ""
         if (first == 1) then
            line(1:8) = label
         end if
         at = 9
         do i = first, last
            write (line(at:at + width - 1), value_format) values(i)
            at = at + width
         end do
         if (last < size(values)) line(at:at + 1) = " >"
         write (unit, "(A)") trim(line)
         first = last + 1
      end do
   end subroutine write_record

   subroutine write_tensor_point(unit, label, xyz, tensor, frequency)
      !! A point's coordinates on its label line, its nine components beneath
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: xyz(3)
      real(dp), intent(in) :: tensor(N_CART_PAIR)
      real(dp), intent(in), optional :: frequency

      if (present(frequency)) then
         ! "CT  1", five characters, then the coordinates -- not the label padded
         ! to eight the way the labelled record sections pad theirs.
         write (unit, "(A,3F15.10,A,F9.6,A)") trim(label), xyz, " -- FOR W=", &
            frequency, "I A.U."
      else
         ! trim, not a fixed width: this path serves both "CT1" from the static
         ! section and "CT  2" from a dynamic block's continuation points, and a
         ! fixed A3 would truncate the second.
         write (unit, "(A,3F15.10)") trim(label), xyz
      end if
      call write_values(unit, tensor, 16, 10, 4)
   end subroutine write_tensor_point

   subroutine write_values(unit, values, width, decimals, per_line)
      !! Unlabelled values, continued with `>`
      integer, intent(in) :: unit
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: width, decimals, per_line

      character(len=MAX_LINE) :: line
      character(len=16) :: value_format
      integer :: i, first, last, at

      write (value_format, "(A,I0,A,I0,A)") "(F", width, ".", decimals, ")"
      first = 1
      do while (first <= size(values))
         last = min(first + per_line - 1, size(values))
         line = ""
         at = 1
         do i = first, last
            write (line(at:at + width - 1), value_format) values(i)
            at = at + width
         end do
         if (last < size(values)) line(at:at + 1) = " >"
         write (unit, "(A)") trim(line)
         first = last + 1
      end do
   end subroutine write_values

   subroutine write_lower_triangle(unit, matrix)
      !! The lower triangle, row by row, four values to a line
      integer, intent(in) :: unit
      real(dp), intent(in) :: matrix(:, :)

      real(dp), allocatable :: packed(:)
      integer :: n, i, j, at

      n = size(matrix, 1)
      allocate (packed(n*(n + 1)/2))
      at = 0
      do i = 1, n
         do j = 1, i
            at = at + 1
            packed(at) = matrix(i, j)
         end do
      end do
      call write_values(unit, packed, 16, 10, 4)
      deallocate (packed)
   end subroutine write_lower_triangle

   subroutine write_wavefunction(unit, orbitals)
      !! `PROJECTION WAVEFUNCTION`: five coefficients a line, orbital by orbital
      !!
      !! The line carries the orbital index and a chunk counter, which is what
      !! GAMESS's reader keys on, so the chunk counter restarts with each orbital.
      integer, intent(in) :: unit
      real(dp), intent(in) :: orbitals(:, :)

      integer :: nao, n_lmo, i, start, last, chunk, k
      character(len=MAX_LINE) :: line
      integer :: at

      nao = size(orbitals, 1)
      n_lmo = size(orbitals, 2)
      do i = 1, n_lmo
         chunk = 0
         start = 1
         do while (start <= nao)
            last = min(start + 4, nao)
            chunk = chunk + 1
            line = ""
            write (line(1:5), "(I2,I3)") i, chunk
            at = 6
            do k = start, last
               write (line(at:at + 14), "(ES15.8E2)") orbitals(k, i)
               at = at + 15
            end do
            write (unit, "(A)") trim(line)
            start = last + 1
         end do
      end do
   end subroutine write_wavefunction

end module mqc_czt_efp_potential
