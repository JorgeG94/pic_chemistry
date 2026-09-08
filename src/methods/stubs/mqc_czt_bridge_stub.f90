!! Stand-in for the CPU SCF when the build has no libcint
module mqc_czt_bridge
   !! Same name and same entry points as the real bridge, declining.
   !!
   !! `czt_backend_available` is how a caller asks in advance, so the
   !! refusal can name the build option.
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   implicit none
   private

   public :: run_czt_hf
   public :: run_czt_mcscf
   public :: run_czt_fmo
   public :: run_czt_efmo
   public :: run_czt_makefp
   public :: run_czt_neo
   public :: run_czt_charges
   public :: run_czt_efp
   public :: run_czt_sapt0
   public :: run_czt_sapt2
   public :: czt_backend_available
   public :: xc_available
   public :: ecp_backend_available
   public :: czt_set_eri_path

contains

   subroutine czt_set_eri_path(name, error)
      !! No-op stand-in: with no CPU integral backend there is no path to choose
      use mqc_error, only: error_t
      character(len=*), intent(in) :: name
      type(error_t), intent(inout) :: error

      ! Nothing to choose and nothing to refuse; the arguments are the contract.
      associate (unused_name => name, unused_error => error)
      end associate
   end subroutine czt_set_eri_path

   pure function ecp_backend_available() result(available)
      !! No libcint means no integral backend at all, so no ECP either.
      logical :: available
      available = .false.
   end function ecp_backend_available

   pure function xc_available() result(available)
      !! No libcint means no `mqc_czt_xc`, and therefore no functionals --
      !! whatever libxc itself was configured to do.
      logical :: available
      available = .false.
   end function xc_available

   pure function czt_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .false.
   end function czt_backend_available

   subroutine run_czt_sapt0(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, basis_name, &
                            charge_a, charge_b, terms, error)
      !! No-op stand-in: SAPT0 needs the CPU integral backend
      use pic_types, only: dp
      use mqc_program_limits, only: N_SAPT_TERMS
      use mqc_error, only: error_t
      integer, intent(in) :: z_a(:), z_b(:)
      character(len=*), intent(in) :: sym_a(:), sym_b(:)
      real(dp), intent(in) :: xyz_a(:, :), xyz_b(:, :)
      character(len=*), intent(in) :: basis_name
      integer, intent(in) :: charge_a, charge_b
      real(dp), intent(out) :: terms(N_SAPT_TERMS)
      type(error_t), intent(inout) :: error

      terms = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "SAPT needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(z_a) < 0 .or. size(z_b) < 0) return
      if (len_trim(sym_a(1))*len_trim(sym_b(1))*len_trim(basis_name) < 0) return
      if (size(xyz_a) < 0 .or. size(xyz_b) < 0) return
      if (charge_a == huge(charge_a) .and. charge_b == huge(charge_b)) return
   end subroutine run_czt_sapt0

   subroutine run_czt_sapt2(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, basis_name, &
                            charge_a, charge_b, terms, error)
      !! No-op stand-in: SAPT2 needs the CPU integral backend
      use pic_types, only: dp
      use mqc_program_limits, only: N_SAPT2_TERMS
      use mqc_error, only: error_t
      integer, intent(in) :: z_a(:), z_b(:)
      character(len=*), intent(in) :: sym_a(:), sym_b(:)
      real(dp), intent(in) :: xyz_a(:, :), xyz_b(:, :)
      character(len=*), intent(in) :: basis_name
      integer, intent(in) :: charge_a, charge_b
      real(dp), intent(out) :: terms(N_SAPT2_TERMS)
      type(error_t), intent(inout) :: error

      terms = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "SAPT needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(z_a) < 0 .or. size(z_b) < 0) return
      if (len_trim(sym_a(1))*len_trim(sym_b(1))*len_trim(basis_name) < 0) return
      if (size(xyz_a) < 0 .or. size(xyz_b) < 0) return
      if (charge_a == huge(charge_a) .and. charge_b == huge(charge_b)) return
   end subroutine run_czt_sapt2

   subroutine run_czt_charges(atomic_numbers, element_symbols, coordinates, &
                              basis_name, scheme, total_charge, charges, error)
      !! No-op stand-in: atomic charges need the CPU integral backend
      use pic_types, only: dp
      use mqc_error, only: error_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: scheme
      integer, intent(in) :: total_charge
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      allocate (charges(0))
      call error%set(ERROR_VALIDATION, &
                     "atomic charges need the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0) return
      if (len_trim(element_symbols(1))*len_trim(basis_name)*len_trim(scheme) < 0) return
      if (total_charge < -huge(1)) return
   end subroutine run_czt_charges

   subroutine run_czt_efp(potentials, fragment_sizes, fragment_atoms, &
                          coordinates, terms, error)
      !! No-op stand-in: an EFP interaction energy needs the CPU backend
      use pic_types, only: dp
      use mqc_program_limits, only: N_EFP_TERMS
      use mqc_error, only: error_t
      character(len=*), intent(in) :: potentials(:)
      integer, intent(in) :: fragment_sizes(:)
      integer, intent(in) :: fragment_atoms(:, :)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: terms(N_EFP_TERMS)
      type(error_t), intent(inout) :: error

      terms = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "EFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (len_trim(potentials(1)) < 0) return
      if (size(fragment_sizes) < 0 .or. size(fragment_atoms) < 0) return
      if (size(coordinates) < 0) return
   end subroutine run_czt_efp
   subroutine run_czt_fmo(atomic_numbers, element_symbols, coordinates, owner, &
                          basis_name, esp, expansion, far_field, resppc, &
                          level, max_outer, outer_tol, scf_max_iter, &
                          scf_energy_tol, scf_density_tol, scf_drive, &
                          bond_breaking, &
                          cap_scale, energy, error, comm)
      !! No-op stand-in: FMO needs the CPU integral backend
      !!
      !! Coordinates are Bohr; `owner(i)` is atom i's fragment, numbered from
      !! one with no gaps.
      use pic_types, only: dp
      use mqc_error, only: error_t
      use pic_mpi_lib, only: comm_t
      use mqc_scf_types, only: scf_numerics_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: owner(:)
      character(len=*), intent(in) :: basis_name, esp, expansion, far_field
      real(dp), intent(in) :: resppc
      integer, intent(in) :: level
      integer, intent(in) :: max_outer
      real(dp), intent(in) :: outer_tol
      integer, intent(in) :: scf_max_iter
      real(dp), intent(in) :: scf_energy_tol, scf_density_tol
      type(scf_numerics_t), intent(in) :: scf_drive
      character(len=*), intent(in) :: bond_breaking
      real(dp), intent(in) :: cap_scale
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm
         !! Present means distribute the fragment work over this communicator.
         !! Absent means one rank does all of it.

      energy = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "FMO needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0 .or. size(owner) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(esp)*len_trim(expansion)*len_trim(far_field) < 0) return
      if (resppc < -huge(1.0_dp) .or. level < 0 .or. max_outer < 0) return
      if (outer_tol < 0.0_dp .or. scf_max_iter < 0 .or. scf_energy_tol < 0.0_dp) return
      if (scf_density_tol < 0.0_dp) return
      if (present(comm)) return
   end subroutine run_czt_fmo

   subroutine run_czt_efmo(atomic_numbers, element_symbols, coordinates, owner, &
                           fragment_charges, basis_name, rcut, charge_transfer, &
                           scf_drive, scf_max_iter, scf_energy_tol, scf_density_tol, &
                           scf_grad_tol, guess, energy, terms, n_qm_pairs, n_efp_pairs, &
                           error, verbose, aux_basis, vdwscl, quadrupole_blocks, &
                           dynamic_tol, dynamic_maxiter, response, &
                           allow_crap_response, response_batch)
      !! No-op stand-in: EFMO needs the CPU integral backend
      !!
      !! Coordinates are Bohr; `owner(i)` is atom i's fragment, numbered from
      !! one with no gaps.
      use pic_types, only: dp
      use mqc_error, only: error_t
      use mqc_scf_types, only: scf_numerics_t
      use mqc_program_limits, only: N_EFMO_TERMS
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: owner(:)
      integer, intent(in) :: fragment_charges(:)
      character(len=*), intent(in) :: basis_name
      real(dp), intent(in) :: rcut
      logical, intent(in) :: charge_transfer
      type(scf_numerics_t), intent(in) :: scf_drive
      integer, intent(in) :: scf_max_iter
      real(dp), intent(in) :: scf_energy_tol, scf_density_tol, scf_grad_tol
      character(len=*), intent(in) :: guess
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: terms(N_EFMO_TERMS)
      integer, intent(out) :: n_qm_pairs, n_efp_pairs
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: aux_basis
      real(dp), intent(in), optional :: vdwscl
      logical, intent(in), optional :: quadrupole_blocks
      real(dp), intent(in), optional :: dynamic_tol
      integer, intent(in), optional :: dynamic_maxiter
      integer, intent(in), optional :: response
      logical, intent(in), optional :: allow_crap_response
      integer, intent(in), optional :: response_batch

      energy = 0.0_dp
      terms = 0.0_dp
      n_qm_pairs = 0
      n_efp_pairs = 0
      call error%set(ERROR_VALIDATION, &
                     "EFMO needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0 .or. size(owner) < 0) return
      if (size(fragment_charges) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(guess) < 0) return
      if (rcut < -huge(1.0_dp) .or. charge_transfer) return
      if (scf_drive%max_iter < 0 .or. scf_max_iter < 0) return
      if (scf_energy_tol < 0.0_dp .or. scf_density_tol < 0.0_dp) return
      if (scf_grad_tol < 0.0_dp) return
      if (present(verbose) .or. present(aux_basis) .or. present(vdwscl)) return
      if (present(quadrupole_blocks) .or. present(dynamic_tol)) return
      if (present(dynamic_maxiter) .or. present(response)) return
      if (present(allow_crap_response) .or. present(response_batch)) return
   end subroutine run_czt_efmo

   subroutine run_czt_makefp(atomic_numbers, element_symbols, coordinates, &
                             basis_name, name, path, error, charge, verbose, &
                             aux_basis, guess, energy_tol, density_tol, grad_tol, &
                             scf_in, max_iter_in, &
                             vdwscl, dynamic_tol, dynamic_maxiter, response, &
                             allow_crap_response, response_batch, quadrupole_blocks)
      !! No-op stand-in: an effective fragment potential needs the CPU backend
      use pic_types, only: dp
      use mqc_error, only: error_t
      use mqc_scf_types, only: scf_numerics_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name, name, path
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge   !! Net charge; a fragment may be an ion
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: aux_basis
      character(len=*), intent(in), optional :: guess
      real(dp), intent(in), optional :: energy_tol, density_tol, grad_tol
      type(scf_numerics_t), intent(in), optional :: scf_in
      integer, intent(in), optional :: max_iter_in
      real(dp), intent(in), optional :: vdwscl, dynamic_tol
      integer, intent(in), optional :: dynamic_maxiter, response
      logical, intent(in), optional :: allow_crap_response
      logical, intent(in), optional :: quadrupole_blocks
      integer, intent(in), optional :: response_batch

      call error%set(ERROR_VALIDATION, &
                     "MAKEFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(name)*len_trim(path) < 0) return
      if (present(charge) .or. present(verbose)) return
      if (present(aux_basis) .or. present(guess)) return
      if (present(energy_tol) .or. present(density_tol)) return
      if (present(vdwscl) .or. present(dynamic_tol)) return
      if (present(dynamic_maxiter) .or. present(response)) return
      if (present(allow_crap_response) .or. present(response_batch)) return
   end subroutine run_czt_makefp

   subroutine run_czt_neo(atomic_numbers, element_symbols, coordinates, basis_name, &
                          nuclear_basis, quantum, charge, energy, error, verbose, &
                          energy_tol, density_tol, max_iter, functional, grid_level, epc, &
                          scf_in)
      !! No-op stand-in: quantum nuclei need the CPU backend
      use pic_types, only: dp
      use mqc_error, only: error_t
      use mqc_scf_types, only: scf_numerics_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name, nuclear_basis
      logical, intent(in) :: quantum(:)
      integer, intent(in) :: charge
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: energy_tol, density_tol
      integer, intent(in), optional :: max_iter
      character(len=*), intent(in), optional :: functional, epc
      integer, intent(in), optional :: grid_level
      type(scf_numerics_t), intent(in), optional :: scf_in
      energy = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "keywords.neo needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_CZT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0 .or. size(quantum) < 0) return
      if (size(element_symbols) < 0 .or. charge < -huge(charge)) return
      if (len_trim(basis_name)*len_trim(nuclear_basis) < 0) return
      if (present(verbose) .or. present(energy_tol)) return
      if (present(density_tol) .or. present(max_iter)) return
      if (present(functional) .or. present(grid_level) .or. present(epc)) return
      if (present(scf_in)) return
   end subroutine run_czt_neo

   subroutine run_czt_hf(settings, fragment, result, want_gradient, want_hessian)
      !! No-op stand-in: report the missing backend, compute nothing
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient
      logical, intent(in), optional :: want_hessian

      call result%error%set(ERROR_VALIDATION, &
                            "This calculation needs an integral backend; build with "// &
                            "-DMQC_ENABLE_CZT=ON for the CPU one, or "// &
                            "-DMQC_ENABLE_CUEST=ON for the GPU one")
      result%has_error = .true.
      result%has_energy = .false.
      if (len_trim(settings%basis_set) == 0 .or. fragment%n_atoms < 0) return
      if (present(want_gradient)) return
      if (present(want_hessian)) return
   end subroutine run_czt_hf

   subroutine run_czt_mcscf(settings, fragment, result, want_gradient)
      !! No-op stand-in: CASSCF and CASCI need the CPU integral backend
      !!
      !! There is no GPU path to fall through to.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      call result%error%set(ERROR_VALIDATION, &
                            "a multiconfigurational calculation needs the CPU integral "// &
                            "backend; build with -DMQC_ENABLE_CZT=ON")
      result%has_error = .true.
      result%has_energy = .false.
      if (len_trim(settings%basis_set) == 0 .or. fragment%n_atoms < 0) return
   end subroutine run_czt_mcscf

end module mqc_czt_bridge
