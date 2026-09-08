!! Centralized JSON output writer
!! Single entry point for all JSON output in metalquicha
module mqc_json_writer
   use pic_types, only: int64, dp
   use pic_logger, only: logger => global_logger
   use mqc_json_output_types, only: json_output_data_t, &
                                    OUTPUT_MODE_UNFRAGMENTED, OUTPUT_MODE_MBE, OUTPUT_MODE_GMBE_PIE
   use mqc_io_helpers, only: get_output_json_filename, get_basename
   use mqc_physical_constants, only: HARTREE_TO_EV
   use mqc_physical_constants, only: HARTREE_TO_CALMOL, R_CALMOLK, AU_TO_DEBYE, CAL_TO_J
   use mqc_program_limits, only: JSON_REAL_FORMAT
   use mqc_mbe_io, only: get_frag_level_name
   use mqc_fragment_table_writer, only: write_fragment_table
   use json_module, only: json_core, json_value
   implicit none
   private

   public :: write_json_output  !! Single entry point for all JSON output

contains

   subroutine add_gradient(json, parent, gradient)
      !! Nuclear gradient, as the norm and as the components it came from
      !!
      !! The norm alone cannot separate a correct gradient from one with a sign
      !! error, a swapped pair of atoms, or a wrong component that preserves the
      !! magnitude, so validation compares the components. Hartree/Bohr, the
      !! internal unit -- unlike the optimization trajectory, which converts to
      !! Angstrom.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      real(dp), intent(in) :: gradient(:, :)  !! (3, natoms)

      type(json_value), pointer :: grad_arr, atom_arr
      integer :: iatom, icomp

      call json%add(parent, "gradient_norm", sqrt(sum(gradient**2)))
      call json%add(parent, "gradient_units", "hartree/bohr")

      ! One [x, y, z] triple per atom, in input order.
      call json%create_array(grad_arr, "gradient")
      call json%add(parent, grad_arr)
      do iatom = 1, size(gradient, 2)
         call json%create_array(atom_arr, "")
         call json%add(grad_arr, atom_arr)
         do icomp = 1, size(gradient, 1)
            call json%add(atom_arr, "", gradient(icomp, iatom))
         end do
      end do

   end subroutine add_gradient

   subroutine add_hessian(json, parent, hessian)
      !! Second derivatives, as the Frobenius norm and as the matrix itself
      !!
      !! The norm is one number out of `9N^2` and invariant under any orthogonal
      !! mixing, so a transposed block, a swapped atom pair or a sign flip in an
      !! off-diagonal element all leave it unchanged. It is kept for a manifest
      !! of many systems; it is not what a Hessian should be validated on.
      !!
      !! One row per array row, in the `(3N, 3N)` order the vibrational
      !! analysis reads, atom-major with `x, y, z` inside each atom. Hartree
      !! per Bohr squared, the internal unit.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      real(dp), intent(in) :: hessian(:, :)  !! (3*natoms, 3*natoms)

      type(json_value), pointer :: hess_arr, row_arr
      integer :: i, j

      call json%add(parent, "hessian_frobenius_norm", sqrt(sum(hessian**2)))
      call json%add(parent, "hessian_units", "hartree/bohr^2")

      call json%create_array(hess_arr, "hessian")
      call json%add(parent, hess_arr)
      do i = 1, size(hessian, 1)
         call json%create_array(row_arr, "")
         call json%add(hess_arr, row_arr)
         do j = 1, size(hessian, 2)
            call json%add(row_arr, "", hessian(i, j))
         end do
      end do

   end subroutine add_hessian

   subroutine write_json_output(output_data)
      !! The single entry point for all JSON output
      !!
      !! Dispatches on `output_mode`. The only place in the codebase where JSON
      !! files are written.
      type(json_output_data_t), intent(in) :: output_data

      select case (output_data%output_mode)
      case (OUTPUT_MODE_UNFRAGMENTED)
         if (output_data%has_vibrational) then
            call write_vibrational_json_impl(output_data)
         else
            call write_unfragmented_json_impl(output_data)
         end if

      case (OUTPUT_MODE_MBE)
         if (output_data%has_vibrational) then
            call write_vibrational_json_impl(output_data)
         else
            ! The per-fragment table goes to exactly one place: in the JSON it
            ! costs six times as much per fragment and forces the consumer to
            ! hold the whole document resident. The summary is written either
            ! way, so "csv" loses nothing but the duplication.
            call write_mbe_breakdown_json_impl(output_data)
            if (trim(output_data%fragment_breakdown) == "csv") then
               call write_fragment_table(output_data)
            end if
         end if

      case (OUTPUT_MODE_GMBE_PIE)
         if (output_data%has_vibrational) then
            call write_vibrational_json_impl(output_data)
         else
            call write_gmbe_pie_json_impl(output_data)
         end if

      case default
         call logger%error("Unknown output mode in write_json_output")
      end select

   end subroutine write_json_output

   subroutine write_unfragmented_json_impl(data)
      !! Write unfragmented calculation results to output JSON file
      type(json_output_data_t), intent(in) :: data

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, dipole_obj
      integer :: iunit, io_stat
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize(real_format=JSON_REAL_FORMAT)
      call json%create_object(root, "")
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)
      ! Stamped by every writer; a restart reads it before it reuses anything.
      if (len_trim(data%fingerprint) > 0) then
         call json%add(main_obj, "fingerprint", trim(data%fingerprint))
      end if

      if (data%has_energy) call json%add(main_obj, "total_energy", data%total_energy)

      call write_sapt_section(json, main_obj, data)
      call write_efmo_section(json, main_obj, data)
      call write_ieda_section(json, main_obj, data)
      call write_charges_section(json, main_obj, data)
      call write_fukui_section(json, main_obj, data)

      ! Only where one SCF covered one system. A fragmented run never sets
      ! this, because a gap assembled from fragment gaps would be arithmetic
      ! without a meaning.
      if (data%has_orbitals) then
         call json%add(main_obj, "homo", data%homo)
         call json%add(main_obj, "lumo", data%lumo)
         call json%add(main_obj, "homo_lumo_gap_ev", (data%lumo - data%homo)*HARTREE_TO_EV)
      end if

      if (data%has_dipole .and. allocated(data%dipole)) then
         call json%create_object(dipole_obj, "dipole")
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, "x", data%dipole(1))
         call json%add(dipole_obj, "y", data%dipole(2))
         call json%add(dipole_obj, "z", data%dipole(3))
         call json%add(dipole_obj, "magnitude_debye", norm2(data%dipole)*AU_TO_DEBYE)
      end if

      if (data%has_gradient .and. allocated(data%gradient)) then
         call add_gradient(json, main_obj, data%gradient)
      end if
      if (data%has_hessian .and. allocated(data%hessian)) then
         call add_hessian(json, main_obj, data%hessian)
      end if

      call logger%info("Writing JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status="replace", action="write", iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("JSON output written successfully to "//trim(output_file))

   end subroutine write_unfragmented_json_impl

   subroutine write_mbe_breakdown_json_impl(data)
      !! Write detailed MBE energy breakdown to JSON file
      type(json_output_data_t), intent(in) :: data

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, levels_arr, level_obj, frags_arr, frag_obj
      type(json_value), pointer :: dipole_obj
      integer(int64) :: i, count_by_level
      ! TODO(mqc): `j` is dead here.
      integer :: fragment_size, j, frag_level, iunit, io_stat
      integer, allocatable :: indices(:)
      character(len=32) :: level_name
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      if (data%max_level > 10) then
         call logger%warning("Fragment levels exceed decamers (10-mers). JSON will use generic N-mers notation.")
      end if

      call json%initialize(real_format=JSON_REAL_FORMAT)
      call json%create_object(root, "")
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)
      ! Stamped by every writer; a restart reads it before it reuses anything.
      if (len_trim(data%fingerprint) > 0) then
         call json%add(main_obj, "fingerprint", trim(data%fingerprint))
      end if

      call json%add(main_obj, "total_energy", data%total_energy)

      call write_unconverged_section(json, main_obj, data)

      ! Build levels array
      call json%create_array(levels_arr, "levels")
      call json%add(main_obj, levels_arr)

      do frag_level = 1, data%max_level
         count_by_level = 0_int64
         do i = 1_int64, data%fragment_count
            fragment_size = count(data%polymers(i, :) > 0)
            if (fragment_size == frag_level) count_by_level = count_by_level + 1_int64
         end do

         if (count_by_level > 0_int64) then
            call json%create_object(level_obj, "")
            call json%add(levels_arr, level_obj)

            level_name = get_frag_level_name(frag_level)
            call json%add(level_obj, "frag_level", frag_level)
            call json%add(level_obj, "name", trim(level_name))
            call json%add(level_obj, "count", int(count_by_level))
            if (allocated(data%sum_by_level)) then
               call json%add(level_obj, "total_energy", data%sum_by_level(frag_level))
            end if

            ! Per-fragment detail only when this is the chosen sink for it
            if (trim(data%fragment_breakdown) /= "json") cycle

            call json%create_array(frags_arr, "fragments")
            call json%add(level_obj, frags_arr)

            do i = 1_int64, data%fragment_count
               fragment_size = count(data%polymers(i, :) > 0)
               if (fragment_size == frag_level) then
                  call json%create_object(frag_obj, "")
                  call json%add(frags_arr, frag_obj)

                  allocate (indices(fragment_size))
                  indices = data%polymers(i, 1:fragment_size)
                  call json%add(frag_obj, "indices", indices)
                  deallocate (indices)

                  if (allocated(data%fragment_energies)) then
                     call json%add(frag_obj, "energy", data%fragment_energies(i))
                  end if

                  if (allocated(data%fragment_distances)) then
                     call json%add(frag_obj, "distance", data%fragment_distances(i))
                  end if

                  if (frag_level > 1 .and. allocated(data%delta_energies)) then
                     call json%add(frag_obj, "delta_energy", data%delta_energies(i))
                  end if
               end if
            end do
         end if
      end do

      ! Add dipole if computed
      if (data%has_dipole .and. allocated(data%dipole)) then
         call json%create_object(dipole_obj, "dipole")
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, "x", data%dipole(1))
         call json%add(dipole_obj, "y", data%dipole(2))
         call json%add(dipole_obj, "z", data%dipole(3))
         call json%add(dipole_obj, "magnitude_debye", norm2(data%dipole)*AU_TO_DEBYE)
      end if

      if (data%has_gradient .and. allocated(data%gradient)) then
         call add_gradient(json, main_obj, data%gradient)
      end if

      if (data%has_hessian .and. allocated(data%hessian)) then
         call add_hessian(json, main_obj, data%hessian)
      end if

      ! Write to file
      call logger%info("Writing JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status="replace", action="write", iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("JSON output written successfully to "//trim(output_file))

   end subroutine write_mbe_breakdown_json_impl

   subroutine write_charges_section(json, parent, data)
      !! Atomic partial charges, per atom, with the scheme that produced them
      !!
      !! The scheme travels with the numbers: two schemes disagree by design, so
      !! a charge without its scheme is not comparable against anything.
      !!
      !! `spin_populations` appears only for an unrestricted reference under
      !! Mulliken. Its absence is meaningful and not a gap: a closed shell has
      !! no spin density, and CHELPG fits a potential the spin density does not
      !! produce.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      type(json_output_data_t), intent(in) :: data

      type(json_value), pointer :: section, arr, entry
      integer :: i, natm
      logical :: with_spin

      if (.not. data%has_charges) return
      if (.not. allocated(data%atomic_charges)) return
      natm = size(data%atomic_charges)
      with_spin = allocated(data%spin_populations)
      if (with_spin) with_spin = size(data%spin_populations) == natm

      call json%create_object(section, "atomic_charges")
      call json%add(parent, section)
      call json%add(section, "scheme", trim(data%charge_scheme))
      ! What the column has to add up to, so a consumer can check rather than
      ! trust. Mulliken satisfies it as a trace identity and CHELPG as the
      ! constraint its fit is solved under, so a mismatch is a bookkeeping bug
      ! and never a property of the molecule.
      call json%add(section, "sum", sum(data%atomic_charges))

      call json%create_array(arr, "atoms")
      call json%add(section, arr)
      do i = 1, natm
         call json%create_object(entry, "")
         call json%add(arr, entry)
         ! Numbered rather than named, as everywhere else here: the payload
         ! carries no element symbols and the order is the input order.
         call json%add(entry, "atom", i)
         call json%add(entry, "charge", data%atomic_charges(i))
         if (with_spin) then
            call json%add(entry, "spin_population", data%spin_populations(i))
         end if
      end do
   end subroutine write_charges_section

   subroutine write_fukui_section(json, parent, data)
      !! Where the molecule reacts, per atom, for something other than a reader
      !!
      !! Emitted per atom rather than reduced to a most-reactive site: which
      !! index to rank on depends on the reaction being asked about -- `f+` for
      !! an incoming nucleophile, `f-` for an electrophile, the dual descriptor
      !! to tell the two kinds of site apart.
      !!
      !! `anion_bound` travels with the numbers. When it is false the `f_plus`
      !! column describes an orbital the basis invented rather than a real
      !! anion, and nothing about the values themselves says so -- a consumer
      !! sorting on `f_plus` would rank sites confidently off a fiction.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      type(json_output_data_t), intent(in) :: data

      type(json_value), pointer :: section, arr, entry
      integer :: i, natm

      if (.not. data%has_fukui) return
      if (.not. allocated(data%fukui_plus)) return
      natm = size(data%fukui_plus)

      call json%create_object(section, "fukui")
      call json%add(parent, section)
      call json%add(section, "population_scheme", trim(data%fukui_scheme))
      call json%add(section, "ionisation_potential", data%fukui_ip)
      call json%add(section, "electron_affinity", data%fukui_ea)
      call json%add(section, "chemical_potential", -0.5_dp*(data%fukui_ip + data%fukui_ea))
      call json%add(section, "hardness", data%fukui_hardness)
      call json%add(section, "electrophilicity", data%fukui_electrophilicity)
      call json%add(section, "anion_bound", data%fukui_anion_bound)

      call json%create_array(arr, "atoms")
      call json%add(section, arr)
      do i = 1, natm
         call json%create_object(entry, "")
         call json%add(arr, entry)
         ! Numbered rather than named: the payload carries no element symbols
         ! and the order is the input order.
         call json%add(entry, "atom", i)
         call json%add(entry, "f_plus", data%fukui_plus(i))
         call json%add(entry, "f_minus", data%fukui_minus(i))
         call json%add(entry, "f_zero", &
                       0.5_dp*(data%fukui_plus(i) + data%fukui_minus(i)))
         call json%add(entry, "dual", data%fukui_dual(i))
      end do
   end subroutine write_fukui_section

   subroutine write_unconverged_section(json, parent, data)
      !! The fragments whose SCF failed, in a form a follow-up job can be built from
      !!
      !! Written whenever the method reports convergence at all, **including
      !! when nothing failed**: a `count` of zero is a statement, so a consumer
      !! reads a missing section as "no information" rather than "all good".
      !!
      !! Each entry carries the monomers the fragment was built from, so a
      !! dimer can be reconstructed without reading back the per-fragment
      !! table -- which for a large run is a separate CSV of millions of rows.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      type(json_output_data_t), intent(in) :: data

      type(json_value), pointer :: section, arr, entry
      integer, allocatable :: members(:)
      integer(int64) :: i
      integer :: level

      if (.not. allocated(data%unconverged_ids)) return

      call json%create_object(section, "unconverged")
      call json%add(parent, section)
      call json%add(section, "count", int(size(data%unconverged_ids)))

      call json%create_array(arr, "fragments")
      call json%add(section, arr)
      do i = 1_int64, int(size(data%unconverged_ids), int64)
         call json%create_object(entry, "")
         call json%add(arr, entry)
         call json%add(entry, "id", int(data%unconverged_ids(i)))
         ! Zero-padded in storage, so the real membership is the non-zero part
         ! and the level is how many of those there are.
         members = pack(data%unconverged_monomers(i, :), &
                        data%unconverged_monomers(i, :) > 0)
         level = size(members)
         call json%add(entry, "level", level)
         call json%add(entry, "monomers", members)
         if (allocated(data%unconverged_deltas)) then
            call json%add(entry, "delta_energy", data%unconverged_deltas(i))
         end if
      end do

      ! How much of the total rests on fragments that never converged.
      if (size(data%unconverged_ids) > 0 .and. allocated(data%unconverged_deltas)) then
         ! Signed, because the failures are a subset of a sum with cancellation
         ! in it. The largest single one is given beside it, since a net near
         ! zero can still be two large terms that happen to oppose.
         call json%add(section, "energy_at_stake", sum(data%unconverged_deltas))
         call json%add(section, "largest_contribution", &
                       maxval(abs(data%unconverged_deltas)))
      end if

      ! Failures cluster on their cause, so the monomers that keep appearing are
      ! the short list worth looking at. Four hundred failed dimers sharing one
      ! monomer is one bad monomer.
      if (allocated(data%culprit_monomers)) then
         if (size(data%culprit_monomers) > 0) then
            call json%create_array(arr, "monomers_involved")
            call json%add(section, arr)
            do i = 1_int64, int(size(data%culprit_monomers), int64)
               call json%create_object(entry, "")
               call json%add(arr, entry)
               call json%add(entry, "monomer", data%culprit_monomers(i))
               call json%add(entry, "fragments", int(data%culprit_counts(i)))
            end do
         end if
      end if
   end subroutine write_unconverged_section

   subroutine write_gmbe_pie_json_impl(data)
      !! Write GMBE PIE calculation results to output JSON file
      type(json_output_data_t), intent(in) :: data

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, pie_obj, terms_arr, term_obj
      ! TODO(mqc): `j` is dead here.
      integer :: j, max_atoms, n_atoms, iunit, io_stat
      integer(int64) :: i, n_nonzero_terms
      integer, allocatable :: atom_indices(:)
      character(len=256) :: output_file, basename

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize(real_format=JSON_REAL_FORMAT)
      call json%create_object(root, "")
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)
      ! Stamped by every writer; a restart reads it before it reuses anything.
      if (len_trim(data%fingerprint) > 0) then
         call json%add(main_obj, "fingerprint", trim(data%fingerprint))
      end if

      call json%add(main_obj, "total_energy", data%total_energy)

      if (data%has_gradient .and. allocated(data%gradient)) then
         call add_gradient(json, main_obj, data%gradient)
      end if
      if (data%has_hessian .and. allocated(data%hessian)) then
         call add_hessian(json, main_obj, data%hessian)
      end if

      ! Count non-zero coefficient terms
      if (allocated(data%pie_coefficients)) then
         n_nonzero_terms = count(data%pie_coefficients(1:data%n_pie_terms) /= 0)
      else
         n_nonzero_terms = 0
      end if

      ! PIE terms section
      call json%create_object(pie_obj, "pie_terms")
      call json%add(main_obj, pie_obj)
      call json%add(pie_obj, "count", int(n_nonzero_terms))

      call json%create_array(terms_arr, "terms")
      call json%add(pie_obj, terms_arr)

      if (allocated(data%pie_atom_sets) .and. allocated(data%pie_coefficients) .and. &
          allocated(data%pie_energies)) then
         max_atoms = size(data%pie_atom_sets, 1)

         do i = 1_int64, data%n_pie_terms
            if (data%pie_coefficients(i) == 0) cycle

            call json%create_object(term_obj, "")
            call json%add(terms_arr, term_obj)

            ! Atoms until the negative sentinel. The bound and the sentinel are
            ! tested in separate statements on purpose: Fortran does not promise
            ! to short-circuit `.and.`, so writing both in one condition lets the
            ! compiler read `pie_atom_sets(max_atoms + 1, i)` on the last pass of
            ! a full set -- one element past the column, which `-fcheck=bounds`
            ! traps and a release build reads out of the next column in silence.
            n_atoms = 0
            do while (n_atoms < max_atoms)
               if (data%pie_atom_sets(n_atoms + 1, i) < 0) exit
               n_atoms = n_atoms + 1
            end do

            allocate (atom_indices(n_atoms))
            atom_indices = data%pie_atom_sets(1:n_atoms, i)
            call json%add(term_obj, "atom_indices", atom_indices)
            deallocate (atom_indices)

            call json%add(term_obj, "coefficient", data%pie_coefficients(i))
            call json%add(term_obj, "energy", data%pie_energies(i))
            call json%add(term_obj, "weighted_energy", real(data%pie_coefficients(i), dp)*data%pie_energies(i))
         end do
      end if

      ! Write to file
      call logger%info("Writing GMBE PIE JSON output to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status="replace", action="write", iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)
      call json%destroy(root)
      call logger%info("GMBE PIE JSON output written successfully to "//trim(output_file))

   end subroutine write_gmbe_pie_json_impl

   subroutine write_vibrational_json_impl(data)
      !! Write vibrational analysis and thermochemistry results to JSON file
      type(json_output_data_t), intent(in) :: data

      type(json_core) :: json
      type(json_value), pointer :: root, main_obj, dipole_obj, vib_obj, thermo_obj
      type(json_value), pointer :: moi_obj, rot_obj, pf_obj, contrib_obj, table_obj
      type(json_value), pointer :: trans_obj, rot_contrib, vib_contrib, elec_obj
      type(json_value), pointer :: vib_row, rot_row, int_row, tr_row, tot_row
      type(json_value), pointer :: thermal_obj, total_e_obj
      integer :: io_stat, iunit
      character(len=256) :: output_file, basename
      real(dp) :: pV_cal, H_vib_cal, H_rot_cal, H_trans_cal, H_total_cal
      real(dp) :: Cv_total, S_total, S_total_J, H_int_cal, Cv_int, S_int, S_int_J, Cp_trans
      ! TODO(mqc): `hess_norm` is computed below and never read; `add_hessian`
      ! writes the norm itself.
      real(dp) :: hess_norm

      output_file = get_output_json_filename()
      basename = get_basename()

      call json%initialize(real_format=JSON_REAL_FORMAT)

      ! Create root object
      call json%create_object(root, "")

      ! Create main object with basename as key
      call json%create_object(main_obj, trim(basename))
      call json%add(root, main_obj)
      ! Stamped by every writer; a restart reads it before it reuses anything.
      if (len_trim(data%fingerprint) > 0) then
         call json%add(main_obj, "fingerprint", trim(data%fingerprint))
      end if

      ! Total energy
      if (data%has_energy) call json%add(main_obj, "total_energy", data%total_energy)

      ! Only where one SCF covered one system. A fragmented run never sets
      ! this, because a gap assembled from fragment gaps would be arithmetic
      ! without a meaning.
      if (data%has_orbitals) then
         call json%add(main_obj, "homo", data%homo)
         call json%add(main_obj, "lumo", data%lumo)
         call json%add(main_obj, "homo_lumo_gap_ev", (data%lumo - data%homo)*HARTREE_TO_EV)
      end if

      ! Dipole
      if (data%has_dipole .and. allocated(data%dipole)) then
         call json%create_object(dipole_obj, "dipole")
         call json%add(main_obj, dipole_obj)
         call json%add(dipole_obj, "x", data%dipole(1))
         call json%add(dipole_obj, "y", data%dipole(2))
         call json%add(dipole_obj, "z", data%dipole(3))
         call json%add(dipole_obj, "magnitude_debye", norm2(data%dipole)*AU_TO_DEBYE)
      end if

      ! Gradient and Hessian norms
      if (data%has_gradient .and. allocated(data%gradient)) then
         call add_gradient(json, main_obj, data%gradient)
      end if
      if (data%has_hessian .and. allocated(data%hessian)) then
         hess_norm = sqrt(sum(data%hessian**2))
         call add_hessian(json, main_obj, data%hessian)
      end if

      ! Vibrational analysis section
      call json%create_object(vib_obj, "vibrational_analysis")
      call json%add(main_obj, vib_obj)
      if (allocated(data%frequencies)) then
         call json%add(vib_obj, "n_modes", size(data%frequencies))
         call json%add(vib_obj, "frequencies_cm1", data%frequencies)
      end if
      if (allocated(data%reduced_masses)) then
         call json%add(vib_obj, "reduced_masses_amu", data%reduced_masses)
      end if
      if (allocated(data%force_constants)) then
         call json%add(vib_obj, "force_constants_mdyne_ang", data%force_constants)
      end if
      if (data%has_ir_intensities .and. allocated(data%ir_intensities)) then
         call json%add(vib_obj, "ir_intensities_km_mol", data%ir_intensities)
      end if

      ! Thermochemistry section
      call json%create_object(thermo_obj, "thermochemistry")
      call json%add(main_obj, thermo_obj)

      ! Conditions
      call json%add(thermo_obj, "temperature_K", data%thermo%temperature)
      call json%add(thermo_obj, "pressure_atm", data%thermo%pressure)
      call json%add(thermo_obj, "molecular_mass_amu", data%thermo%total_mass)
      call json%add(thermo_obj, "symmetry_number", data%thermo%symmetry_number)
      call json%add(thermo_obj, "spin_multiplicity", data%thermo%spin_multiplicity)
      call json%add(thermo_obj, "is_linear", data%thermo%is_linear)
      call json%add(thermo_obj, "n_real_frequencies", data%thermo%n_real_freqs)
      call json%add(thermo_obj, "n_imaginary_frequencies", data%thermo%n_imag_freqs)

      ! Moments of inertia
      call json%create_object(moi_obj, "moments_of_inertia_amu_ang2")
      call json%add(thermo_obj, moi_obj)
      call json%add(moi_obj, "Ia", data%thermo%moments(1))
      call json%add(moi_obj, "Ib", data%thermo%moments(2))
      call json%add(moi_obj, "Ic", data%thermo%moments(3))

      ! Rotational constants
      call json%create_object(rot_obj, "rotational_constants_GHz")
      call json%add(thermo_obj, rot_obj)
      call json%add(rot_obj, "A", data%thermo%rot_const(1))
      call json%add(rot_obj, "B", data%thermo%rot_const(2))
      call json%add(rot_obj, "C", data%thermo%rot_const(3))

      ! Partition functions
      call json%create_object(pf_obj, "partition_functions")
      call json%add(thermo_obj, pf_obj)
      call json%add(pf_obj, "translational", data%thermo%q_trans)
      call json%add(pf_obj, "rotational", data%thermo%q_rot)
      call json%add(pf_obj, "vibrational", data%thermo%q_vib)

      ! Thermodynamic contributions
      call json%create_object(contrib_obj, "contributions")
      call json%add(thermo_obj, contrib_obj)

      call json%create_object(trans_obj, "translational")
      call json%add(contrib_obj, trans_obj)
      call json%add(trans_obj, "energy_hartree", data%thermo%E_trans)
      call json%add(trans_obj, "entropy_cal_mol_K", data%thermo%S_trans)
      call json%add(trans_obj, "Cv_cal_mol_K", data%thermo%Cv_trans)

      call json%create_object(rot_contrib, "rotational")
      call json%add(contrib_obj, rot_contrib)
      call json%add(rot_contrib, "energy_hartree", data%thermo%E_rot)
      call json%add(rot_contrib, "entropy_cal_mol_K", data%thermo%S_rot)
      call json%add(rot_contrib, "Cv_cal_mol_K", data%thermo%Cv_rot)

      call json%create_object(vib_contrib, "vibrational")
      call json%add(contrib_obj, vib_contrib)
      call json%add(vib_contrib, "energy_hartree", data%thermo%E_vib)
      call json%add(vib_contrib, "entropy_cal_mol_K", data%thermo%S_vib)
      call json%add(vib_contrib, "Cv_cal_mol_K", data%thermo%Cv_vib)

      call json%create_object(elec_obj, "electronic")
      call json%add(contrib_obj, elec_obj)
      call json%add(elec_obj, "energy_hartree", data%thermo%E_elec)
      call json%add(elec_obj, "entropy_cal_mol_K", data%thermo%S_elec)

      ! Contribution table
      pV_cal = R_CALMOLK*data%thermo%temperature
      H_vib_cal = data%thermo%E_vib*HARTREE_TO_CALMOL
      H_rot_cal = data%thermo%E_rot*HARTREE_TO_CALMOL
      H_trans_cal = data%thermo%E_trans*HARTREE_TO_CALMOL + pV_cal
      H_total_cal = H_vib_cal + H_rot_cal + H_trans_cal
      H_int_cal = H_vib_cal + H_rot_cal
      Cp_trans = data%thermo%Cv_trans + R_CALMOLK
      Cv_int = data%thermo%Cv_vib + data%thermo%Cv_rot
      Cv_total = Cp_trans + data%thermo%Cv_rot + data%thermo%Cv_vib
      S_int = data%thermo%S_vib + data%thermo%S_rot
      S_int_J = S_int*CAL_TO_J
      S_total = data%thermo%S_trans + data%thermo%S_rot + data%thermo%S_vib + data%thermo%S_elec
      S_total_J = S_total*CAL_TO_J

      call json%create_object(table_obj, "contribution_table")
      call json%add(thermo_obj, table_obj)

      call json%create_object(vib_row, "VIB")
      call json%add(table_obj, vib_row)
      call json%add(vib_row, "H_cal_mol", H_vib_cal)
      call json%add(vib_row, "Cp_cal_mol_K", data%thermo%Cv_vib)
      call json%add(vib_row, "S_cal_mol_K", data%thermo%S_vib)
      call json%add(vib_row, "S_J_mol_K", data%thermo%S_vib*CAL_TO_J)

      call json%create_object(rot_row, "ROT")
      call json%add(table_obj, rot_row)
      call json%add(rot_row, "H_cal_mol", H_rot_cal)
      call json%add(rot_row, "Cp_cal_mol_K", data%thermo%Cv_rot)
      call json%add(rot_row, "S_cal_mol_K", data%thermo%S_rot)
      call json%add(rot_row, "S_J_mol_K", data%thermo%S_rot*CAL_TO_J)

      call json%create_object(int_row, "INT")
      call json%add(table_obj, int_row)
      call json%add(int_row, "H_cal_mol", H_int_cal)
      call json%add(int_row, "Cp_cal_mol_K", Cv_int)
      call json%add(int_row, "S_cal_mol_K", S_int)
      call json%add(int_row, "S_J_mol_K", S_int_J)

      call json%create_object(tr_row, "TR")
      call json%add(table_obj, tr_row)
      call json%add(tr_row, "H_cal_mol", H_trans_cal)
      call json%add(tr_row, "Cp_cal_mol_K", Cp_trans)
      call json%add(tr_row, "S_cal_mol_K", data%thermo%S_trans)
      call json%add(tr_row, "S_J_mol_K", data%thermo%S_trans*CAL_TO_J)

      call json%create_object(tot_row, "TOT")
      call json%add(table_obj, tot_row)
      call json%add(tot_row, "H_cal_mol", H_total_cal)
      call json%add(tot_row, "Cp_cal_mol_K", Cv_total)
      call json%add(tot_row, "S_cal_mol_K", S_total)
      call json%add(tot_row, "S_J_mol_K", S_total_J)

      ! Zero-point energy
      call json%add(thermo_obj, "zero_point_energy_hartree", data%thermo%zpe_hartree)
      call json%add(thermo_obj, "zero_point_energy_kcal_mol", data%thermo%zpe_kcalmol)

      ! Thermal corrections
      call json%create_object(thermal_obj, "thermal_corrections_hartree")
      call json%add(thermo_obj, thermal_obj)
      call json%add(thermal_obj, "to_energy", data%thermo%thermal_correction_energy)
      call json%add(thermal_obj, "to_enthalpy", data%thermo%thermal_correction_enthalpy)
      call json%add(thermal_obj, "to_gibbs", data%thermo%thermal_correction_gibbs)

      ! Total energies
      call json%create_object(total_e_obj, "total_energies_hartree")
      call json%add(thermo_obj, total_e_obj)
      call json%add(total_e_obj, "electronic", data%total_energy)
      call json%add(total_e_obj, "electronic_plus_zpe", data%total_energy + data%thermo%zpe_hartree)
      call json%add(total_e_obj, "electronic_plus_thermal_E", data%total_energy + data%thermo%thermal_correction_energy)
      call json%add(total_e_obj, "electronic_plus_thermal_H", data%total_energy + data%thermo%thermal_correction_enthalpy)
      call json%add(total_e_obj, "electronic_plus_thermal_G", data%total_energy + data%thermo%thermal_correction_gibbs)

      ! Write to file
      call logger%info("Writing vibrational/thermochemistry JSON to "//trim(output_file))
      open (newunit=iunit, file=trim(output_file), status="replace", action="write", iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         call json%destroy(root)
         return
      end if

      call json%print(root, iunit)
      close (iunit)

      call json%destroy(root)
      call logger%info("Vibrational/thermochemistry JSON written successfully to "//trim(output_file))

   end subroutine write_vibrational_json_impl

   subroutine write_ieda_section(json, parent, data)
      !! The energy resolved onto atoms and atom pairs, under `bonding_energy`
      !!
      !! Written whenever `properties.bonding_analysis.energy_decomposition`
      !! ran. The tables are printed too, but a printed table is not a number a
      !! script can act on.
      !!
      !! **Atom indices are 0-based**, as everywhere else in the interfaces,
      !! and the pairs are written once each with `i < j` carrying the full
      !! pair energy -- the (A,B)/(B,A) doubling of the internal matrices is an
      !! internal convention and there is no reason to export it. So the total
      !! is `sum(atoms) + sum(pairs)` here, with no factor of a half.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      type(json_output_data_t), intent(in) :: data

      type(json_value), pointer :: ieda_obj, atoms_arr, atom_obj, pairs_arr, pair_obj
      type(json_value), pointer :: indices_arr
      integer :: i, j, natm

      if (.not. data%has_ieda) return
      if (.not. allocated(data%ieda_atom)) return
      natm = size(data%ieda_atom)

      call json%create_object(ieda_obj, "bonding_energy")
      call json%add(parent, ieda_obj)
      call json%add(ieda_obj, "energy_of_formation", data%ieda_formation)

      call json%create_array(atoms_arr, "atoms")
      call json%add(ieda_obj, atoms_arr)
      do i = 1, natm
         call json%create_object(atom_obj, "")
         call json%add(atoms_arr, atom_obj)
         call json%add(atom_obj, "index", i - 1)
         call json%add(atom_obj, "energy", data%ieda_atom(i))
         if (allocated(data%ieda_free_atom)) then
            call json%add(atom_obj, "free_atom", data%ieda_free_atom(i))
            ! What it cost to prepare this atom in the shape the molecule
            ! needs. Written out rather than left as a subtraction the reader
            ! has to know to make.
            call json%add(atom_obj, "adaptation", &
                          data%ieda_atom(i) - data%ieda_free_atom(i))
         end if
      end do

      if (.not. allocated(data%ieda_pair)) return
      call json%create_array(pairs_arr, "pairs")
      call json%add(ieda_obj, pairs_arr)
      do i = 1, natm
         do j = i + 1, natm
            call json%create_object(pair_obj, "")
            call json%add(pairs_arr, pair_obj)
            call json%create_array(indices_arr, "atoms")
            call json%add(pair_obj, indices_arr)
            call json%add(indices_arr, "", i - 1)
            call json%add(indices_arr, "", j - 1)
            call json%add(pair_obj, "energy", data%ieda_pair(i, j))
            if (allocated(data%ieda_classical)) then
               call json%add(pair_obj, "classical", data%ieda_classical(i, j))
               ! The part no electrostatic model can produce, which is the
               ! claim the analysis exists to make.
               call json%add(pair_obj, "interference", &
                             data%ieda_pair(i, j) - data%ieda_classical(i, j))
            end if
         end do
      end do
   end subroutine write_ieda_section

   subroutine write_efmo_section(json, parent, data)
      !! The EFMO breakdown, under `efmo`, when there is one
      !!
      !! Named rather than a bare array, as the SAPT one is: the order is fixed
      !! in `EFMO_TERM_NAMES` and nothing downstream should have to know it. The
      !! dimer counts go alongside, because the same system at two cutoffs gives
      !! two different splits and the totals alone do not say which ran.
      use mqc_program_limits, only: N_EFMO_TERMS, EFMO_TERM_NAMES
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      type(json_output_data_t), intent(in) :: data

      type(json_value), pointer :: efmo_obj
      integer :: i

      if (.not. data%has_efmo) return
      if (.not. allocated(data%efmo_terms)) return
      if (size(data%efmo_terms) /= N_EFMO_TERMS) return

      call json%create_object(efmo_obj, "efmo")
      call json%add(parent, efmo_obj)
      do i = 1, N_EFMO_TERMS
         call json%add(efmo_obj, trim(EFMO_TERM_NAMES(i)), data%efmo_terms(i))
      end do
      call json%add(efmo_obj, "qm_dimers", data%efmo_qm_dimers)
      call json%add(efmo_obj, "efp_dimers", data%efmo_efp_dimers)
   end subroutine write_efmo_section

   subroutine write_sapt_section(json, parent, data)
      !! The decomposition, under `sapt`, when there is one
      !!
      !! Named rather than a bare array: the order is fixed in
      !! `SAPT2_TERM_NAMES` and nothing downstream should have to know it.
      ! TODO(mqc): `SAPT_TERM_NAMES` is imported and never used -- both levels
      ! are labelled from `SAPT2_TERM_NAMES`, whose first terms it duplicates.
      use mqc_program_limits, only: N_SAPT_TERMS, SAPT_TERM_NAMES, &
                                    N_SAPT2_TERMS, SAPT2_TERM_NAMES
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: parent
      type(json_output_data_t), intent(in) :: data

      type(json_value), pointer :: sapt_obj
      integer :: i

      if (.not. data%has_sapt) return
      if (.not. allocated(data%sapt_terms)) return
      ! The term count says which level ran: SAPT2 reports every SAPT0 term
      ! in the same slots and appends its own, so the first twelve names are
      ! shared and a SAPT0 consumer can read a SAPT2 output unchanged.
      if (size(data%sapt_terms) /= N_SAPT_TERMS .and. &
          size(data%sapt_terms) /= N_SAPT2_TERMS) return

      call json%create_object(sapt_obj, "sapt")
      call json%add(parent, sapt_obj)
      do i = 1, size(data%sapt_terms)
         call json%add(sapt_obj, trim(SAPT2_TERM_NAMES(i)), data%sapt_terms(i))
      end do
   end subroutine write_sapt_section

end module mqc_json_writer
