!! Atomic radii, one module for every hand-maintained table
module mqc_atomic_radii
   !! The atomic radii the program uses, named for their source so that a caller
   !! picks a parametrisation deliberately.
   !!
   !! There is no single "atomic radius". These tables were fitted to different
   !! things and disagree with each other by tenths of an Angstrom, so swapping
   !! one for another silently changes results.
   !!
   !! | Accessor                  | Source                     | Fitted for            |
   !! |---------------------------|----------------------------|-----------------------|
   !! | `covalent_radius_cordero` | Cordero (2008), Z <= 96    | bond perception       |
   !! | `covalent_radius_emsley`  | Emsley, as GAMESS, Z <= 36 | DMA bond midpoints    |
   !! | `vdw_radius_bondi`        | Bondi/Mantina, Z <= 18     | PCM cavity            |
   !! | `vdw_radius_geodesic`     | GAMESS geodesic, Z <= 17   | ESP screening grid    |
   !! | `vdw_radius_fmo`          | GAMESS `$FMO VDWRAD`       | FMO/EFMO pair cutoffs |
   !!
   !! The Bragg-Slater radii that size the DFT partition cells are *not* here:
   !! they live in the generated `mqc_dft_radial_data`, reached through
   !! `mqc_dft_radial`'s `bragg_radius`.
   !!
   !! Every accessor returns Angstrom, and reports an unknown element rather
   !! than substituting a number -- an invented radius shows up as an invented
   !! bond or an invented cavity, never as an error.
   use pic_types, only: dp
   implicit none
   private

   public :: covalent_radius_cordero   !! Cordero covalent radius, Angstrom
   public :: covalent_radius_emsley    !! Emsley/GAMESS covalent radius, Angstrom
   public :: vdw_radius_bondi          !! Bondi/Mantina van der Waals radius, Angstrom
   public :: vdw_radius_geodesic       !! GAMESS geodesic van der Waals radius, Angstrom
   public :: vdw_radius_fmo            !! GAMESS FMO van der Waals radius, Angstrom
   public :: MAX_Z_CORDERO, MAX_Z_EMSLEY, MAX_Z_BONDI, MAX_Z_GEODESIC, MAX_Z_FMO
   public :: GEODESIC_RADIUS_DEFAULT, FMO_RADIUS_DEFAULT

   integer, parameter :: MAX_Z_CORDERO = 96
      !! Cordero's set stops at curium; past it there are no consensus radii.

   integer, parameter :: MAX_Z_EMSLEY = 36
      !! Highest element with an Emsley covalent radius here, so a bond can be
      !! decided.

   integer, parameter :: MAX_Z_BONDI = 18

   integer, parameter :: MAX_Z_GEODESIC = 17

   integer, parameter :: MAX_Z_FMO = 35
      !! Past bromine GAMESS's FMO table is `FMO_RADIUS_DEFAULT` throughout, so
      !! there is nothing to tabulate.

   real(dp), parameter :: GEODESIC_RADIUS_DEFAULT = 1.8_dp
   !! What an element absent from the geodesic table takes, as GAMESS does.

   real(dp), parameter :: CORDERO_RADII(MAX_Z_CORDERO) = [ &
      !! Covalent radii in Angstrom, from Cordero et al., Dalton Trans. 2008, 2832
      !! -- the set most codes use for distance-based bond perception.
                          0.31_dp, 0.28_dp, 1.28_dp, 0.96_dp, 0.84_dp, 0.76_dp, 0.71_dp, 0.66_dp, &  ! Z = 1-8
                          0.57_dp, 0.58_dp, 1.66_dp, 1.41_dp, 1.21_dp, 1.11_dp, 1.07_dp, 1.05_dp, &  ! Z = 9-16
                          1.02_dp, 1.06_dp, 2.03_dp, 1.76_dp, 1.70_dp, 1.60_dp, 1.53_dp, 1.39_dp, &  ! Z = 17-24
                          1.39_dp, 1.32_dp, 1.26_dp, 1.24_dp, 1.32_dp, 1.22_dp, 1.22_dp, 1.20_dp, &  ! Z = 25-32
                          1.19_dp, 1.20_dp, 1.20_dp, 1.16_dp, 2.20_dp, 1.95_dp, 1.90_dp, 1.75_dp, &  ! Z = 33-40
                          1.64_dp, 1.54_dp, 1.47_dp, 1.46_dp, 1.42_dp, 1.39_dp, 1.45_dp, 1.44_dp, &  ! Z = 41-48
                          1.42_dp, 1.39_dp, 1.39_dp, 1.38_dp, 1.39_dp, 1.40_dp, 2.44_dp, 2.15_dp, &  ! Z = 49-56
                          2.07_dp, 2.04_dp, 2.03_dp, 2.01_dp, 1.99_dp, 1.98_dp, 1.98_dp, 1.96_dp, &  ! Z = 57-64
                          1.94_dp, 1.92_dp, 1.92_dp, 1.89_dp, 1.90_dp, 1.87_dp, 1.87_dp, 1.75_dp, &  ! Z = 65-72
                          1.70_dp, 1.62_dp, 1.51_dp, 1.44_dp, 1.41_dp, 1.36_dp, 1.36_dp, 1.32_dp, &  ! Z = 73-80
                          1.45_dp, 1.46_dp, 1.48_dp, 1.40_dp, 1.50_dp, 1.50_dp, 2.60_dp, 2.21_dp, &  ! Z = 81-88
                          2.15_dp, 2.06_dp, 2.00_dp, 1.96_dp, 1.90_dp, 1.87_dp, 1.80_dp, 1.69_dp]  ! Z = 89-96

   real(dp), parameter :: EMSLEY_RADII(MAX_Z_EMSLEY) = [ &
      !! Covalent radii in Angstrom, Emsley's table as GAMESS uses it, indexed by Z.
      !!
      !! A bond is `|Ri - Rj| <= rad_i + rad_j`, evaluated in Angstrom, which decides
      !! where the bond midpoints go and therefore how many expansion points there
      !! are. Anything past krypton is refused rather than given an invented radius.
                          0.32_dp, 0.93_dp, &                                     ! H  He
                          1.23_dp, 0.90_dp, 0.82_dp, 0.77_dp, 0.75_dp, 0.73_dp, 0.72_dp, 0.71_dp, &  ! Li..Ne
                          1.54_dp, 1.36_dp, 1.18_dp, 1.11_dp, 1.06_dp, 1.02_dp, 0.99_dp, 0.98_dp, &  ! Na..Ar
                          2.03_dp, 1.74_dp, 1.44_dp, 1.32_dp, 1.22_dp, 1.18_dp, 1.17_dp, 1.17_dp, &  ! K..Fe
                          1.16_dp, 1.15_dp, 1.17_dp, 1.25_dp, 1.26_dp, 1.22_dp, 1.20_dp, 1.16_dp, &  ! Co..Se
                          1.14_dp, 1.12_dp]                                       ! Br Kr

   real(dp), parameter :: BONDI_RADII(MAX_Z_BONDI) = [ &
      !! Van der Waals radii in Angstrom, H through Ar, indexed by atomic number.
      !!
      !! Bondi (1964) except Li, Be, B, Na, Mg and Al, which are Mantina (2009).
                          1.20_dp, &   ! H
                          1.40_dp, &   ! He
                          1.82_dp, &   ! Li
                          1.53_dp, &   ! Be   (Mantina)
                          1.92_dp, &   ! B    (Mantina)
                          1.70_dp, &   ! C
                          1.55_dp, &   ! N
                          1.52_dp, &   ! O
                          1.47_dp, &   ! F
                          1.54_dp, &   ! Ne
                          2.27_dp, &   ! Na   (Mantina)
                          1.73_dp, &   ! Mg   (Mantina)
                          1.84_dp, &   ! Al   (Mantina)
                          2.10_dp, &   ! Si
                          1.80_dp, &   ! P
                          1.80_dp, &   ! S
                          1.75_dp, &   ! Cl
                          1.88_dp]     ! Ar

   real(dp), parameter :: GEODESIC_RADII(MAX_Z_GEODESIC) = [ &
      !! Geodesic van der Waals radii, Angstrom, indexed by Z. Zero means absent,
      !! and an absent element takes `GEODESIC_RADIUS_DEFAULT`, as GAMESS does.
                          1.20_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.85_dp, 1.50_dp, 1.50_dp, &
                          1.40_dp, 1.35_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.07_dp, 2.05_dp, &
                          1.96_dp, 1.89_dp, 1.80_dp]

   real(dp), parameter :: FMO_RADIUS_DEFAULT = 2.5_dp
      !! What GAMESS's `vdwr0` carries for every element it does not name.

   real(dp), parameter :: FMO_RADII(MAX_Z_FMO) = [ &
      !! GAMESS's `$FMO VDWRAD` defaults (`vdwr0` in `fmoio.f90`), Angstrom.
      !!
      !! **Not Bondi**, and the difference is not cosmetic: Bondi puts hydrogen
      !! at 1.10 and oxygen at 1.52, this table at 1.20 and 1.40, so the same
      !! water pair comes out at a different multiple of contact and a cutoff of
      !! 2.0 falls on the other side of it. The FMO literature's `RESPPC`,
      !! `RESDIM` and `R_cut` values are all quoted against *this* table.
                          1.20_dp, 1.20_dp, 1.37_dp, 1.45_dp, 1.45_dp, 1.50_dp, 1.50_dp, &
                          1.40_dp, 1.35_dp, 1.30_dp, 1.57_dp, 1.36_dp, 1.24_dp, 1.17_dp, &
                          1.80_dp, 1.75_dp, 1.70_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, &
                          2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, &
                          2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.50_dp, 2.30_dp]

contains

   pure function covalent_radius_cordero(atomic_number) result(radius)
      !! Cordero covalent radius in Angstrom, or 0 where none is tabulated
      !!
      !! Zero is a refusal, not a default: a caller doing bond perception has to
      !! decide for itself what an element with no radius means.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      select case (atomic_number)
      case (1:MAX_Z_CORDERO)
         radius = CORDERO_RADII(atomic_number)
      case default
         radius = 0.0_dp
      end select
   end function covalent_radius_cordero

   pure function covalent_radius_emsley(atomic_number) result(radius)
      !! Emsley/GAMESS covalent radius in Angstrom, or 0 where none is tabulated
      !!
      !! Zero means the caller must refuse the element rather than place a bond
      !! midpoint from a radius that was never measured.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      select case (atomic_number)
      case (1:MAX_Z_EMSLEY)
         radius = EMSLEY_RADII(atomic_number)
      case default
         radius = 0.0_dp
      end select
   end function covalent_radius_emsley

   pure function vdw_radius_bondi(atomic_number) result(radius)
      !! Bondi/Mantina van der Waals radius in Angstrom, or 0 if untabulated
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      select case (atomic_number)
      case (1:MAX_Z_BONDI)
         radius = BONDI_RADII(atomic_number)
      case default
         radius = 0.0_dp
      end select
   end function vdw_radius_bondi

   pure function vdw_radius_geodesic(atomic_number) result(radius)
      !! GAMESS geodesic van der Waals radius in Angstrom
      !!
      !! Unlike the others this one substitutes `GEODESIC_RADIUS_DEFAULT`, as
      !! GAMESS does: the screening grid is a fitting aid, so a slightly wrong
      !! sphere moves sample points and not an energy. Elements inside the
      !! table's range but absent from it (the noble gases and Li..Be) are zero
      !! there, and take the default too.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      radius = GEODESIC_RADIUS_DEFAULT
      if (atomic_number < 1 .or. atomic_number > MAX_Z_GEODESIC) return
      if (GEODESIC_RADII(atomic_number) <= 0.0_dp) return
      radius = GEODESIC_RADII(atomic_number)
   end function vdw_radius_geodesic

   pure function vdw_radius_fmo(atomic_number) result(radius)
      !! GAMESS FMO van der Waals radius in Angstrom, for the fragment-pair cutoff
      !!
      !! Substitutes `FMO_RADIUS_DEFAULT` for an element outside the table, as
      !! GAMESS does. A separation in units of contact is a screening decision,
      !! not an energy, so a coarse radius costs a pair its class and never a
      !! number; refusing the element instead would refuse every heavy-atom
      !! system GAMESS runs.
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      radius = FMO_RADIUS_DEFAULT
      if (atomic_number < 1 .or. atomic_number > MAX_Z_FMO) return
      radius = FMO_RADII(atomic_number)
   end function vdw_radius_fmo

end module mqc_atomic_radii
