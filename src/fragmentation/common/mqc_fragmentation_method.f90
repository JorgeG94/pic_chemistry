!! Which fragmentation method a deck asked for
module mqc_fragmentation_method
   !! **`keywords.fragmentation.method` was required, validated for presence,
   !! and then never read.**
   !!
   !! The expansion was chosen by two other keys instead: `expansion`, taking
   !! `mbe`/`ee-mbe`/`fmo`, and the separate boolean
   !! `allow_overlapping_fragments` for GMBE. Neither was checked against a
   !! list of values either, so a misspelled `expansion` fell through to a
   !! plain MBE without a word.
   !!
   !! What that produced is visible in this repository's own decks: all 33
   !! fragmented cases say `method: "MBE"`, including the five that run FMO or
   !! EE-MBE and the four that run GMBE. The file says one thing and the run
   !! does another, and `method: "wibble"` was accepted just as readily --
   !! measured, on `fmo3_water3.json`, which returns the same FMO energy to
   !! twelve figures whatever that key says.
   !!
   !! So the key means something now. It names the method, it is checked
   !! against the list below, and `expansion` and `allow_overlapping_fragments`
   !! are derived from it rather than read independently.
   !!
   !! ### Adding one
   !!
   !! A new expansion needs an integer here, a `case` in each of the three
   !! procedures, and nothing else -- the adapter reads the mapping rather than
   !! branching on a name. A name that is known and not yet built stays in
   !! `parse_fragmentation_method` while failing `fragmentation_method_implemented`,
   !! so a deck asking for it gets "not implemented yet" rather than "unknown
   !! method", which are different facts. `efmo` was that until Phase 2 built
   !! the orchestrator; nothing is in that state today.
   use pic_types, only: dp
   implicit none
   private

   public :: FRAG_METHOD_UNKNOWN, FRAG_METHOD_MBE, FRAG_METHOD_EE_MBE
   public :: FRAG_METHOD_GMBE, FRAG_METHOD_FMO, FRAG_METHOD_EFMO
   public :: parse_fragmentation_method, fragmentation_method_name
   public :: fragmentation_method_expansion, fragmentation_method_overlapping
   public :: fragmentation_method_implemented, fragmentation_method_list

   integer, parameter :: FRAG_METHOD_UNKNOWN = 0
      !! A spelling this module does not know.
   integer, parameter :: FRAG_METHOD_MBE = 1
      !! The plain many-body expansion over disjoint fragments.
   integer, parameter :: FRAG_METHOD_EE_MBE = 2
      !! Electrostatically embedded MBE: each term sees the point charges of
      !! the fragments it does not contain.
   integer, parameter :: FRAG_METHOD_GMBE = 3
      !! The generalized MBE, over OVERLAPPING fragments, with the
      !! intersection terms supplied by inclusion-exclusion.
   integer, parameter :: FRAG_METHOD_FMO = 4
      !! The fragment molecular orbital method: fragment SCFs converged in each
      !! other's electrostatic potential, then n-mer corrections on top.
   integer, parameter :: FRAG_METHOD_EFMO = 5
      !! Effective fragment molecular orbital: in-vacuo monomers and near
      !! dimers, effective fragment potentials for the far pairs, and one
      !! many-body induction over every fragment. It takes neither the MBE nor
      !! the FMO path -- `expansion_kind` is `"efmo"` and the driver dispatches
      !! it whole to `run_czt_efmo`.

contains

   pure function fragmentation_method_list() result(text)
      !! The accepted spellings, for a refusal message
      !!
      !! One place, so that a message naming the options cannot drift from the
      !! options -- which is the failure mode of writing the list out at each
      !! `error%set`.
      character(len=:), allocatable :: text

      text = "mbe, ee-mbe, gmbe, fmo, efmo"
   end function fragmentation_method_list

   pure subroutine parse_fragmentation_method(name, method, ok)
      !! A deck spelling to a `FRAG_METHOD_*` kind
      !!
      !! Case-insensitive, and `_` reads as `-`. Every deck in this repository
      !! writes `MBE` or `mbe`, and the hyphen in `ee-mbe` is the kind of thing
      !! a person types either way; refusing over that would be pedantry.
      character(len=*), intent(in) :: name
      integer, intent(out) :: method
      logical, intent(out) :: ok

      character(len=len(name)) :: lower
      integer :: i, c

      lower = name
      do i = 1, len_trim(lower)
         c = iachar(lower(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) lower(i:i) = achar(c + 32)
         if (lower(i:i) == "_") lower(i:i) = "-"
      end do

      ok = .true.
      select case (trim(adjustl(lower)))
      case ("mbe")
         method = FRAG_METHOD_MBE
      case ("ee-mbe", "eembe")
         method = FRAG_METHOD_EE_MBE
      case ("gmbe")
         method = FRAG_METHOD_GMBE
      case ("fmo")
         method = FRAG_METHOD_FMO
      case ("efmo")
         method = FRAG_METHOD_EFMO
      case default
         method = FRAG_METHOD_UNKNOWN
         ok = .false.
      end select
   end subroutine parse_fragmentation_method

   pure function fragmentation_method_implemented(method) result(yes)
      !! Whether an expansion exists behind the name
      integer, intent(in) :: method
      logical :: yes

      yes = method == FRAG_METHOD_MBE .or. method == FRAG_METHOD_EE_MBE .or. &
            method == FRAG_METHOD_GMBE .or. method == FRAG_METHOD_FMO .or. &
            method == FRAG_METHOD_EFMO
   end function fragmentation_method_implemented

   pure function fragmentation_method_name(method) result(name)
      !! The canonical spelling, for a message or an echo
      integer, intent(in) :: method
      character(len=:), allocatable :: name

      select case (method)
      case (FRAG_METHOD_MBE)
         name = "mbe"
      case (FRAG_METHOD_EE_MBE)
         name = "ee-mbe"
      case (FRAG_METHOD_GMBE)
         name = "gmbe"
      case (FRAG_METHOD_FMO)
         name = "fmo"
      case (FRAG_METHOD_EFMO)
         name = "efmo"
      case default
         name = "unknown"
      end select
   end function fragmentation_method_name

   pure function fragmentation_method_expansion(method) result(kind)
      !! The `expansion_kind` the driver branches on
      !!
      !! GMBE maps to `"mbe"` and is distinguished by the overlapping flag
      !! below, which is how the driver has always told them apart -- the two
      !! share a term builder and differ in whether primaries may overlap.
      integer, intent(in) :: method
      character(len=16) :: kind

      select case (method)
      case (FRAG_METHOD_EE_MBE)
         kind = "ee-mbe"
      case (FRAG_METHOD_FMO)
         kind = "fmo"
      case (FRAG_METHOD_EFMO)
         ! Not an expansion the term builders can run: EFMO has no n-mer list
         ! and no subset closure, so the driver sends it whole to the backend
         ! before the fragmented path is entered.
         kind = "efmo"
      case default
         kind = "mbe"
      end select
   end function fragmentation_method_expansion

   pure function fragmentation_method_overlapping(method) result(yes)
      !! Whether primaries may overlap, which is what makes it GMBE
      integer, intent(in) :: method
      logical :: yes

      yes = method == FRAG_METHOD_GMBE
   end function fragmentation_method_overlapping

end module mqc_fragmentation_method
