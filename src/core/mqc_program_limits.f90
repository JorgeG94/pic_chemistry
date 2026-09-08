!! Program limits and default parameter, publics
module mqc_program_limits
   !! Compile-time limits and defaults that control memory allocation and
   !! algorithm behaviour.
   use pic_types, only: dp
   implicit none
   private

   !---------------------------------------------------------------------------
   ! Many-Body Expansion Limits
   !---------------------------------------------------------------------------

   integer, parameter, public :: MAX_MBE_LEVEL = 10
   !! Maximum MBE truncation order (1-body, 2-body, ..., N-body)
   !! Higher orders require factorial growth in fragment combinations

   integer, parameter, public :: N_EFP_TERMS = 6
   !! Terms in an EFP2 interaction energy: electrostatics, polarization,
   !! exchange repulsion, dispersion, charge transfer, and their total.

   integer, parameter, public :: N_EFMO_TERMS = 8
   !! The sums an EFMO energy is reported as; `EFMO_TERM_NAMES` names them in
   !! slot order. Eight rather than the six of eq 6 because the far half is
   !! reported term by term: the four effective-fragment terms are what a
   !! comparison against another code is made of, and adding them up first
   !! would hide which one disagreed.

   character(len=*), parameter, public :: EFMO_TERM_NAMES(N_EFMO_TERMS) = &
                                          [character(len=24) :: "monomer_sum", &
                                                                 "qm_dimer_correction", "pair_polarization", &
                                                                 "efp_electrostatics", "efp_dispersion", &
                                                                 "efp_exchange_repulsion", "efp_charge_transfer", &
                                                                 "polarization_total"]
   !! The slot each EFMO sum occupies, and the name it is written out under.
   !!
   !! `pair_polarization` is `sum_IJ E_IJ^pol` over the quantum dimers and is
   !! **subtracted** from the total; it is reported with the sign it has as a
   !! sum, not with the sign it enters with, so that it can be compared against
   !! another code's pair induction directly.

   integer, parameter, public :: N_SAPT_TERMS = 12
   !! Terms a SAPT0 interaction energy is reported as; `SAPT_TERM_NAMES` names
   !! them in slot order.

   character(len=*), parameter, public :: SAPT_TERM_NAMES(N_SAPT_TERMS) = &
                                          [character(len=12) :: "elst10", "exch10_s2", "exch10", &
                                                                 "ind20_u", "ind20_r", "exch_ind20_u", "exch_ind20_r", &
                                                              "disp20", "exch_disp20", "delta_hf", "e_int_hf_cp", "total"]
   !! The slot each SAPT term occupies, and the name it is written out under.
   !! The names are the literature's rather than the console's prose, since
   !! these exist to be compared against another code's output.
   !!
   !! `_u` is uncoupled and `_r` is response (coupled); `_s2` is the single
   !! exchange approximation against the full `S^inf`. A code that reports only
   !! one of each pair cannot be checked against one that reports the other.

   integer, parameter, public :: N_SAPT2_TERMS = 18
   !! SAPT2 reports every SAPT0 term in the same slots -- "total" stays the
   !! SAPT0 total, so the two levels' outputs line up term for term -- and
   !! appends the four intramonomer-correlation corrections, the scaled
   !! exchange-induction, and its own total.

   character(len=*), parameter, public :: SAPT2_TERM_NAMES(N_SAPT2_TERMS) = &
                                          [SAPT_TERM_NAMES, "elst12      ", "exch11      ", "exch12      ", &
                                           "ind22       ", "exch_ind22  ", "total_sapt2 "]

   integer, parameter, public :: GROUP_RESULT_BATCH_SIZE = 256
   !! Group-global result batching size for MPI MBE (multi-global coordinator)

   !---------------------------------------------------------------------------
   ! In-Core Integral Budgets
   !---------------------------------------------------------------------------

   ! How much memory a rank will spend on stored two-electron integrals.

   real(dp), parameter, public :: ERI_CORE_BUDGET_CAP = 2.0e9_dp
      !! Ceiling on what one rank may spend on stored two-electron integrals.
      !!
      !! The tensor is the full n^4, not the eightfold-unique set, so the cap is
      !! reached at 128 functions.

   real(dp), parameter, public :: ERI_CORE_BUDGET_SHARE = 0.25_dp
      !! Fraction of *currently available* memory a rank will claim.
      !!
      !! Available rather than installed, so several ranks on one node degrade
      !! toward smaller claims -- a rank deciding later sees what the earlier
      !! ones already took -- without needing a count that is not discoverable
      !! here. Ranks deciding at the same instant see the same number, which is
      !! why the share is a quarter and not a half. Claiming too much is fatal;
      !! claiming too little falls back to a direct build, which is slower and
      !! correct.

   real(dp), parameter, public :: ERI_CORE_BUDGET_BLIND = 5.0e8_dp
      !! What to allow when available memory cannot be read at all.
      !!
      !! /proc/meminfo is Linux; macOS and anything else land here. Half a
      !! gigabyte is 94 functions, which keeps the small fragments a
      !! fragmented run is made of while refusing to guess on a machine this
      !! cannot measure.

   real(dp), parameter, public :: SAPT_CORE_BUDGET_SHARE = 0.8_dp
      !! Fraction of available memory a SAPT run may size itself against.
      !!
      !! Far above `ERI_CORE_BUDGET_SHARE` because SAPT has no direct path:
      !! every term is a contraction over the stored dimer tensor, so refusing
      !! the memory refuses the calculation, where an SCF only gets slower. The
      !! guard turns an OOM kill into a message naming the basis, and should
      !! refuse only what genuinely cannot run.

   !---------------------------------------------------------------------------
   ! Density Fitting
   !---------------------------------------------------------------------------

   real(dp), parameter, public :: DF_METRIC_PANEL_BYTES = 8.0e6_dp
      !! Working size of one row panel of the fitted-tensor metric contraction.
      !!
      !! `B = (mn|Q) J^(-1/2)` is split over the pair index so that each thread
      !! reads its own slice of the three-centre tensor rather than all of it,
      !! and this is how tall a slice is. The panel is packed per thread and
      !! lives for the whole call, so the cost is this times the thread count.
   ! TODO(mqc): 8 MB here, where the sizing argument that came with it was
   ! written for two -- a hundred threads hold 800 MB of panel, not the 200 MB
   ! the choice was justified by. Value and reasoning disagree by four times.

   real(dp), parameter, public :: DF_PAIR_SCREEN = 1.0e-12_dp
      !! Below this, a shell pair contributes no three-centre integral.
      !!
      !! Schwarz: `|(mn|P)| <= sqrt((mn|mn)) sqrt((P|P))`, so a pair whose
      !! bound times the largest auxiliary diagonal falls under this cannot
      !! reach it, for any P, and the whole shell triplet is skipped.
      !!
      !! 1e-12 rather than the 1e-10 the literature usually quotes, because the
      !! validation suite compares total energies at 1e-9 and a fitted energy
      !! has already spent its error budget on the fit.

   integer, parameter, public :: DF_AUX_CHUNK = 32
      !! Auxiliary functions per thread chunk in the fitted Coulomb build.
      !!
      !! Both halves of J are BLAS-2 over a column block of B, and this is how
      !! wide a block is: small enough that a few thousand auxiliary functions
      !! still make many more chunks than there are threads, large enough that
      !! the per-call BLAS overhead is not paid n_aux times.

   !---------------------------------------------------------------------------
   ! Numerical Differentiation Defaults
   !---------------------------------------------------------------------------

   real(dp), parameter, public :: DEFAULT_FD_DISPLACEMENT = 0.005_dp
   !! Default step size for finite difference calculations (Bohr)
   !! ~0.005 Bohr = ~0.0026 Angstrom, suitable for Hessian/gradient FD

   !---------------------------------------------------------------------------
   ! I/O Limits
   !---------------------------------------------------------------------------

   integer, parameter, public :: MAX_LINE_LENGTH = 1024
   !! Maximum length for input file lines

   integer, parameter, public :: MAX_ORBITAL_LABEL_LEN = 8
      !! "Cr 3d" and its longest relatives. An atomic orbital label is an
      !! element symbol, a space, a principal quantum number and a subshell
      !! letter, so eight characters is generous.

   integer, parameter, public :: MAX_ELEMENT_SYMBOL_LEN = 4
   !! Maximum length for element symbols (e.g., "He", "Uue")
   ! TODO(mqc): redefined privately, and identically, in `mqc_geometry` and
   ! `mqc_xyz_reader` rather than used from here -- three spellings of one
   ! constant, which is what MQC002 exists to prevent.

   character(len=*), parameter, public :: JSON_REAL_FORMAT = "ES"
   !! JSON output format for real numbers (scientific notation)
   !! Valid values: 'G', 'E', 'EN', 'ES' (json-fortran uses machine precision)

   !---------------------------------------------------------------------------
   ! Geometry/Structure Limits
   !---------------------------------------------------------------------------

   real(dp), parameter, public :: MIN_ATOM_DISTANCE = 0.01_dp
   !! Minimum allowed distance between atoms (Bohr)
   !! Atoms closer than this are considered overlapping (error condition)

end module mqc_program_limits
