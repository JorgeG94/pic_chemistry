!! Configuration types an input file parses into
module mqc_config_types
   !! `mqc_config_t` and the pieces it is made of: one molecule, one input
   !! fragment, and the bonds between them.
   !!
   !! Every default is declared here, once. A reader that finds no value for a
   !! key leaves the component alone, so there is no second table of defaults.
   use pic_types, only: dp, int32
   use mqc_scf_types, only: guess_step_t, scf_numerics_t, deltascf_options_t
   use mqc_method_types, only: METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY
   use mqc_geometry, only: geometry_type
   use mqc_physical_fragment, only: bond_t
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE, DEFAULT_RESPONSE_TOL, &
                                       DEFAULT_RESPONSE_MAX_ITER, DEFAULT_AIMD_DT, &
                                       DEFAULT_AIMD_NSTEPS, DEFAULT_AIMD_TEMPERATURE, &
                                       DEFAULT_AIMD_OUTPUT_FREQ, DEFAULT_SCF_CONV, &
                                       DEFAULT_SCF_DENSITY_CONV, DEFAULT_VDW_SCALE, &
                                       DEFAULT_DYNAMIC_TOL, DEFAULT_DYNAMIC_MAXITER, &
                                       DEFAULT_RESPONSE_BATCH, &
                                       EFP_RESPONSE_AUTO, &
                                       DEFAULT_OPT_MAX_STEPS, DEFAULT_OPT_GRADIENT_TOLERANCE, &
                                       DEFAULT_OPT_ENERGY_TOLERANCE, DEFAULT_OPT_MAX_STEP, &
                                       DEFAULT_OPT_LBFGS_MEMORY, DEFAULT_OPT_PRINT_LEVEL, &
                                       DEFAULT_CPCM_NANG, DEFAULT_CPCM_RSCALE, &
                                       DEFAULT_PCM_NANG, DEFAULT_PCM_RSCALE, &
                                       DEFAULT_PCM_ZETA, &
                                       DEFAULT_FRAG_LEVEL, DEFAULT_MAX_INTERSECTION
   implicit none
   private

   public :: mqc_config_t       !! Everything one input file specifies
   public :: molecule_t         !! One molecule of a multi-molecule input
   public :: input_fragment_t   !! One fragment definition, as written in the input
   public :: bond_t             !! Re-exported so consumers need only this module

   type :: input_fragment_t
      !! One fragment as the input file declares it: charge, multiplicity and
      !! atom indices
      !!
      !! The parsed representation, not the computational fragment.
      integer :: charge = 0
      integer :: multiplicity = 1
      integer, allocatable :: indices(:)  !! Atom indices in this fragment
      character(len=:), allocatable :: potential
         !! The effective fragment potential describing this fragment, if it has
         !! one. A fragment carrying one is an EFP fragment, one without stays
         !! quantum, which is how a mixed calculation is spelled. The path is
         !! resolved relative to the deck, like `xyz`.
         !!
         !! **A path, not a format.** Nothing here may require an extension or
         !! infer a format from one; whoever loads it dispatches on what the
         !! file actually is.
   contains
      procedure :: destroy => input_fragment_destroy
   end type input_fragment_t

   type :: molecule_t
      !! Single molecule definition with structure, geometry, fragments, and connectivity
      character(len=:), allocatable :: name  !! Optional molecule name

      ! Structure information
      integer :: charge = 0
      integer :: multiplicity = 1

      ! Geometry
      type(geometry_type) :: geometry

      ! Fragments
      integer :: nfrag = 0
      type(input_fragment_t), allocatable :: fragments(:)
      logical :: uniform_system = .false.
         !! Every fragment is the same species, so one potential describes all
         !! of them, which is what lets `fragment_potentials` be a single name.
         !! Checked rather than taken on trust: the fragments must agree in
         !! size, so a deck claiming uniformity without it is refused.

      ! Connectivity
      integer :: nbonds = 0
      integer :: nbroken = 0
      type(bond_t), allocatable :: bonds(:)

   contains
      procedure :: destroy => molecule_destroy
   end type molecule_t

   type :: mqc_config_t
      !! Complete configuration parsed from one input deck

      ! Schema information
      character(len=:), allocatable :: schema_name
      character(len=:), allocatable :: schema_version
      integer :: index_base = 0  !! 0-based or 1-based indexing
      character(len=:), allocatable :: units  !! angstrom or bohr

      ! Properties asked for alongside the energy. Not settings for *how* to
      ! compute -- those are `keywords` -- but requests for something extra to
      ! be reported once the wave function exists.
      character(len=:), allocatable :: bonding_analysis
         !! From `properties.bonding_analysis.type`. Absent or "none" means no
         !! analysis; "gms_quao" means the quasi-atomic bonding picture.
      character(len=:), allocatable :: fukui_population
         !! From `properties.fukui.population`. Allocated means the deck asked
         !! for Fukui indices; "chelpg" (the default) or "mulliken" chooses how
         !! the density difference is condensed onto atoms. The object may omit
         !! it -- `properties: {"fukui": {}}` is a valid request.
      character(len=:), allocatable :: fukui_guess
         !! `properties.fukui.guess`: "neutral" (the default) or "independent".
         !!
         !! "neutral" seeds each ion from the converged neutral's orbitals,
         !! filled Aufbau at that ion's own occupation. Usually much the better
         !! start, but it puts the extra electron in the *neutral's* LUMO, and
         !! where that orbital is a poor guess for the ion's own the SCF can
         !! settle on a different stationary point rather than climb out.
         !! "independent" lets each ion take the ordinary guess instead.
      type(deltascf_options_t) :: fukui_scf
         !! `properties.fukui.scf`: how the ion SCFs converge.
         !!
         !! Seeded from the resolved `keywords.scf` by the reader unless
         !! `inherit_scf` is false, then overwritten by whatever the object
         !! names.
      character(len=:), allocatable :: charges_scheme
         !! From `properties.charges.scheme`. Allocated means the deck asked for
         !! atomic partial charges; "mulliken" (the default) or "chelpg" chooses
         !! how the converged density is partitioned.
      logical :: bonding_energy = .false.
         !! From `properties.bonding_analysis.energy_decomposition`. Off by
         !! default: the two-electron term needs the dense `n_ao^4` integral
         !! array, where the bonding tables on their own need only one-electron
         !! integrals.
      logical :: bonding_no_sharing = .false.
         !! From `properties.bonding_analysis.no_sharing`. Asks for the
         !! no-sharing wave function, which needs a full valence CI over the
         !! quasi-atomic orbitals, so it is off by default.
      logical :: bonding_restrict_localization = .false.
         !! From `properties.bonding_analysis.restrict_localization`. Confine
         !! the quasi-atomic localization to the wave function's
         !! occupation-restricted subspaces, so that no rotation mixes two of
         !! them and the wave function stays invariant under the
         !! transformation. Off by default.
      character(len=:), allocatable :: bonding_no_sharing_ci
         !! From `properties.bonding_analysis.no_sharing_ci`. How the CI
         !! expansion over quasi-atomic orbitals is obtained: `"transform"`
         !! solves in the molecular orbital basis and carries the vector across
         !! with the orbital transformation, `"resolve"` runs a second Davidson
         !! in the quasi-atomic basis. The two give the same wave function and
         !! this chooses how it is computed, not which one it is.
      real(dp) :: bonding_threshold = 1.0_dp
         !! From `properties.bonding_analysis.energy_threshold`, in kcal/mol.
         !! Orbital pairs whose kinetic bond order is weaker than this are
         !! counted and then not printed.

      ! Model information
      integer(int32) :: method = METHOD_TYPE_GFN2
      character(len=:), allocatable :: basis
      character(len=:), allocatable :: ecp
         !! `model.ecp`. Unallocated means no potential, which is not the
         !! same as an empty one: a deck that names a set gets it applied to
         !! whichever atoms the file covers, and a deck that names none gets
         !! an all-electron calculation.
      character(len=:), allocatable :: aux_basis
      logical :: cartesian = .false.
         !! Read the basis in Cartesian form whatever its file declares.
         !!
         !! Not a preference: for a basis with a shell above p the two forms
         !! span different spaces and give different energies, so this changes
         !! the model rather than its representation. Defaults false, which
         !! honours the file.
      character(len=:), allocatable :: functional
         !! XC functional name, only meaningful when method = dft
      logical :: scf_incremental_fock = .true.
         !! `keywords.scf.incremental_fock`. False forces a full Fock build every
         !! iteration; see `scf_config_t` for why anyone would want that.
      character(len=:), allocatable :: scf_accelerator
         !! `keywords.scf.accelerator`: 'diis' (the default), 'adiis' or
         !! 'ediis'. The energy-based pair runs only while the error is large
         !! and hands over to DIIS below `ACCEL_SWITCH`, so naming one asks for
         !! a different opening, not a different endgame.
      character(len=:), allocatable :: scf_eri_path
         !! `keywords.scf.eri_path`: 'rys', 'rotaxis' or 'auto'. Unallocated
         !! when the deck says nothing, which leaves the backend's default.
      character(len=:), allocatable :: scf_guess
         !! Initial guess name from `keywords.scf.guess`.
         !!
         !! Superseded by `keywords.guess.type`. Still read so existing decks
         !! keep working; naming both is an error rather than a precedence rule.
      character(len=:), allocatable :: guess_type
         !! Initial guess name from `keywords.guess.type`, the current spelling.
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The projection ladder, from `keywords.guess.subscf.steps`.
         !!
         !! One entry per preliminary SCF, in order. The target basis is
         !! `model.basis` and is not repeated here, so three SCFs -- STO-3G, then
         !! 6-31G, then cc-pVTZ -- is two steps and a model basis.
      logical :: unrestricted = .false.
         !! `model.unrestricted`. Force UHF/UKS even when the shell is closed.
      logical :: allow_crap_scf = .false.
         !! Let a non-converged SCF into the expansion instead of stopping. Off
         !! by default: the energy of an SCF that ran out of iterations is of
         !! the right magnitude and nothing downstream can tell. On, every
         !! offender is named at the end of the run and marked in the fragment
         !! table.
      logical :: scf_density_fitting = .false.
         !! Fit J and K against the auxiliary basis, on the CPU backend
      ! keywords.correlation -- post-Hartree-Fock, kept apart from the SCF ones
      logical :: corr_freeze_core = .true.
      integer :: corr_n_frozen_core = -1   !! -1 counts the core from the elements
      logical :: corr_density_fitting = .false.  !! RI for the correlation step
      logical :: corr_scs = .false.        !! Spin-component scaling
      real(dp) :: corr_scs_ss = 1.0_dp/3.0_dp
      real(dp) :: corr_scs_os = 1.2_dp
      ! keywords.cc -- what only an iterative correlation method needs
      integer :: cc_maxiter = 100
      real(dp) :: cc_tolerance = 1.0e-8_dp
      logical :: cc_triples_set = .false.
         !! Whether a deck named `triples` explicitly. Without this there is no
         !! way to tell "triples: false" from "the key was absent", and the
         !! method name has to lose to an explicit keyword while winning over an
         !! absent one.
      logical :: cc_triples = .false.
      logical :: cc_diis = .true.
      integer :: cc_diis_size = 8
      logical :: cc_spin_adapted = .true.
         !! Which coupled-cluster formulation runs: spin-adapted (spatial
         !! orbitals) when true, spin orbitals when false.
         !!
         !! The two agree to machine precision for the closed-shell reference
         !! either formulation supports, so this chooses how a number is
         !! computed and not which number. Spin-adapted is the default: about
         !! sixteen times less memory and several times faster.
      ! keywords.mcscf -- the active space, for CASSCF and CASCI
      character(len=:), allocatable :: mcscf_avas_orbitals(:)
         !! Atomic orbital labels from `keywords.mcscf.avas.orbitals`, e.g.
         !! "N 2p". Unallocated means the active space was given by counts.
      logical :: mcscf_full_valence = .false.
         !! `keywords.mcscf.full_valence` -- take the whole valence shell as the
         !! active space, sized by counting the free-atom minimal basis rather
         !! than by any threshold.
      real(dp) :: mcscf_avas_threshold = 0.2_dp
      integer, allocatable :: mcscf_ormas_subspaces(:)
         !! Active orbital each subspace starts at, from
         !! `keywords.mcscf.ormas.subspaces`. Absent leaves a complete active
         !! space, which is the one subspace covering everything.
      integer, allocatable :: mcscf_ormas_min_electrons(:)
      integer, allocatable :: mcscf_ormas_max_electrons(:)
         !! Electrons allowed in each subspace, both spins together
      integer :: mcscf_n_active_electrons = 0
      integer :: mcscf_n_active_orbitals = 0
      integer :: mcscf_n_inactive_orbitals = -1
         !! -1 derives it from the electron count: every electron the active
         !! space does not hold is in a doubly occupied orbital.
      logical :: mcscf_optimize_orbitals = .true.
         !! CASSCF when true, CASCI when false. "casci" and "casscf" parse to
         !! the same method type, so the spelling is the only place the
         !! distinction survives; the reader defaults this from it and lets
         !! `keywords.mcscf.optimize_orbitals` override.
      integer :: mcscf_max_macro_iter = 100
      real(dp) :: mcscf_orbital_convergence = 1.0e-6_dp
      ! keywords.dft -- the quadrature, not the functional
      real(dp) :: dft_screening_tolerance = 1.0e-12_dp
         !! From `keywords.dft.screening_tolerance`. The AO value below which a
         !! shell is dropped from a grid block. Zero or negative evaluates the
         !! whole basis, which is how a run measures what the screen is worth.
      integer :: dft_block_size = -1
         !! From `keywords.dft.block_size`. Grid points per block; -1 keeps the
         !! backend default. Sets both how tightly the screen bites and how many
         !! blocks there are to spread over threads.
      integer :: dft_grid_level = 3
      integer :: dft_nlc_grid_level = -1
         !! `keywords.dft.nlc_grid_level`, the quadrature VV10's double integral
         !! uses. Separate from `grid_level` because the non-local term costs
         !! the *product* of two grid sizes while everything else costs one, so
         !! the level that is right for exchange is an order of magnitude too
         !! expensive here. Negative means the backend's own default.
      integer :: dft_radial_points = -1    !! -1 leaves grid_level in charge
      integer :: dft_angular_points = -1

      ! keywords.pcm -- the polarizable continuum on the ab initio backends.
      !
      ! Not the xTB solvation keys below: those configure tblite's own CPCM,
      ! which builds its own cavity with its own defaults -- cpcm_rscale is 1.0
      ! there, against the 1.2 conventional for a van der Waals cavity.
      ! Which integral backend to run on: "auto" (default), "cuest"/"gpu", or
      ! "libcint"/"cpu". A request that the build or the method cannot honour is
      ! refused, not substituted.
      character(len=16) :: backend = "auto"
      ! system.gpu -- the same choice said the way people ask for it. It
      ! resolves into `backend` above, so there is one answer downstream however
      ! the deck spelled the question; naming both and disagreeing is refused.
      logical :: gpu = .false.
      logical :: gpu_set = .false.
         !! Whether the deck named `system.gpu`. "Absent" and "false" have to be
         !! told apart: absent leaves `backend` alone, false pins it to the CPU.

      logical :: pcm_enabled = .false.
      character(len=16) :: pcm_method = "cpcm"
         !! Continuum model: "cpcm" or "iefpcm". See `pcm_config_t`.
      real(dp) :: pcm_dielectric = -1.0_dp
         !! Solvent dielectric constant. Required: there is no solvent-name
         !! table on this path.
      integer :: pcm_angular_points = DEFAULT_PCM_NANG
      real(dp) :: pcm_radii_scale = DEFAULT_PCM_RSCALE
      real(dp) :: pcm_zeta = DEFAULT_PCM_ZETA
      real(dp) :: pcm_tolerance = 1.0e-8_dp
      integer :: pcm_max_iter = 100

      ! XTB solvation settings
      character(len=:), allocatable :: solvent  !! Solvent name (e.g., "water", "ethanol") or empty for gas phase
      character(len=:), allocatable :: solvation_model  !! Solvation model: "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.               !! Include non-polar CDS terms in solvation (not for CPCM)
      logical :: use_shift = .true.             !! Include solution state shift in solvation (not for CPCM)
      ! CPCM-specific settings
      real(dp) :: dielectric = -1.0_dp              !! Direct dielectric constant (-1 = use solvent lookup)
      integer :: cpcm_nang = DEFAULT_CPCM_NANG      !! Number of angular grid points for CPCM cavity
      real(dp) :: cpcm_rscale = DEFAULT_CPCM_RSCALE  !! Radii scaling factor for CPCM cavity

      ! Driver information
      integer(int32) :: calc_type = CALC_TYPE_ENERGY

      ! Multiple molecules support
      integer :: nmol = 0  !! Number of molecules (0 = single molecule mode for backward compatibility)
      type(molecule_t), allocatable :: molecules(:)  !! Array of molecules (if nmol > 0)

      ! Single molecule fields (backward compatibility - used if nmol == 0)
      ! Structure information
      integer :: charge = 0
      integer :: multiplicity = 1

      ! Geometry
      type(geometry_type) :: geometry

      ! Fragments
      integer :: nfrag = 0
      type(input_fragment_t), allocatable :: fragments(:)
      !! As on `molecule_t`: every fragment the same species, which is what lets
      !! `fragment_potentials` be one name rather than one per fragment.
      logical :: uniform_system = .false.

      ! Connectivity
      integer :: nbonds = 0
      integer :: nbroken = 0
      type(bond_t), allocatable :: bonds(:)

      ! SCF settings
      ! TODO(mqc): 300 here against `DEFAULT_SCF_MAXITER` = 100 in
      ! `mqc_calculation_defaults`, which `scf_config_t` uses. How many
      ! iterations an SCF is allowed depends on which of the two built its
      ! config, and neither number is derived from the other.
      integer :: scf_maxiter = 300
      real(dp) :: scf_tolerance = DEFAULT_SCF_CONV
      logical :: scf_maxiter_set = .false.
         !! Whether the deck named `keywords.scf.maxiter`; see `scf_config_t`.
      character(len=:), allocatable :: scf_convergence_metric
         !! `keywords.scf.convergence_metric`. Unallocated leaves the default,
         !! which is today's energy-and-commutator pair.
      logical :: scf_tolerance_set = .false.
         !! Whether the deck named it. "The default" and "the user asked for
         !! the default" are different requests to a caller with a stricter
         !! default of its own -- MAKEFP converges to 1e-10/1e-8, because the
         !! multipoles are taken from that density, unless a deck says otherwise.
      real(dp) :: scf_density_tolerance = DEFAULT_SCF_DENSITY_CONV
      real(dp) :: scf_gradient_tolerance = 0.0_dp
         !! `keywords.scf.gradient_tolerance`; zero derives `sqrt(tolerance)`
      real(dp) :: scf_linear_dependence = 0.0_dp
         !! `keywords.scf.linear_dependence_threshold`. Overlap eigenvalues at
         !! or below this are dropped as linearly dependent. Raise it to shed
         !! more of a diffuse basis, lower it to keep modes the default would
         !! discard.
         !!
         !! Zero means "not set", and the orthogonaliser then uses
         !! `LINEAR_DEPENDENCE_TOL` from `mqc_scf_common`.
      logical :: scf_density_tolerance_set = .false.
      real(dp) :: scf_level_shift = 0.0_dp
         !! Hartree added to the virtual orbitals before each diagonalisation.
         !! Zero is off, which is the default: a shift costs iterations near the
         !! solution and only earns them far from it.
      logical :: scf_use_diis = .true.
         !! `keywords.scf.diis`. Off is a diagnostic rather than a setting: an
         !! SCF without DIIS is how you find out whether DIIS is the thing
         !! hiding a problem, not a way to run one.
      integer :: scf_diis_size = 8
         !! `keywords.scf.diis_size`, the subspace the extrapolation is drawn
         !! from. Widening it to 12-20 is the first thing to try on an
         !! open-shell SCF that converges monotonically but slowly, where a
         !! level shift would only make it slower.

      ! EFP settings, read from `keywords.efp` and used only by the MAKEFP
      ! driver. These are the stages after the SCF; the SCF a potential runs is
      ! configured by `keywords.scf` and has no second spelling here.
      real(dp) :: efp_dynamic_tolerance = DEFAULT_DYNAMIC_TOL
      integer :: efp_dynamic_maxiter = DEFAULT_DYNAMIC_MAXITER
      logical :: efp_allow_crap_response = .false.
      integer :: efp_response_batch = 0
         !! `keywords.efp.response_batch`: densities per integral pass in the
         !! matrix-free response. Zero, the default, is as wide as the memory
         !! budget allows, up to every system at once; a positive figure is
         !! taken as given.
      integer :: efp_response = EFP_RESPONSE_AUTO
         !! Held as a code rather than as the spelling, like `calc_type`; the
         !! reader is the only place that knows the three words it can take.
      real(dp) :: efp_vdw_scale = DEFAULT_VDW_SCALE
      logical :: efp_quadrupole_blocks = .true.
         !! `keywords.efp.dispersion`: "all" writes the dipole-quadrupole and
         !! quadrupole-quadrupole dynamic blocks as well as the dipole-dipole
         !! one, which is what GAMESS does by default (`DISP7`, `DISP8`);
         !! "dipole" writes the dipole-dipole block alone and solves the
         !! response for the three dipole operators only.

      ! Quantum nuclei, read from `keywords.neo`. Absent means every nucleus is
      ! classical and nothing below is looked at.
      logical :: neo_active = .false.
      character(len=64) :: neo_nuclear_basis = "pb4-d"
      integer, allocatable :: neo_quantum_indices(:)   !! 1-based, from a 0-based deck
      character(len=8), allocatable :: neo_quantum_symbols(:)
      character(len=8) :: neo_epc = ""   !! "17-1" or "17-2"; empty is none

      ! Hessian settings
      real(dp) :: hessian_displacement = DEFAULT_DISPLACEMENT  !! Finite difference displacement (Bohr)
      real(dp) :: hessian_temperature = DEFAULT_TEMPERATURE    !! Temperature for thermochemistry (K)
      real(dp) :: hessian_pressure = DEFAULT_PRESSURE          !! Pressure for thermochemistry (atm)
      real(dp) :: hessian_response_tol = DEFAULT_RESPONSE_TOL
         !! Residual at which the analytic Hessian's coupled-perturbed solve stops
      integer :: hessian_response_max_iter = DEFAULT_RESPONSE_MAX_ITER
         !! Cycles that solve may take
      integer :: hessian_response_batch = DEFAULT_RESPONSE_BATCH
         !! Perturbations solved together, sharing each integral pass

      ! AIMD settings
      real(dp) :: aimd_dt = DEFAULT_AIMD_DT                          !! Timestep (femtoseconds)
      integer :: aimd_nsteps = DEFAULT_AIMD_NSTEPS                   !! Number of MD steps (0 = no AIMD)
      real(dp) :: aimd_initial_temperature = DEFAULT_AIMD_TEMPERATURE  !! Initial temperature for velocity init (K)
      integer :: aimd_output_frequency = DEFAULT_AIMD_OUTPUT_FREQ    !! Write output every N steps

      ! Geometry optimization settings. The spelled ones stay strings here and
      ! are parsed in the adapter, where a misspelling is refused by name.
      integer :: opt_max_steps = DEFAULT_OPT_MAX_STEPS
      real(dp) :: opt_gradient_tolerance = DEFAULT_OPT_GRADIENT_TOLERANCE  !! Hartree/Bohr
      real(dp) :: opt_energy_tolerance = DEFAULT_OPT_ENERGY_TOLERANCE      !! Hartree; <0 = engine default
      real(dp) :: opt_max_step = DEFAULT_OPT_MAX_STEP                      !! Bohr; <0 = engine default
      integer :: opt_lbfgs_memory = DEFAULT_OPT_LBFGS_MEMORY               !! <0 = engine default
      integer :: opt_print_level = DEFAULT_OPT_PRINT_LEVEL                 !! <0 = engine default
      character(len=:), allocatable :: opt_coordinates  !! cartesian, hdlc, dlc
      character(len=:), allocatable :: opt_algorithm    !! lbfgs, cg, sd, prfo
      logical :: opt_trajectory = .true.               !! Record the path taken
      logical :: opt_freeze_terms = .true.             !! Fix the MBE term list for the run
      logical :: opt_hess_end = .false.                !! Hessian at the converged geometry
      character(len=:), allocatable :: opt_hessian_update  !! none, powell, bofill, auto
      character(len=:), allocatable :: opt_target  !! minimum or saddle
      character(len=:), allocatable :: opt_endpoint  !! xyz of the product, for NEB
      character(len=:), allocatable :: opt_neb_ends  !! frozen, perpendicular, free
      integer :: opt_n_images = 0                      !! Images along the path
      real(dp) :: opt_neb_spring = -1.0_dp             !! <0 = engine default
      character(len=:), allocatable :: opt_saddle_method  !! prfo or dimer
      real(dp) :: opt_dimer_separation = -1.0_dp       !! <0 = engine default
      integer :: opt_dimer_max_rotations = -1          !! <0 = engine default
      real(dp) :: opt_dimer_rot_tol = -1.0_dp          !! <0 = engine default
      integer :: opt_zero_modes = -1                   !! <0 = engine default
      real(dp) :: opt_soft_modes = -1.0_dp             !! <0 = engine default
      logical :: opt_connect = .false.                 !! Relax downhill from the saddle
      real(dp) :: opt_connect_distort = -1.0_dp        !! <0 = engine default
      real(dp) :: opt_timestep = -1.0_dp               !! Damped dynamics step
      real(dp) :: opt_friction = -1.0_dp               !! Damped dynamics start friction
      real(dp) :: opt_friction_factor = -1.0_dp        !! Friction decay while descending
      real(dp) :: opt_friction_rising = -1.0_dp        !! Friction on an uphill step
      integer, allocatable :: opt_frozen_atoms(:)      !! 0-based in the deck, 1-based here
      integer, allocatable :: opt_constraint_kinds(:)  !! One per constraint
      integer, allocatable :: opt_constraint_atoms(:, :)
         !! (4, n_constraints), 1-based, unused slots zero

      ! Fragmentation settings
      character(len=:), allocatable :: frag_method  !! MBE, etc.
      integer :: frag_level = DEFAULT_FRAG_LEVEL
         !! Which expansion runs. `keywords.fragmentation.method` is the only
         !! input: `allow_overlapping_fragments` and `expansion` were two more
         !! ways to say the same thing, and having three meant the resolver had
         !! to reconcile them -- which it got wrong for an explicit
         !! `{"method": "gmbe", "allow_overlapping_fragments": false}`,
         !! silently flipping the flag back. Deleting them removes the class
         !! rather than fixing the instance. The schema refuses unknown keys,
         !! so a deck still naming them fails at parse time and says which.
      character(len=:), allocatable :: counterpoise
         !! How basis-set superposition error is handled: "none" (default) or
         !! "vmfc".
         !!
         !! Under "none" each subfragment is solved in its own basis, so a pair
         !! borrows functions its monomers did not have and looks more bound
         !! than it is; that error lands in the pair term and survives
         !! truncation. "vmfc" solves every subfragment in its parent's basis
         !! instead -- valiron-mayer function counterpoise -- so the borrowing
         !! appears on both sides of each difference and cancels.
         !!
         !! It is not free: a monomer's energy becomes one number per pair
         !! rather than one reusable across every pair it belongs to, so a
         !! level-2 expansion goes from N + C(N,2) subcalculations to 3*C(N,2).

      character(len=:), allocatable :: fmo_far_field
         !! What a distant fragment contributes: mulliken, chelpg or ignore
      real(dp) :: fmo_resppc = 2.0_dp
         !! Separation past which a fragment becomes point charges; negative
         !! turns the approximation off and makes every neighbour exact
      integer :: fmo_max_outer = 50
      real(dp) :: fmo_tolerance = 1.0e-7_dp
         !! Outer (monomer) SCF convergence, on the monomer energy sum
      integer :: fmo_scf_max_iter = 100
      real(dp) :: fmo_scf_energy_tol = 1.0e-9_dp
      real(dp) :: fmo_scf_density_tol = 1.0e-7_dp
         !! The inner per-fragment SCF: iteration cap and the energy/density
         !! convergence each monomer and n-mer SCF is held to. Independent of the
         !! outer loop above, and of a top-level `keywords.scf`, so a fragment run
         !! can be converged more loosely than a whole-system one would be.
      real(dp) :: efmo_rcut = 2.0_dp
         !! `keywords.fragmentation.rcut`: where EFMO stops solving a dimer.
         !!
         !! **Unitless, not an Angstrom distance.** Each interatomic separation
         !! is divided by the two van der Waals radii, so 1 is contact and the
         !! 2.0 default is twice that -- a threshold that means the same thing
         !! for two waters and for two aromatic rings, which
         !! `keywords.fragmentation.cutoffs` (Angstrom, and MBE's) does not.
         !! It is the `resppc` of FMO measured the same way, on a different
         !! decision: a pair inside it is a quantum dimer, a pair beyond it is
         !! four effective-fragment terms.
      logical :: efmo_charge_transfer = .true.
         !! `keywords.efmo.charge_transfer`: include `E_IJ^CT` in the far pairs.
         !! GAMESS's EFMO has it; the 2012 method left it out.
      integer :: max_intersection_level = DEFAULT_MAX_INTERSECTION  !! Maximum k-way intersection depth for GMBE
      character(len=:), allocatable :: bond_breaking
         !! How a fragment represents a covalent bond the partition cut.
         !!
         !! `"none"` refuses a partition that cuts one. `"caps"` closes it with
         !! a hydrogen. `"afo"` is the adjusted frozen orbital -- see
         !! `mqc_czt_afo.f90` and `mqc_docs/source/fmo.rst`. An FMO
         !! expansion accepts `"none"` and `"afo"`; it refuses `"caps"` as not
         !! implemented for that expansion.
         !!
         !! Independent of `embedding`, but not every pairing of the two is
         !! sound: a cap puts an electron where an exact embedding already
         !! supplies the neighbour's density, so that combination is refused.
      real(dp) :: cap_scale = 1.0_dp
         !! Where the cap goes along the cut bond: `R_H = R_kept + s (R_gone - R_kept)`.
         !!
         !! `1.0`, the default, puts it exactly where the removed atom was, so
         !! redistributing the cap's gradient onto that atom is the whole chain
         !! rule rather than an approximation. A physical C-H wants roughly
         !! `0.71` of a C-C, and the gradient then splits `s` and `1 - s`
         !! between the two atoms -- still exact, one term longer.
      character(len=:), allocatable :: embedding
      character(len=:), allocatable :: cutoff_method
      character(len=:), allocatable :: distance_metric
      real(dp), allocatable :: fragment_cutoffs(:)  !! Distance cutoffs indexed by n-mer level (2=dimer, 3=trimer, etc.)
      integer :: global_groups = 0
      integer :: nodes_per_group = 0

      ! Logger settings (kept for compatibility)
      character(len=:), allocatable :: log_level

      real(dp) :: memory_gb = -1.0_dp
         !! `system.memory_gb`: what this run may plan on, in gigabytes, for
         !! the decisions that trade memory for time -- which route builds the
         !! response Hessian, whether an operator is formed at all, how many
         !! frequencies factorize at once. Negative, the default, means the
         !! machine decides: a share of what it reports available.

      ! Output control
      logical :: skip_json_output = .false.  !! Skip JSON output for large calculations
      character(len=:), allocatable :: checkpoint_file
         !! Where to append each fragment energy as it finishes, and where to
         !! resume from. Absent means neither. An existing file whose
         !! fingerprint disagrees with this run is refused, not overwritten.
      logical :: unchecked_input = .false.
         !! Downgrade semantic validation to warnings.
         !!
         !! Off by default. The checks it disables are for inputs that produce a
         !! converged, plausible, wrong answer -- monomer charges that do not
         !! sum, atoms in no fragment, a term list the expansion cannot be
         !! assembled from.
      character(len=:), allocatable :: fragment_breakdown
         !! Where the per-fragment MBE table goes: "csv", "json" or "none"

   contains
      procedure :: destroy => config_destroy
   end type mqc_config_t

contains

   subroutine molecule_destroy(this)
      !! Clean up allocated memory in molecule_t
      class(molecule_t), intent(inout) :: this
      integer :: i

      if (allocated(this%name)) deallocate (this%name)

      call this%geometry%destroy()

      if (allocated(this%fragments)) then
         do i = 1, size(this%fragments)
            call this%fragments(i)%destroy()
         end do
         deallocate (this%fragments)
      end if

      if (allocated(this%bonds)) deallocate (this%bonds)

   end subroutine molecule_destroy

   subroutine config_destroy(this)
      !! Clean up allocated memory in mqc_config_t
      class(mqc_config_t), intent(inout) :: this
      integer :: i

      if (allocated(this%schema_name)) deallocate (this%schema_name)
      if (allocated(this%schema_version)) deallocate (this%schema_version)
      if (allocated(this%units)) deallocate (this%units)
      if (allocated(this%basis)) deallocate (this%basis)
      if (allocated(this%aux_basis)) deallocate (this%aux_basis)
      if (allocated(this%functional)) deallocate (this%functional)
      if (allocated(this%charges_scheme)) deallocate (this%charges_scheme)
      if (allocated(this%log_level)) deallocate (this%log_level)
      if (allocated(this%fragment_breakdown)) deallocate (this%fragment_breakdown)
      if (allocated(this%frag_method)) deallocate (this%frag_method)
      if (allocated(this%embedding)) deallocate (this%embedding)
      if (allocated(this%bond_breaking)) deallocate (this%bond_breaking)
      if (allocated(this%cutoff_method)) deallocate (this%cutoff_method)
      if (allocated(this%distance_metric)) deallocate (this%distance_metric)
      if (allocated(this%fragment_cutoffs)) deallocate (this%fragment_cutoffs)

      call this%geometry%destroy()

      if (allocated(this%fragments)) then
         do i = 1, size(this%fragments)
            call this%fragments(i)%destroy()
         end do
         deallocate (this%fragments)
      end if

      if (allocated(this%bonds)) deallocate (this%bonds)

      ! Clean up molecules array (multi-molecule mode)
      if (allocated(this%molecules)) then
         do i = 1, size(this%molecules)
            call this%molecules(i)%destroy()
         end do
         deallocate (this%molecules)
      end if

   end subroutine config_destroy

   subroutine input_fragment_destroy(this)
      !! Clean up allocated memory in input_fragment_t
      class(input_fragment_t), intent(inout) :: this
      if (allocated(this%indices)) deallocate (this%indices)
   end subroutine input_fragment_destroy

end module mqc_config_types
