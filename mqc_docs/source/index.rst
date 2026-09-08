.. metalquicha documentation master file, created by
   sphinx-quickstart on Mon Dec 29 15:08:50 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================
User guide to Metalquicha
=========================

This is the user guide to Metalquicha, a Fortran program for quantum chemistry
calculations. It is designed to be modular and extensible, allowing users to
easily add new features and methods.

The API docs for the code itself can be found here: https://jorgeg94.github.io/metalquicha/

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   history
   installation
   building_on_clusters
   getting_started
   capabilities
   input_files
   python_interface
   json_output
   validation
   scf_guess
   scf_convergence
   vibrational_analysis
   analytic_hessians
   geometry_optimization
   conformer_sampling
   fmo
   efmo
   counterpoise
   continuum_solvation
   charges_and_bond_orders
   bonding_analysis
   makefp
   sapt
   neo

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide:

   developer_input_parameters
   developer_json_output
   developer_method_config
