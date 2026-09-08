target_sources(
  ${main_lib}
  PRIVATE mqc_czt_efp_potential.f90
          mqc_czt_efp_read.f90
          mqc_czt_efp_interaction.f90
          mqc_czt_efp_pair.f90
          mqc_czt_efp_rotate.f90
          mqc_czt_efp_energy.f90
          mqc_czt_efp_convert.f90)
