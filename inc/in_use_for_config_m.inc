!STARTOFREGISTRYGENERATEDINCLUDE 'inc/in_use_for_config_m.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF(TRIM(vname).EQ.'mh_urb2d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_urban_physics(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_urban_physics(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mosaic_cat_index')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%sf_surface_mosaic.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%sf_surface_mosaic.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'massflux_edkf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mf_temf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bl_pbl_physics(id).EQ.10
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bl_pbl_physics(id).EQ.10
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'md3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mu3d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'maxg2d')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.7
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.7
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mfup_cup')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.10
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.10
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mfup_ent_cup')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.10
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.10
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mfdn_cup')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.10
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.10
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mfdn_ent_cup')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%cu_physics(id).EQ.10
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%cu_physics(id).EQ.10
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'max_p')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%nest_interp_coord.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%nest_interp_coord.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'min_p')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%nest_interp_coord.EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%nest_interp_coord.EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_isop')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_apin')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_bcar')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_acet')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_mbo')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_no')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'msebio_isop')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mlai')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mtsa')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mswdown')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_isop')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_apin')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_bcar')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_acet')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_mbo')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mebio_no')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'msebio_isop')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mlai')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mtsa')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mswdown')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%bio_emiss_opt(id).EQ.4
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%bio_emiss_opt(id).EQ.4
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agtf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agsv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_aggr')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agtf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agsv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_aggr')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.2
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.2
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agtf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agsv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_aggr')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.3
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.3
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agtf')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agef')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_agsv')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'mean_fct_aggr')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%biomass_burn_opt(id).EQ.5
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%biomass_burn_opt(id).EQ.5
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'midrh_min')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%afwa_severe_opt(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%afwa_severe_opt(id).EQ.1
  ENDIF
ENDIF
IF(TRIM(vname).EQ.'midrh_min_old')THEN
  IF(uses.EQ.0)THEN
    in_use = model_config_rec%afwa_severe_opt(id).EQ.1
    uses = 1
  ELSE
    in_use = in_use.OR.model_config_rec%afwa_severe_opt(id).EQ.1
  ENDIF
ENDIF
!ENDOFREGISTRYGENERATEDINCLUDE
