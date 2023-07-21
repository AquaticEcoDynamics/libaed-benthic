!###############################################################################
!#                                                                             #
!# aed_habitat_galaxiid.F90                                                    #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2021 - 2022 -  The University of Western Australia               #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2021                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_habitat_galaxiid
!-------------------------------------------------------------------------------
! aed_habitat_glaxiid --- galaxiid habitat model
!
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_habitat_galaxiid_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_habitat_galaxiid_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_bird, id_mtox
      INTEGER :: id_chsi, id_chpl, id_chfl, id_chsd
      INTEGER :: id_fssi, id_fsdp, id_fssd, id_fsst, id_fsls, id_fsdw, id_fsmt
      INTEGER :: id_rhsi, id_rhpl, id_rhfl, id_rhsd, id_rhtr, id_rhsp, id_rhtd
      INTEGER :: id_chhsi, id_chhpl, id_chhfl, id_chhsd, id_chhtr, id_chhsp
      INTEGER :: id_crhsi, id_crhpl, id_crhfl, id_crhsd, id_crhtr
      INTEGER :: id_pshsi, id_pshsd, id_pshpl, id_pshfl
      INTEGER :: id_wettime, id_drytime
      INTEGER, ALLOCATABLE :: id_d_rupfs(:),id_d_rupft(:),id_d_rupfl(:),id_d_rupfa(:),id_d_rupfd(:)
      INTEGER, ALLOCATABLE :: id_d_chafs(:),id_d_chaft(:),id_d_chafl(:),id_d_chafa(:),id_d_chafd(:),id_d_chafv(:)
      INTEGER, ALLOCATABLE :: id_d_crcfm(:),id_d_crcft(:),id_d_crcfv(:),id_d_crcfs(:),id_d_crcfd(:),id_d_crcfp(:)
      INTEGER, ALLOCATABLE :: id_d_pssfm(:), id_d_pssft(:),id_d_pssfv(:),id_d_pssfs(:),id_d_pssfd(:), id_d_pssfu(:)

      !# Dependencies
      INTEGER :: id_l_ph, id_l_hab, id_l_aass, id_l_rveg, id_l_bveg
      INTEGER :: id_l_salg, id_l_falg, id_d_turb, id_l_ncs1, id_l_ncs2, id_l_tau0
      INTEGER :: id_l_otrc, id_l_oxy, id_l_sav
      INTEGER :: id_l_svwc, id_l_stmp25, id_l_stmp, id_l_veg1, id_l_veg2, id_l_pass
      INTEGER, ALLOCATABLE :: id_l_mtox(:)

      !# Environment variables
      INTEGER :: id_E_temp, id_E_salt, id_E_bathy, id_E_matz, id_E_depth
      INTEGER :: id_E_nearlevel, id_E_extc, id_E_Io, id_E_stress, id_E_airtemp

      !# Model switches
      LOGICAL :: simBirdForaging,simBenthicProd,simFishTolerance,simGalaxiidSpawning
      LOGICAL :: simCrabHabitat,simCharaHabitat
      LOGICAL :: simMosquitoRisk,simCyanoRisk
      LOGICAL :: simMetalTox,simClearWater
      LOGICAL :: simCrocEggs,simPassiflora
      INTEGER :: simRuppiaHabitat

      !# Model parameters
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      INTEGER :: num_mtox, n_zones_chara, n_zones_fishspawn
      INTEGER :: n_zones_crocs, n_zones_pass
      AED_REAL,ALLOCATABLE :: active_zones_chara(:), active_zones_fishspawn(:)
      AED_REAL,ALLOCATABLE :: active_zones_crocs(:), active_zones_pass(:)


     CONTAINS
         PROCEDURE :: define             => aed_define_habitat_galaxiid
!        PROCEDURE :: calculate          => aed_calculate_habitat_galaxiid
!        PROCEDURE :: calculate_benthic  => aed_calculate_benthic_habitat_galaxiid
         PROCEDURE :: calculate_riparian => aed_calculate_riparian_habitat_galaxiid
!        PROCEDURE :: mobility           => aed_mobility_habitat_galaxiid
!        PROCEDURE :: light_extinction   => aed_light_extinction_habitat_galaxiid
!        PROCEDURE :: delete             => aed_delete_habitat_galaxiid

   END TYPE

!-------------------------------------------------------------------------------
!MODULE VARIABLES
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep
   LOGICAL :: extra_diag = .false.
   INTEGER :: diag_level = 10                ! 0 = no diagnostic outputs
                                             ! 1 = basic diagnostic outputs
                                             ! 2 = flux rates, and supporitng
                                             ! 3 = other metrics
                                             !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_habitat_galaxiid(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED HABITAT module
!
!  Here, the aed namelist is read and the variables exported
!  are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_habitat_galaxiid_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, z, status, num_mtox

!  %% NAMELIST   %%  /aed_habitat_galaxiid/
   INTEGER           :: n_zones_chara = 0
   INTEGER           :: active_zones_chara(MAX_ZONES)
   INTEGER           :: n_zones_fishspawn = 0
   INTEGER           :: active_zones_fishspawn(MAX_ZONES)
   LOGICAL           :: simFishTolerance = .FALSE.
   LOGICAL           :: simGalaxiidSpawning = .FALSE.
   LOGICAL           :: simBenthicProd = .FALSE.
   LOGICAL           :: simCyanoRisk = .FALSE.
   LOGICAL           :: simMosquitoRisk = .FALSE.
   LOGICAL           :: simCrabHabitat = .FALSE.
   LOGICAL           :: simClearWater = .FALSE.
   LOGICAL           :: simCharaHabitat = .FALSE.
   LOGICAL           :: simMetalTox = .FALSE.
   INTEGER           :: simRuppiaHabitat = 0
   AED_REAL          :: mtox_lims(10)
   CHARACTER(len=40) :: mtox_vars(10)

! %% From Module Globals
!  LOGICAL :: extra_diag = .false.      !## Obsolete Use diag_level = 10
!  INTEGER :: diag_level = 10
!  %% END NAMELIST   %%  /aed_habitat_galaxiid/

   CHARACTER(len=64) :: bird_acid_link, bird_habs_link, bird_aass_link, bird_rveg_link, bird_bveg_link
   CHARACTER(len=64) :: fshsi_veg_link, fshsi_oxy_link, fshsi_otrc_link
   CHARACTER(len=64) :: chsi_otrc_link, chsi_oxy_link, chsi_veg_link
   CHARACTER(len=64) :: chhsi_salg_link,chhsi_falg_link
   CHARACTER(len=64) :: chhsi_ncs1_link,chhsi_ncs2_link,chhsi_tau0_link
   CHARACTER(len=64) :: crhsi_ncs1_link,crhsi_ncs2_link
   CHARACTER(len=64) :: crhsi_stmp_link,crhsi_svwc_link,crhsi_pass_link
   CHARACTER(len=64) :: pshsi_stmp_link,pshsi_svwc_link,pshsi_veg1_link,pshsi_veg2_link
   CHARACTER(len=64) :: pshsi_ncs1_link, pshsi_ncs2_link
   CHARACTER(len=64) :: rhsi_salg_link, rhsi_falg_link
   CHARACTER(len=40) :: mtox_acid_link, mtox_aass_link

   NAMELIST /aed_habitat_galaxiid/ &
                           simFishTolerance, &
                           simGalaxiidSpawning,n_zones_fishspawn,active_zones_fishspawn,  &
                           simBenthicProd,   &
                           simCyanoRisk,     &
                           simMosquitoRisk,  &
                           simCrabHabitat,   &
                           simRuppiaHabitat, &
                           simCharaHabitat,n_zones_chara,active_zones_chara, &
                           simClearWater,    &
                           simMetalTox, mtox_vars, mtox_lims,   &
                           rhsi_falg_link, rhsi_salg_link,      &
                           extra_diag, diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_habitat_galaxiid initialization"
   print *,"          WARNING! aed_habitat model is under development"

   ! Read the namelist
   read(namlst,nml=aed_habitat_galaxiid,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_habitat'

   ! Update module level switches
   data%num_habitats = 0
   data%simFishTolerance = simFishTolerance ; IF(simFishTolerance) data%num_habitats=data%num_habitats+1
   data%simGalaxiidSpawning  = simGalaxiidSpawning  ; IF(simGalaxiidSpawning) data%num_habitats=data%num_habitats+1
   data%simBenthicProd   = simBenthicProd   ; IF(simBenthicProd) data%num_habitats=data%num_habitats+1
   data%simMetalTox      = simMetalTox      ; IF(simMetalTox) data%num_habitats=data%num_habitats+1
   data%simMosquitoRisk  = simMosquitoRisk  ; IF(simMosquitoRisk) data%num_habitats=data%num_habitats+1
   data%simCyanoRisk     = simCyanoRisk     ; IF(simCyanoRisk) data%num_habitats=data%num_habitats+1
   data%simCrabHabitat   = simCrabHabitat   ; IF(simCrabHabitat) data%num_habitats=data%num_habitats+1
   data%simCharaHabitat  = simCharaHabitat  ; IF(simCharaHabitat) data%num_habitats=data%num_habitats+1
   data%simClearWater    = simClearWater    ; IF(simClearWater) data%num_habitats=data%num_habitats+1
   data%simRuppiaHabitat = simRuppiaHabitat ; IF(simRuppiaHabitat>0) data%num_habitats=data%num_habitats+1

   print *,"          ... # habitat templates simulated: ",data%num_habitats

   IF( extra_diag )   diag_level = 10           ! legacy use of extra_debug

   !----------------------------------------------------------------------------
   ! Define variables and dependencies


   !-- CONTAMINATION
   IF( simMetalTox ) THEN
     data%id_mtox =  aed_define_sheet_diag_variable('toxicity','-', 'Suitability')

     mtox_acid_link = 'CAR_pH'
     mtox_aass_link = 'ASS_uzaass'

     mtox_vars = '' ;  mtox_lims = 1.0
     num_mtox = 0
     DO i=1,10 ; IF (mtox_vars(i)  .EQ. '' ) THEN ; num_mtox = i-1 ; EXIT ; ENDIF ; ENDDO
     ALLOCATE(data%id_l_mtox(num_mtox)); ALLOCATE(data%mtox_lims(num_mtox))
     data%num_mtox = num_mtox
     DO i=1,data%num_mtox
       data%id_l_mtox(i) =  aed_locate_variable(mtox_vars(i))
       data%mtox_lims(i) =  mtox_lims(i)
       !print*,'Tox : ', TRIM(tfe_vars(i)), ' * ', data%tfe_varscale(i)
     ENDDO
   ENDIF

   !-- BLUE SWIMMER CRABS
   IF( simCrabHabitat ) THEN
     data%id_chsi =  aed_define_sheet_diag_variable('crab_hsi','-', 'Crab Habitat Suitability Index')
     data%id_chpl =  aed_define_sheet_diag_variable('crab_hsi_larvae','-', 'Crab Habitat Suitability - larval connectivity')
     data%id_chfl =  aed_define_sheet_diag_variable('crab_hsi_wq','-', 'Crab Habitat Suitability - water quality')
     data%id_chsd =  aed_define_sheet_diag_variable('crab_hsi_ben','-', 'Crab Habitat Suitability - benthic habitat')

     chsi_otrc_link = 'TRC_tr1'    ! ocean larvae tracer
     chsi_oxy_link  = 'OXY_oxy'    ! oxygen
     chsi_veg_link  = 'MAC_mac'    ! submerged aquatic vegetation

     data%id_l_otrc  = aed_locate_variable(TRIM(chsi_otrc_link))
     data%id_l_oxy  = aed_locate_variable(TRIM(chsi_oxy_link))
     data%id_l_sav  = aed_locate_sheet_variable(TRIM(chsi_veg_link))
   ENDIF

   !-- SEAGRASS : RUPPIA
   IF( simRuppiaHabitat>0 ) THEN
     data%id_rhsi =  aed_define_sheet_diag_variable('ruppia_hsi','-', 'Ruppia Habitat Suitability Index')
     data%id_rhpl =  aed_define_sheet_diag_variable('ruppia_hsi_plant',  '-', 'Ruppia Habitat Suitability - plant')
     data%id_rhfl =  aed_define_sheet_diag_variable('ruppia_hsi_flower', '-', 'Ruppia Habitat Suitability - flowering')
     data%id_rhsd =  aed_define_sheet_diag_variable('ruppia_hsi_seed',   '-', 'Ruppia Habitat Suitability - seed germination')
     data%id_rhtr =  aed_define_sheet_diag_variable('ruppia_hsi_turion', '-', 'Ruppia Habitat Suitability - turion formation')
     data%id_rhsp =  aed_define_sheet_diag_variable('ruppia_hsi_sprout', '-', 'Ruppia Habitat Suitability - turion sprouting')
     data%id_rhtd =  aed_define_sheet_diag_variable('ruppia_hsi_dormant','-', 'Ruppia Habitat Suitability - turion viability')
    !data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
    !data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')

!    rhsi_falg_link = 'MAG_ulva_ben'
!    rhsi_salg_link = 'MAG_ulva'
     IF (rhsi_falg_link .EQ. "") THEN
         STOP 'need to set rhsi_falg_link and rhsi_salg_link'
     ENDIF

     data%id_l_salg  = aed_locate_variable(TRIM(rhsi_salg_link))
     data%id_l_falg  = aed_locate_sheet_variable(TRIM(rhsi_falg_link))

     IF (diag_level>1) THEN
       ALLOCATE(data%id_d_rupfs(6))
       ALLOCATE(data%id_d_rupft(6))
       ALLOCATE(data%id_d_rupfl(6))
       ALLOCATE(data%id_d_rupfa(6))
       ALLOCATE(data%id_d_rupfd(6))
       DO i =1,6
        data%id_d_rupfs(i) = aed_define_sheet_diag_variable('ruppia_hsi_fsal_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Ruppia Habitat Suitability - fSal')
        data%id_d_rupft(i) = aed_define_sheet_diag_variable('ruppia_hsi_ftem_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Ruppia Habitat Suitability - fTem')
        data%id_d_rupfl(i) = aed_define_sheet_diag_variable('ruppia_hsi_flgt_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Ruppia Habitat Suitability - fLgt')
        data%id_d_rupfa(i) = aed_define_sheet_diag_variable('ruppia_hsi_falg_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Ruppia Habitat Suitability - fAlg')
        data%id_d_rupfd(i) = aed_define_sheet_diag_variable('ruppia_hsi_fdep_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Ruppia Habitat Suitability - fDep')
       ENDDO
     ENDIF

   ENDIF

   !-- CHARAPHYTES
   IF( simCharaHabitat ) THEN

     data%n_zones_chara = n_zones_chara
     IF (n_zones_chara > 0) THEN
        ALLOCATE(data%active_zones_chara(n_zones_chara))
        DO z=1,n_zones_chara
           data%active_zones_chara(z) = active_zones_chara(z)
        ENDDO
     ENDIF

     data%id_chhsi =  aed_define_sheet_diag_variable('chara_hsi','-', 'Chara Habitat Suitability Index')
     data%id_chhpl =  aed_define_sheet_diag_variable('chara_hsi_plant','-', 'Chara Habitat Suitability - plant growth')
     data%id_chhfl =  aed_define_sheet_diag_variable('chara_hsi_oospore','-', 'Chara Habitat Suitability - oospore production')
     data%id_chhsd =  aed_define_sheet_diag_variable('chara_hsi_seed','-', 'Chara Habitat Suitability - seed germination')
     data%id_chhtr =  aed_define_sheet_diag_variable('chara_hsi_turion','-', 'Chara Habitat Suitability - turion production')
     data%id_chhsp =  aed_define_sheet_diag_variable('chara_hsi_sprout','-', 'Chara Habitat Suitability - sprout establishment')

     chhsi_falg_link = 'MAG_cladophora_ben'
     chhsi_salg_link = 'MAG_cladophora'
     data%id_l_salg  = 0 !aed_locate_variable(TRIM(rhsi_salg_link))
     data%id_l_falg  = 0 !aed_locate_sheet_variable(TRIM(rhsi_falg_link))

     chhsi_ncs1_link = 'NCS_fs1'
     chhsi_ncs2_link = 'NCS_fs2'
     chhsi_tau0_link = 'NCS_tau_0'
     data%id_l_ncs1  = aed_locate_sheet_variable(TRIM(chhsi_ncs1_link))
     data%id_l_ncs2  = aed_locate_sheet_variable(TRIM(chhsi_ncs2_link))
     data%id_l_tau0  = aed_locate_sheet_variable(TRIM(chhsi_tau0_link))

     if( diag_level>1 )then
       ALLOCATE(data%id_d_chafs(6))
       ALLOCATE(data%id_d_chaft(6))
       ALLOCATE(data%id_d_chafl(6))
       ALLOCATE(data%id_d_chafa(6))
       ALLOCATE(data%id_d_chafd(6))
       ALLOCATE(data%id_d_chafv(6))
       DO i =1,6
        data%id_d_chafs(i) = aed_define_sheet_diag_variable('chara_hsi_fsal_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Chara Habitat Suitability - fSal')
        data%id_d_chaft(i) = aed_define_sheet_diag_variable('chara_hsi_ftem_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Chara Habitat Suitability - fTem')
        data%id_d_chafl(i) = aed_define_sheet_diag_variable('chara_hsi_flgt_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Chara Habitat Suitability - fLgt')
        data%id_d_chafa(i) = aed_define_sheet_diag_variable('chara_hsi_falg_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Chara Habitat Suitability - fAlg')
        data%id_d_chafd(i) = aed_define_sheet_diag_variable('chara_hsi_fdep_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Chara Habitat Suitability - fDep')
        data%id_d_chafv(i) = aed_define_sheet_diag_variable('chara_hsi_fstr_'//CHAR(ICHAR('0') + i),'-', &
                                                                            'Chara Habitat Suitability - fStr')
       ENDDO
     endif
   ENDIF


   !-- FISH : GALAXIID SPAWNING
   IF( simGalaxiidSpawning ) THEN

     data%n_zones_fishspawn = n_zones_fishspawn
     IF (n_zones_fishspawn > 0) THEN
        ALLOCATE(data%active_zones_fishspawn(n_zones_fishspawn))
        DO z=1,n_zones_fishspawn
           data%active_zones_fishspawn(z) = active_zones_fishspawn(z)
        ENDDO
     ENDIF

     data%id_fssi =  aed_define_sheet_diag_variable('spawning_hsi','-',           &
                                          'Galaxiid Spawning Habitat Suitability Index')

     data%id_fsdp =  aed_define_sheet_diag_variable('spawning_hsi_depth','-',    &
                                          'Galaxiid Spawning Habitat Suitability - depth trigger')
     data%id_fssd =  aed_define_sheet_diag_variable('spawning_hsi_substrate','-',        &
                                          'Galaxiid Spawning Habitat Suitability - substrate quality')
     data%id_fsst =  aed_define_sheet_diag_variable('spawning_hsi_temp','-',       &
                                          'Galaxiid Spawning Habitat Suitability - temperature trigger')
     data%id_fsls =  aed_define_sheet_diag_variable('spawning_hsi_stress','-',    &
                                          'Galaxiid Spawning Habitat Suitability - egg stability')
     data%id_fsdw =  aed_define_sheet_diag_variable('spawning_hsi_dewatering','-',        &
                                          'Galaxiid Spawning Habitat Suitability - egg dewatering')
     data%id_fsmt =  aed_define_sheet_diag_variable('spawning_hsi_maturation','-',       &
                                          'Galaxiid Spawning Habitat Suitability - egg maturation')

     fshsi_otrc_link = 'TOT_turbidity'   ! turbidity
     chhsi_ncs1_link = 'NCS_fs1'
     chhsi_ncs2_link = 'NCS_fs2'
     chhsi_tau0_link = 'NCS_tau_0'
     data%id_d_turb  = aed_locate_variable(TRIM(fshsi_otrc_link))
     data%id_l_ncs1  = aed_locate_sheet_variable(TRIM(chhsi_ncs1_link))
     data%id_l_ncs2  = aed_locate_sheet_variable(TRIM(chhsi_ncs2_link))
     data%id_l_tau0  = aed_locate_sheet_variable(TRIM(chhsi_tau0_link))

   ENDIF

   !-- GENERAL
   IF( simGalaxiidSpawning .OR. simCharaHabitat .OR. simRuppiaHabitat>0) THEN
     data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
     data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')
   ENDIF

   ! Register environmental dependencies
   data%id_E_salt  = aed_locate_global('salinity')
   data%id_E_extc  = aed_locate_global('extc_coef')
   data%id_E_temp  = aed_locate_global('temperature')
   data%id_E_depth = aed_locate_global('layer_ht')
   data%id_E_bathy     = aed_locate_sheet_global('bathy')
   data%id_E_matz      = aed_locate_sheet_global('material')
   data%id_E_Io        = aed_locate_sheet_global('par_sf')
   data%id_E_airtemp   = aed_locate_sheet_global('air_temp')
   data%id_E_stress    = aed_locate_sheet_global('taub')
  !data%id_E_nearlevel = aed_locate_sheet_global('nearest_depth')
END SUBROUTINE aed_define_habitat_galaxiid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian_habitat_galaxiid(data,column,layer_idx,pc_wet)
!-------------------------------------------------------------------------------
! Calculate galaxiid habitat calculations
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_galaxiid_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wlevel, extc, bathy, matz, Io, Ig, vel, stress, tau0, stem05, stem25

   ! State
   AED_REAL :: depth, ph, hab, sdepth, uzaass, aass, conc, mtox, turb, grav, stem, svwc

   ! Temporary variables
   INTEGER  :: i
   AED_REAL :: bird_acid, bird_soil, bird_rveg, bird_bveg, bird_salt, bird_dept, bird_habs
   AED_REAL :: euphotic, drytime, light

   ! Parameters
   AED_REAL, PARAMETER :: crit_soil_acidity = 100.000
   AED_REAL, PARAMETER :: crit_water_ph = 6.0
   AED_REAL, PARAMETER :: crit_soil_type = 5.5 ! (matz 6&7 is sand)
   AED_REAL, PARAMETER :: crit_rveg_depth = 0.6
   AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
   AED_REAL, PARAMETER :: crit_salinity = 50.0
   AED_REAL, PARAMETER :: crit_leg_depth = 0.12
   AED_REAL, PARAMETER :: crit_hab_conc = 500.

   AED_REAL :: fs_sdepth , fs_substr, fs_spntem, fs_stress, fs_dewatr, fs_mattem

   AED_REAL :: rhpl=0.,rhfl=0.,rhsd=0.,rhtr=0.,rhsp =0.,falg,rhtd=1.0
   AED_REAL :: pshpl, pshfl, pshsd, pass, height
   AED_REAL :: crns = 0.,creg = 0.,crht = 0.,crml = 0.
   AED_REAL :: limitation(6,6)

!-------------------------------------------------------------------------------
!BEGIN
   matz = 0.0 ; salt = 0.0 ; euphotic = 0.0 ; bathy = 0.0  !## CAB [-Wmaybe-uninitialized]
   limitation = 0.                                         !## CAB [-Wmaybe-uninitialized]


      !---------------------------------------------------------------------------+
      !-- HABITAT TEMPLATE 4: Ruppia Habitat
      IF( data%simRuppiaHabitat>0 ) THEN

        depth = _STATE_VAR_(data%id_E_depth)  ! metres
        salt  = _STATE_VAR_(data%id_E_salt)   ! salinity g/L
        temp  = _STATE_VAR_(data%id_E_temp)   ! degC
        extc  = _STATE_VAR_(data%id_E_extc)   ! /m
        Io    = _STATE_VAR_S_(data%id_E_Io)   ! W/m2
        falg  =(_STATE_VAR_S_(data%id_l_falg) +   &
                _STATE_VAR_(data%id_l_salg)*depth) * 12. * 1e-3/0.5  ! convert mmolC/m2 to gDW/m2
        vel   = 0. !

        ! Inundation time counter and wetness checker
        IF( pc_wet < 0.1 ) THEN
          _DIAG_VAR_S_(data%id_drytime) = _DIAG_VAR_S_(data%id_drytime) + DDT
          IF( _DIAG_VAR_S_(data%id_drytime)>2. ) _DIAG_VAR_S_(data%id_wettime) = zero_
        ELSE
          _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT
          IF( _DIAG_VAR_S_(data%id_wettime)>2. )_DIAG_VAR_S_(data%id_drytime) = zero_
        ENDIF

        ! Overall HSI : Habitat Suitability Index (Issue here re time integration)
        _DIAG_VAR_S_(data%id_rhsi) = (rhpl+rhfl+rhsd+rhtr+rhsp)/5.

        iF( diag_level>1 ) THEN
          DO i = 1,6
           _DIAG_VAR_S_(data%id_d_rupfs(i)) = limitation(i,1)
           _DIAG_VAR_S_(data%id_d_rupft(i)) = limitation(i,2)
           _DIAG_VAR_S_(data%id_d_rupfl(i)) = limitation(i,3)
           _DIAG_VAR_S_(data%id_d_rupfa(i)) = limitation(i,4)
           _DIAG_VAR_S_(data%id_d_rupfd(i)) = limitation(i,5)
         !_DIAG_VAR_S_(data%id_d_rupfm(i)) = limitation(i,6)
         ENDDO
        ENDiF
      ENDIF

END SUBROUTINE aed_calculate_riparian_habitat_galaxiid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_habitat_galaxiid
