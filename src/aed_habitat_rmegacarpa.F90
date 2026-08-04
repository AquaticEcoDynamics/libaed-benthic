!###############################################################################
!#                                                                             #
!# aed_habitat_rmegacarpa.F90                                               #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2021-2026 : The University of Western Australia                  #
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
!# Created Aug 2026 by Sherry Zhai                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_habitat_rmegacarpa
!-------------------------------------------------------------------------------
! aed_habitat_rmegacarpa --- ruppia megacarpa habitat model
! (adult, flower, seed germination only)
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_habitat_rmegacarpa_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_habitat_rmegacarpa_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_mtox
      INTEGER :: id_rhsi, id_rhpl, id_rhfl, id_rhsd
      INTEGER :: id_wettime, id_drytime
      INTEGER, ALLOCATABLE :: id_d_rupfs(:),id_d_rupft(:),id_d_rupfl(:),id_d_rupfa(:),id_d_rupfd(:)

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
      !     LOGICAL :: simBirdForaging,simBenthicProd,simFishTolerance,simGalaxiidSpawning
!     LOGICAL :: simCrabHabitat,simCharaHabitat
!     LOGICAL :: simMosquitoRisk,simCyanoRisk
!     LOGICAL :: simMetalTox,simClearWater
!     LOGICAL :: simCrocEggs,simPassiflora
      INTEGER :: simRuppiaHabitat

      !# Model parameters
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      INTEGER :: num_mtox


     CONTAINS
         PROCEDURE :: define             => aed_define_habitat_rmegacarpa
!        PROCEDURE :: calculate          => aed_calculate_habitat_rmegacarpa
!        PROCEDURE :: calculate_rmegacarpa  => aed_calculate_benthic_habitat_rmegacarpa
         PROCEDURE :: calculate_riparian => aed_calculate_riparian_habitat_rmegacarpa
!        PROCEDURE :: mobility           => aed_mobility_habitat_rmegacarpa
!        PROCEDURE :: light_extinction   => aed_light_extinction_habitat_rmegacarpa
!        PROCEDURE :: delete             => aed_delete_habitat_rmegacarpa 

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
SUBROUTINE aed_define_habitat_rmegacarpa(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED HABITAT module
!
!  Here, the aed namelist is read and the variables exported
!  are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_habitat_rmegacarpa_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, z, status, num_mtox

!  %% NAMELIST   %%  /aed_habitat_rmegacarpa/
   LOGICAL           :: simBenthicProd = .FALSE.
   LOGICAL           :: simCyanoRisk = .FALSE.
   LOGICAL           :: simMetalTox = .FALSE.
   INTEGER           :: simRuppiaHabitat = 1
   AED_REAL          :: mtox_lims(10)
   CHARACTER(len=40) :: mtox_vars(10)

! %% From Module Globals
!  LOGICAL :: extra_diag = .false.      !## Obsolete Use diag_level = 10
!  INTEGER :: diag_level = 10
!  %% END NAMELIST   %%  /aed_habitat_rmegacarpa/

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

   NAMELIST /aed_habitat_rmegacarpa/ &
                           simBenthicProd,   &
                           simCyanoRisk,     &
                           simRuppiaHabitat, &
                           simMetalTox, mtox_vars, mtox_lims,   &
                           rhsi_falg_link, rhsi_salg_link,      &
                           extra_diag, diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_habitat_rmegacarpa initialization"
   print *,"          WARNING! aed_habitat model is under development"

   ! Read the namelist
   read(namlst,nml=aed_habitat_rmegacarpa,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_habitat'

   ! Update module level switches
   data%num_habitats = 0
   !  data%simBenthicProd   = simBenthicProd   ; IF(simBenthicProd) data%num_habitats=data%num_habitats+1
   !  data%simMetalTox      = simMetalTox      ; IF(simMetalTox) data%num_habitats=data%num_habitats+1
   !  data%simCyanoRisk     = simCyanoRisk     ; IF(simCyanoRisk) data%num_habitats=data%num_habitats+1
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

   !-- SEAGRASS : RUPPIA MEGACARPA (adult / flower / seed germination)
   data%id_rhsi =  aed_define_sheet_diag_variable('rmegacarpa_hsi','-', 'Ruppia Habitat Suitability Index')
   data%id_rhpl =  aed_define_sheet_diag_variable('rmegacarpa_hsi_plant',  '-', 'Ruppia Habitat Suitability - plant')
   data%id_rhfl =  aed_define_sheet_diag_variable('rmegacarpa_hsi_flower', '-', 'Ruppia Habitat Suitability - flowering')
   data%id_rhsd =  aed_define_sheet_diag_variable('rmegacarpa_hsi_seed',   '-', 'Ruppia Habitat Suitability - seed germination')
  !data%id_rhtr =  aed_define_sheet_diag_variable('rmegacarpa_hsi_turion', '-', 'Ruppia Habitat Suitability - turion formation')
  !data%id_rhsp =  aed_define_sheet_diag_variable('rmegacarpa_hsi_sprout', '-', 'Ruppia Habitat Suitability - turion sprouting')
  !data%id_rhtd =  aed_define_sheet_diag_variable('rmegacarpa_hsi_dormant','-', 'Ruppia Habitat Suitability - turion viability')
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
     ! 3 stages: (1) adult, (2) flower, (3) seed germination
     ALLOCATE(data%id_d_rupfs(3))
     ALLOCATE(data%id_d_rupft(3))
     ALLOCATE(data%id_d_rupfl(3))
     ALLOCATE(data%id_d_rupfa(3))
     ALLOCATE(data%id_d_rupfd(3))
     DO i =1,3
      data%id_d_rupfs(i) = aed_define_sheet_diag_variable('rmegacarpa_hsi_fsal_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fSal')
      data%id_d_rupft(i) = aed_define_sheet_diag_variable('rmegacarpa_hsi_ftem_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fTem')
      data%id_d_rupfl(i) = aed_define_sheet_diag_variable('rmegacarpa_hsi_flgt_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fLgt')
      data%id_d_rupfa(i) = aed_define_sheet_diag_variable('rmegacarpa_hsi_falg_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fAlg')
      data%id_d_rupfd(i) = aed_define_sheet_diag_variable('rmegacarpa_hsi_fdep_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fDep')
     ENDDO
   ENDIF

   !-- GENERAL
   data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
   data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')

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
END SUBROUTINE aed_define_habitat_rmegacarpa
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian_habitat_rmegacarpa(data,column,layer_idx,pc_wet)
!-------------------------------------------------------------------------------
! Calculate benthic habitat calculations
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_rmegacarpa_data_t),INTENT(in) :: data
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

   AED_REAL :: rhpl,rhfl,rhsd
   AED_REAL :: pshpl, pshfl, pshsd, pass, height
   AED_REAL :: crns = 0.,creg = 0.,crht = 0.,crml = 0.
   AED_REAL :: limitation(3,6)

!-------------------------------------------------------------------------------
!BEGIN
    matz = 0.0 ; salt = 0.0 ; euphotic = 0.0 ; bathy = 0.0  !## CAB [-Wmaybe-uninitialized]

    depth = _STATE_VAR_(data%id_E_depth)  ! metres
    salt  = _STATE_VAR_(data%id_E_salt)   ! salinity g/L
    temp  = _STATE_VAR_(data%id_E_temp)   ! degC
    extc  = _STATE_VAR_(data%id_E_extc)   ! /m
    Io    = _STATE_VAR_S_(data%id_E_Io)   ! W/m2
    falg  =(_STATE_VAR_S_(data%id_l_falg) +   &
            _STATE_VAR_(data%id_l_salg)*depth) * 12. * 1e-3/0.5  ! convert mmolC/m2 to gDW/m2
    vel   = 0. !

    CALL rmegacarpa_habitat_suitability(data,&
                                    rhpl,rhfl,rhsd,&
                                    depth,salt,temp,extc,falg,Io,vel,pc_wet,&
                                    limitation)

    _DIAG_VAR_S_(data%id_rhpl) = rhpl
    _DIAG_VAR_S_(data%id_rhfl) = rhfl
    _DIAG_VAR_S_(data%id_rhsd) = rhsd

    ! ------------------------------------------------------------------
    ! Seed germination HSI is the instantaneous rate-style output from the
    ! current timestep. If a cumulative germination probability is wanted,
    ! it should be calculated outside this routine (for example using a
    ! start/end time window) so that the integration interval is explicit.
    ! ------------------------------------------------------------------
 
    ! Kept here for reference; this was the previous cumulative accumulator.
    ! germ_progress = _DIAG_VAR_S_(data%id_rhsd_accum)
    ! IF( pc_wet < 0.1 .OR. salt > data%sal_ge_high ) THEN
    !   germ_progress = zero_
    ! ELSE
    !   germ_rate = rhsd_inst
    !   germ_progress = MIN(one_, germ_progress + germ_rate*DDT)
    ! ENDIF
    ! _DIAG_VAR_S_(data%id_rhsd_accum) = germ_progress
    ! _DIAG_VAR_S_(data%id_rhsd)       = germ_progress
    ! ------------------------------------------------------------------

    ! Inundation time counter and wetness checker
    IF( pc_wet < 0.1 ) THEN
      _DIAG_VAR_S_(data%id_drytime) = _DIAG_VAR_S_(data%id_drytime) + DDT
      IF( _DIAG_VAR_S_(data%id_drytime)>2. ) _DIAG_VAR_S_(data%id_wettime) = zero_
    ELSE
      _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT
      IF( _DIAG_VAR_S_(data%id_wettime)>2. )_DIAG_VAR_S_(data%id_drytime) = zero_
    ENDIF

    ! Overall HSI : Habitat Suitability Index (Issue here re time integration)
    _DIAG_VAR_S_(data%id_rhsi) = (rhpl+rhfl+rhsd)/3. !for testing only. not for scientific interpretation

    IF( diag_level>1 ) THEN
      DO i = 1,3
       _DIAG_VAR_S_(data%id_d_rupfs(i)) = limitation(i,1)
       _DIAG_VAR_S_(data%id_d_rupft(i)) = limitation(i,2)
       _DIAG_VAR_S_(data%id_d_rupfl(i)) = limitation(i,3)
       _DIAG_VAR_S_(data%id_d_rupfa(i)) = limitation(i,4)
       _DIAG_VAR_S_(data%id_d_rupfd(i)) = limitation(i,5)
     ! _DIAG_VAR_S_(data%id_d_rupfm(i)) = limitation(i,6)
     ENDDO
  ENDIF

END SUBROUTINE aed_calculate_riparian_habitat_rmegacarpa
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE rmegacarpa_habitat_suitability(data,rhpl,rhfl,rhsd,depth,salt,temp,extc,fa,Io,vel,pc_wet,limitation)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_rmegacarpa_data_t),INTENT(in) :: data
   AED_REAL :: rhpl,rhfl,rhsd
   AED_REAL :: depth,salt,temp,extc,fa,Io,vel,pc_wet
   AED_REAL :: limitation(:,:)
!
!LOCALS
   AED_REAL :: rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess
   AED_REAL :: light
  !INTEGER  :: model
!
!-----------------------------------------------------------------------
!BEGIN

   rupp_salt=one_; rupp_temp=one_; rupp_lght=one_;
   rupp_falg=one_; rupp_matz=one_; rupp_dess=one_;

  !model = data%simRuppiaHabitat ! select either Gen 0 or Gen II model

   IF( depth<0.1 ) THEN
      light = 100.
   ELSE
      light = 100. * exp(-extc*(depth-0.08))
   ENDIF

   !-- First do ADULT tolerance --!
   IF ( pc_wet < 0.1 ) THEN
      ! Dry cell - set dessication factor
      rupp_dess = zero_    ! maybe need a time counter here.

   ELSE
   ! Wet cell
      rupp_salt = rmegacarpa_salinity(salt, "adult")
      rupp_temp = rmegacarpa_temp    (temp, "adult")
      rupp_lght = rmegacarpa_light   (light,"adult")
      rupp_falg = rmegacarpa_filalgae(fa,   "adult")
      rupp_dess = rmegacarpa_depth   (depth,"adult")
   ENDIF
   ! Adult plant habitat suitability
   rhpl = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)

   limitation(1,1) = rupp_salt
   limitation(1,2) = rupp_temp
   limitation(1,3) = rupp_lght
   limitation(1,4) = rupp_falg
   limitation(1,5) = rupp_dess
   limitation(1,6) = rupp_matz


   !-- Second do FLOWER tolerance --!
   IF( pc_wet < 0.1 ) THEN
     ! Dry cell - set dessication factor
     rupp_dess = zero_    ! maybe need a time counter here.

   ELSE
     ! Wet cell
     rupp_salt = rmegacarpa_salinity(salt, "flower")
     rupp_temp = rmegacarpa_temp    (temp, "flower")
     rupp_lght = rmegacarpa_light   (light,"flower")
     rupp_falg = rmegacarpa_filalgae(fa,   "flower")
     rupp_dess = rmegacarpa_depth   (depth,"flower")
   ENDIF
   ! Habitat suitability for flowering
   rhfl = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
   limitation(2,1) = rupp_salt
   limitation(2,2) = rupp_temp
   limitation(2,3) = rupp_lght
   limitation(2,4) = rupp_falg
   limitation(2,5) = rupp_dess
   limitation(2,6) = rupp_matz

   !-- Third do seed germination instantaneous rate contribution --!
   ! NOTE: unlike adult/flower, this is NOT a 0-1 suitability snapshot;
   ! it is a daily fractional germination rate (1/days-to-germinate),
   ! which can be integrated outside this routine if a cumulative probability
   ! is desired over a specific start/end time window.
   IF( pc_wet < 0.1 ) THEN
     rupp_dess = zero_
  
   ELSE
    ! Wet cell
     rupp_salt = rmegacarpa_salinity(salt, "seed")
     rupp_temp = rmegacarpa_temp    (temp, "seed")
     rupp_lght = one_
     rupp_falg = one_
     rupp_dess = rmegacarpa_depth   (depth,"seed")
   ENDIF
  ! Habitat suitability for seed germination
   !rhsd = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
     
   IF( rupp_dess <= zero_ ) THEN
     rhsd = zero_
   ELSE
     rhsd = rupp_salt
   ENDIF
   limitation(3,1) = rupp_salt
   limitation(3,2) = rupp_temp
   limitation(3,3) = rupp_lght
   limitation(3,4) = rupp_falg
   limitation(3,5) = rupp_dess
   limitation(3,6) = rupp_matz


  !---------------------------------------------------------------------
  CONTAINS

  !#############################################################################
  AED_REAL FUNCTION rmegacarpa_salinity(salt,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function 
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: salt
    CHARACTER(len=*),INTENT(in) :: stage
    AED_REAL :: optsal
  !
  !---------------------------------------------------------------------
  !BEGIN

     rmegacarpa_salinity = one_

     IF( TRIM(stage)=="adult" .OR. TRIM(stage)=="flower" ) THEN
       ! Use the supplied adult salinity curve for both adult and flower
         ! <1 unsuitable
         !  1-12 suboptimal
         !  12-40 optimal
         !  40-50 suboptimal
         ! >50 unsuitable
       IF( salt<=1. ) THEN
         rmegacarpa_salinity = zero_
       ELSE IF ( salt>1. .AND. salt<=12.  ) THEN
           rmegacarpa_salinity = 0. + ( (salt-1.)/(12.-1.) )
         ELSE IF ( salt>12. .AND. salt<=40. ) THEN
           rmegacarpa_salinity = one_
         ELSE IF ( salt>40. .AND. salt<=50. ) THEN
           rmegacarpa_salinity = 1. - ( (salt-40.)/(50.-40.) )
         ELSE IF ( salt>50. ) THEN
         rmegacarpa_salinity = zero_
       ENDIF    

     ELSEIF( TRIM(stage)=="seed" ) THEN
        !This function does not calculate HSI, but a daily fractional germination rate/progress then converted to per second
        !(i.e. 1/days-to-germinate) as a function of salinity. 
        ! which can be integrated outside this routine if a cumulative probability
        ! is desired over a specific time window.
        !<=10g/L, seed germinate in 5 days (daily germination progress 1/5 =0.2); 
        !=30g/L, seed germinate in 10 days (daily germination progress 1/10=0.1);
        !10-30g/L, germination days is linearly interpolated 
        !days-to-germination = time_low + (time_high - time_low) * (sal - sal_low) / (sal_high - sal_low);
        !>30g/L, germination progress=0
      IF (salt<=10.) THEN
        rmegacarpa_salinity = 1. /5. / 86400.
      ELSE IF (salt>10. .AND. salt <= 30.) THEN
        rmegacarpa_salinity = 1. /(5. + (10. - 5.) * (salt - 10.) / (30. - 10.)) / 86400.
      ELSE IF (salt>30.) THEN
        rmegacarpa_salinity = zero_
      ENDIF
    ENDIF

  END FUNCTION rmegacarpa_salinity
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  AED_REAL FUNCTION rmegacarpa_temp(temp,stage)
  !-----------------------------------------------------------------------------
  ! Temperature function 
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: temp
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN
     rmegacarpa_temp = one_

     IF( TRIM(stage)=="adult"  .OR. &
         TRIM(stage)=="flower" .OR. &
         TRIM(stage)=="seed"   ) THEN
           ! thresholds unknown, assume same as Ruppia tuberosa
           !  <4 unsuitable
           !   4 - 20 suboptimal
           !  20 - 23 optimal
           !  23-30 suboptimal
           ! >30 unsuitable
          IF( temp<=4. ) THEN
            rmegacarpa_temp = zero_
          ELSE IF ( temp>4. .AND. temp<=12.  ) THEN
            rmegacarpa_temp = 0. + ( (temp-4.)/(12.-4.) )
          ELSE IF ( temp>12. .AND. temp<=23. ) THEN
            rmegacarpa_temp = one_
          ELSE IF ( temp>23. .AND. temp<=30. ) THEN
            rmegacarpa_temp = 1. - ( (temp-23.)/(30.-23.) )
          ELSE IF ( temp>30. ) THEN
            rmegacarpa_temp = zero_
       ENDIF
     ENDIF

  END FUNCTION rmegacarpa_temp
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !#############################################################################
  AED_REAL FUNCTION rmegacarpa_light(light,stage)
  !-----------------------------------------------------------------------------
  ! Light function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: light
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

    rmegacarpa_light = one_

    IF( TRIM(stage)=="adult" .OR. &
        TRIM(stage)=="flower" ) THEN
       !   0 - 5 unsuitable
       !   5 - 15 suboptimal
       !   >15 optimal
       IF( light<=5.0 ) THEN
         rmegacarpa_light = zero_
       ELSE IF ( light>5.0 .AND. light<=15.  ) THEN
         rmegacarpa_light = 0. + ( (light-5.0)/(15.-5.0) )
       ELSE IF ( light>15. ) THEN
         rmegacarpa_light = one_
      ENDIF
    ENDIF

  END FUNCTION rmegacarpa_light
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  AED_REAL FUNCTION rmegacarpa_filalgae(fa,stage)
  !-----------------------------------------------------------------------------
  ! Filamentous algae function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: fa
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

    rmegacarpa_filalgae = one_

    ! The algae HSI is currently kept at 1.0 for all stages until species-
    ! specific thresholds are available.

    ! R. tuberosa threshold-based code retained for reference:
!    IF( TRIM(stage)=="adult" ) THEN
!
!      IF ( fa>=0. .AND. fa<=100. ) THEN
!        rmegacarpa_filalgae = one_
!      ELSE IF ( fa>100. .AND. fa<=368. ) THEN
!        rmegacarpa_filalgae = 1. - ( (fa-100.)/(368.-100.) )
!      ELSE IF ( fa>368. ) THEN
!        rmegacarpa_filalgae = zero_
!      ENDIF
!
!    ELSEIF( TRIM(stage)=="flower" ) THEN
!
!      IF ( fa>=0. .AND. fa<=100. ) THEN
!        rmegacarpa_filalgae = one_
!      ELSE IF ( fa>100. .AND. fa<=184. ) THEN
!        rmegacarpa_filalgae = 1. - ( (fa-100.)/(184.-100.) )
!      ELSE IF ( fa>184. ) THEN
!        rmegacarpa_filalgae = zero_
!      ENDIF
!    ENDIF

  END FUNCTION rmegacarpa_filalgae
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  AED_REAL FUNCTION rmegacarpa_depth(depth,stage)
  !-----------------------------------------------------------------------------
  ! Depth function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: depth
    CHARACTER(len=*),INTENT(in) :: stage
    AED_REAL :: maxdep
  !
  !---------------------------------------------------------------------
  !BEGIN

     rmegacarpa_depth = one_

     IF( TRIM(stage)=="adult" ) THEN
           !  <0.1 unsuitable
           !  0.1 - 0.5 suboptimal
           !  0.5 - 2 optimal
           !  2-3.5 suboptimal
           !  >3.5 unsuitable
       IF ( depth<=0.1 ) THEN
         rmegacarpa_depth = zero_
       ELSE IF ( depth>0.1 .AND. depth<=0.5 ) THEN
         rmegacarpa_depth = 1. - ( (depth-0.1)/(0.5-0.1) )
       ELSE IF ( depth>0.5 .AND. depth<=2 ) THEN
         rmegacarpa_depth = one_
       ELSE IF ( depth>2 .AND. depth<=3.5 ) THEN
         rmegacarpa_depth = 1. - ( (depth-2)/(3.5-2) )
       ELSE IF ( depth>3.5 ) THEN
         rmegacarpa_depth = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="seed" ) THEN
       ! Seeds simply need to be wet.
       IF( depth<=0.01 ) THEN
         rmegacarpa_depth = zero_
       ELSE
         rmegacarpa_depth = one_
       ENDIF

     ELSEIF( TRIM(stage)=="flower" ) THEN
           !  <0.1 unsuitable
           !  0.1 - 0.5 suboptimal
           !  0.5 - 1 optimal
           !  >1 unsuitable
       IF ( depth<=0.1 ) THEN
         rmegacarpa_depth = zero_
       ELSE IF ( depth>0.1 .AND. depth<=0.5 ) THEN
         rmegacarpa_depth = 1. - ( (depth-0.1)/(0.5-0.1) )
       ELSE IF ( depth>0.5 .AND. depth<=1 ) THEN
         rmegacarpa_depth = one_
       ELSE IF ( depth>1 ) THEN
         rmegacarpa_depth = zero_
       ENDIF
     ENDIF

  END FUNCTION rmegacarpa_depth
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE rmegacarpa_habitat_suitability
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_habitat_rmegacarpa
