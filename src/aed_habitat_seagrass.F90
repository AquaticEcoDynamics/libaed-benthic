!###############################################################################
!#                                                                             #
!# aed_habitat_seagrass.F90                                                    #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2023 - The University of Western Australia                       #
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
!# Created March 2023                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_habitat_seagrass
!-------------------------------------------------------------------------------
! aed_habitat_seagrass --- seagrass habitat model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_habitat_seagrass_data_t
!
   TYPE :: env_t
      AED_REAL :: env
      AED_REAL :: value
   END TYPE

   TYPE :: env_var_t
      CHARACTER(40) :: name
      INTEGER :: id_env
      INTEGER :: n_env
      TYPE(env_t),ALLOCATABLE :: envs(:)
   END TYPE

   TYPE :: seagrass_params_t
      CHARACTER(40) :: name
      INTEGER :: n_funcs
      TYPE(env_var_t),ALLOCATABLE :: funcs(:)
   END TYPE

!
   TYPE,extends(aed_model_data_t) :: aed_habitat_seagrass_data_t
      INTEGER :: num_habitats
      !# Variable identifiers
      INTEGER :: id_mtox
      INTEGER :: id_rhsi, id_rhpl, id_rhfl, id_rhsd, id_rhtr, id_rhsp, id_rhtd
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
      INTEGER :: id_E_nearlevel, id_E_extc, id_E_Io !, id_E_stress, id_E_airtemp

      !# Model switches
!     LOGICAL :: simBirdForaging,simBenthicProd,simFishTolerance,simGalaxiidSpawning
!     LOGICAL :: simCrabHabitat,simCharaHabitat
!     LOGICAL :: simMosquitoRisk,simCyanoRisk
!     LOGICAL :: simMetalTox,simClearWater
!     LOGICAL :: simCrocEggs,simPassiflora
      INTEGER :: simRuppiaHabitat

      INTEGER :: n_grasses
      TYPE(seagrass_params_t),ALLOCATABLE :: grasses(:)

      !# Model parameters
      AED_REAL, ALLOCATABLE :: mtox_lims(:)
      INTEGER :: num_mtox


     CONTAINS
         PROCEDURE :: define             => aed_define_habitat_seagrass
!        PROCEDURE :: calculate          => aed_calculate_habitat_seagrass
!        PROCEDURE :: calculate_seagrass => aed_calculate_benthic_habitat_seagrass
         PROCEDURE :: calculate_riparian => aed_calculate_riparian_habitat_seagrass
!        PROCEDURE :: mobility           => aed_mobility_habitat_seagrass
!        PROCEDURE :: light_extinction   => aed_light_extinction_habitat_seagrass
!        PROCEDURE :: delete             => aed_delete_habitat_seagrass

   END TYPE

!-------------------------------------------------------------------------------
!MODULE VARIABLES
   AED_REAL, PARAMETER :: DDT = 0.25/24.    ! Currently assuming 15 min timestep
   INTEGER :: diag_level = 10               ! 0 = no diagnostic outputs
                                            ! 1 = basic diagnostic outputs
                                            ! 2 = flux rates, and supporitng
                                            ! 3 = other metrics
                                            !10 = all debug & checking outputs

!===============================================================================
CONTAINS


!###############################################################################
INTEGER FUNCTION load_csv(dbase, seagrass_param, dbsize)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(seagrass_params_t),INTENT(out) :: seagrass_param(:)
   INTEGER,INTENT(out) :: dbsize
!
!LOCALS
   INTEGER :: unit, nccols
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   LOGICAL :: meh, start
   INTEGER :: n_grs = 0, g_no, i, n, env_no = 0
!
!BEGIN
!-------------------------------------------------------------------------------
   dbsize = 0
   unit = aed_rcsv_open(dbase, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   ALLOCATE(values(nccols))

   !# Do one pass to count the number of entries
   start = .true.
   DO WHILE ( aed_rcsv_read_row(unit, values) > 0 )
      IF (values(1)%length > 0) THEN
         CALL copy_name(values(1), name)
         IF ( start ) THEN
            print *,"Counting: ", name
            n_grs = n_grs + 1
         ENDIF
         start = .false.
      ELSE
         start = .true.
      ENDIF
   ENDDO
   print*,"load_csv found ",n_grs," entries"

   meh = aed_rcsv_close(unit)
   !# don't care if close fails

   !# Now a second pass to actually read the data
   unit = aed_rcsv_open(dbase, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   start = .true. ; g_no = 0
   DO WHILE ( g_no < n_grs )
      n = aed_rcsv_read_row(unit, values)
      IF (values(1)%length > 0) THEN
         ! do some stuff ....
         CALL copy_name(values(1), name)
         IF ( start ) THEN
            print *,"Decoding: ", name
            g_no = g_no + 1
            n = aed_rcsv_read_row(unit, values) ! probbaly a comment, should check
            CALL copy_name(values(1), name)
            IF ( name(1:1) == '!' ) THEN
               n = aed_rcsv_read_row(unit, values)
            ENDIF
            CALL copy_name(values(1), seagrass_param(g_no)%name)
            seagrass_param(g_no)%n_funcs = (n / 2)
            ALLOCATE(seagrass_param(g_no)%funcs(seagrass_param(g_no)%n_funcs))
            DO i = 1, seagrass_param(g_no)%n_funcs
               CALL copy_name(values((i*2)-1), seagrass_param(g_no)%funcs(i)%name)
               seagrass_param(g_no)%funcs(i)%id_env = 0
            ENDDO
            seagrass_param(g_no)%funcs(i)%n_env = extract_integer(values(i*2))
            ALLOCATE(seagrass_param(g_no)%funcs(i)%envs(seagrass_param(g_no)%funcs(i)%n_env))
            env_no = 0
         ELSE
            env_no = env_no + 1
            DO i = 1, seagrass_param(g_no)%n_funcs
               IF ( env_no <= seagrass_param(g_no)%funcs(i)%n_env ) THEN
                  seagrass_param(g_no)%funcs(i)%envs(env_no)%env   = extract_double(values((i*2)-1))
                  seagrass_param(g_no)%funcs(i)%envs(env_no)%value = extract_double(values((i*2)))
               ENDIF
            ENDDO
         ENDIF
         start = .false.
      ELSE
         start = .true.
      ENDIF
   ENDDO

   meh = aed_rcsv_close(unit)
   !# don't care if close fails

   IF (ALLOCATED(values)) DEALLOCATE(values)

   dbsize = n_grs
   load_csv = 0
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_seagrass_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_seagrass_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count   !Number of seagrass groups
   CHARACTER(40),INTENT(in)    :: list(*) !List of seagrass groups to simulate
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i, j, k, dbsize

   TYPE(seagrass_params_t),ALLOCATABLE :: seagrass_param(:)
!-------------------------------------------------------------------------------
!BEGIN
    dbsize = MAX_BVLV_TYPES !# nml cant give us a dbsize (maybe need to check name?)
    ALLOCATE(seagrass_param(MAX_BVLV_TYPES))

    IF (param_file_type(dbase) == CSV_TYPE ) THEN
       status = load_csv(dbase, seagrass_param, dbsize)
    ELSE
       print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    ENDIF

    IF (status /= 0) STOP 'Error reading seagrass database'

    ALLOCATE(data%grasses(count))
    DO i=1,count
       DO j=1,dbsize
          IF ( list(i) == seagrass_param(j)%name ) THEN
             data%grasses(i) = seagrass_param(j)
             seagrass_param(j)%name = ''
             DO K=1,seagrass_param(j)%n_funcs
                data%grasses(i)%funcs(k)%id_env = &
                              aed_locate_global(data%grasses(i)%funcs(k)%name)
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    !# Now cleanup the unused database entries
    DO j=1,dbsize
       IF ( seagrass_param(j)%name /= '' ) THEN
          DO K=1,seagrass_param(j)%n_funcs
             DEALLOCATE(seagrass_param(j)%funcs(k)%envs)
          ENDDO
          DEALLOCATE(seagrass_param(j)%funcs)
       ENDIF
    ENDDO

    DO i=1,count
       print*,"Grass Name: ", TRIM(data%grasses(i)%name)
       print*,"Funcs:"
       DO k=1,data%grasses(i)%n_funcs
          print*,"   name: ",data%grasses(i)%funcs(k)%name
       ENDDO
    ENDDO

    DEALLOCATE(seagrass_param)
END SUBROUTINE aed_seagrass_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_habitat_seagrass(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED HABITAT module
!
!  Here, the aed namelist is read and the variables exported
!  are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_habitat_seagrass_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: i, status, num_mtox

!  %% NAMELIST   %%  /aed_habitat_seagrass/
   LOGICAL            :: simBenthicProd = .FALSE.
   LOGICAL            :: simCyanoRisk = .FALSE.
   LOGICAL            :: simMetalTox = .FALSE.
   INTEGER            :: simRuppiaHabitat = 1
   AED_REAL           :: mtox_lims(10)
   CHARACTER(len=40)  :: mtox_vars(10)

   CHARACTER(len=30)  :: seagrass_list(10)
   CHARACTER(len=128) :: dbase = 'seagrass_threshold.csv'

! %% From Module Globals
!  INTEGER :: diag_level = 10
!  %% END NAMELIST   %%  /aed_habitat_seagrass/

   INTEGER :: n_grasses

   CHARACTER(len=64) :: rhsi_salg_link, rhsi_falg_link

   NAMELIST /aed_habitat_seagrass/ &
                           simBenthicProd,   &
                           simCyanoRisk,     &
                           simRuppiaHabitat, &
                           simMetalTox, mtox_vars, mtox_lims,   &
                           rhsi_falg_link, rhsi_salg_link,      &
                           dbase, n_grasses, seagrass_list,     &
                           diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_habitat_seagrass initialization"
   print *,"          WARNING! aed_habitat model is under development"

   ! Read the namelist
   read(namlst,nml=aed_habitat_seagrass,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_habitat'

   ! Update module level switches
   data%num_habitats = 0
!  data%simBenthicProd   = simBenthicProd   ; IF(simBenthicProd) data%num_habitats=data%num_habitats+1
!  data%simMetalTox      = simMetalTox      ; IF(simMetalTox) data%num_habitats=data%num_habitats+1
!  data%simCyanoRisk     = simCyanoRisk     ; IF(simCyanoRisk) data%num_habitats=data%num_habitats+1
   data%simRuppiaHabitat = simRuppiaHabitat ; IF(simRuppiaHabitat>0) data%num_habitats=data%num_habitats+1

   print *,"          ... # habitat templates simulated: ",data%num_habitats

   !----------------------------------------------------------------------------
   ! Define variables and dependencies


   !-- CONTAMINATION
   IF( simMetalTox ) THEN
     data%id_mtox =  aed_define_sheet_diag_variable('toxicity','-', 'Suitability')

!    mtox_acid_link = 'CAR_pH'
!    mtox_aass_link = 'ASS_uzaass'

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

   !-- SEAGRASS : RUPPIA
   data%id_rhsi =  aed_define_sheet_diag_variable('seagrass_hsi','-', 'Ruppia Habitat Suitability Index')
   data%id_rhpl =  aed_define_sheet_diag_variable('seagrass_hsi_plant',  '-', 'Ruppia Habitat Suitability - plant')
   data%id_rhfl =  aed_define_sheet_diag_variable('seagrass_hsi_flower', '-', 'Ruppia Habitat Suitability - flowering')
   data%id_rhsd =  aed_define_sheet_diag_variable('seagrass_hsi_seed',   '-', 'Ruppia Habitat Suitability - seed germination')
   data%id_rhtr =  aed_define_sheet_diag_variable('seagrass_hsi_turion', '-', 'Ruppia Habitat Suitability - turion formation')
   data%id_rhsp =  aed_define_sheet_diag_variable('seagrass_hsi_sprout', '-', 'Ruppia Habitat Suitability - turion sprouting')
   data%id_rhtd =  aed_define_sheet_diag_variable('seagrass_hsi_dormant','-', 'Ruppia Habitat Suitability - turion viability')

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
      data%id_d_rupfs(i) = aed_define_sheet_diag_variable('seagrass_hsi_fsal_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fSal')
      data%id_d_rupft(i) = aed_define_sheet_diag_variable('seagrass_hsi_ftem_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fTem')
      data%id_d_rupfl(i) = aed_define_sheet_diag_variable('seagrass_hsi_flgt_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fLgt')
      data%id_d_rupfa(i) = aed_define_sheet_diag_variable('seagrass_hsi_falg_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fAlg')
      data%id_d_rupfd(i) = aed_define_sheet_diag_variable('seagrass_hsi_fdep_'//CHAR(ICHAR('0') + i),'-', &
                                                                          'Ruppia Habitat Suitability - fDep')
     ENDDO
   ENDIF

   !-- GENERAL
   data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
   data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')

   ! Register variables and store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed_seagrass_load_params to values per second.
   data%n_grasses = 0
   CALL aed_seagrass_load_params(data, dbase, n_grasses, seagrass_list)

   ! Register environmental dependencies
   data%id_E_salt  = aed_locate_global('salinity')
   data%id_E_extc  = aed_locate_global('extc_coef')
   data%id_E_temp  = aed_locate_global('temperature')
   data%id_E_depth = aed_locate_global('layer_ht')
   data%id_E_bathy     = aed_locate_sheet_global('bathy')
   data%id_E_matz      = aed_locate_sheet_global('material')
   data%id_E_Io        = aed_locate_sheet_global('par_sf')
!  data%id_E_airtemp   = aed_locate_sheet_global('air_temp')
!  data%id_E_stress    = aed_locate_sheet_global('taub')
  !data%id_E_nearlevel = aed_locate_sheet_global('nearest_depth')
END SUBROUTINE aed_define_habitat_seagrass
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian_habitat_seagrass(data,column,layer_idx,pc_wet)
!-------------------------------------------------------------------------------
! Calculate benthic habitat calculations
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_seagrass_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   ! Environment
!  AED_REAL :: temp, salt, wlevel, extc, bathy, matz, Io, vel, stress, tau0, stem05, stem25
   AED_REAL :: temp, salt, extc, Io, matz, vel

   ! State
!  AED_REAL :: depth, ph, hab, sdepth, uzaass, aass, conc, mtox, turb, grav, stem, svwc
   AED_REAL :: depth

   ! Temporary variables
   INTEGER  :: i
!  AED_REAL :: bird_acid, bird_soil, bird_rveg, bird_bveg, bird_salt, bird_dept, bird_habs
!  AED_REAL :: euphotic, drytime, light
   AED_REAL :: euphotic

   ! Parameters
   AED_REAL, PARAMETER :: crit_soil_acidity = 100.000
   AED_REAL, PARAMETER :: crit_water_ph = 6.0
   AED_REAL, PARAMETER :: crit_soil_type = 5.5 ! (matz 6&7 is sand)
   AED_REAL, PARAMETER :: crit_rveg_depth = 0.6
   AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
   AED_REAL, PARAMETER :: crit_salinity = 50.0
   AED_REAL, PARAMETER :: crit_leg_depth = 0.12
   AED_REAL, PARAMETER :: crit_hab_conc = 500.

!  AED_REAL :: fs_sdepth , fs_substr, fs_spntem, fs_stress, fs_dewatr, fs_mattem

   AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp =0.,falg,rhtd=1.0
!  AED_REAL :: pshpl, pshfl, pshsd, pass, height
!  AED_REAL :: crns = 0.,creg = 0.,crht = 0.,crml = 0.
   AED_REAL :: limitation(6,6)

!-------------------------------------------------------------------------------
!BEGIN
!   matz = 0.0 ; salt = 0.0 ; euphotic = 0.0 ; bathy = 0.0  !## CAB [-Wmaybe-uninitialized]
    matz = 0.0 ; euphotic = 0.0  !## CAB [-Wmaybe-uninitialized]

    depth = _STATE_VAR_(data%id_E_depth)  ! metres
    salt  = _STATE_VAR_(data%id_E_salt)   ! salinity g/L
    temp  = _STATE_VAR_(data%id_E_temp)   ! degC
    extc  = _STATE_VAR_(data%id_E_extc)   ! /m
    Io    = _STATE_VAR_S_(data%id_E_Io)   ! W/m2
    falg  =(_STATE_VAR_S_(data%id_l_falg) +   &
            _STATE_VAR_(data%id_l_salg)*depth) * 12. * 1e-3/0.5  ! convert mmolC/m2 to gDW/m2
    vel   = 0. !

    CALL seagrass_habitat_suitability(data,&
                                    rhpl,rhfl,rhsd,rhtr,rhsp,rhtd,&
                                  ! depth,salt,temp,extc,falg,Io,vel,pc_wet,&
                                    depth,salt,temp,extc,falg,pc_wet,&
                                    limitation)

    _DIAG_VAR_S_(data%id_rhpl) = rhpl
    _DIAG_VAR_S_(data%id_rhfl) = rhfl
    _DIAG_VAR_S_(data%id_rhsd) = rhsd
    _DIAG_VAR_S_(data%id_rhtr) = rhtr
    _DIAG_VAR_S_(data%id_rhsp) = rhsp
    _DIAG_VAR_S_(data%id_rhsp) = rhtd

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
     ! _DIAG_VAR_S_(data%id_d_rupfm(i)) = limitation(i,6)
     ENDDO
  ENDIF

END SUBROUTINE aed_calculate_riparian_habitat_seagrass
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!SUBROUTINE seagrass_habitat_suitability(data,rhpl,rhfl,rhsd,rhtr,rhsp,rhtd,depth,salt,temp,extc,fa,Io,vel,pc_wet,limitation)
SUBROUTINE seagrass_habitat_suitability(data,rhpl,rhfl,rhsd,rhtr,rhsp,rhtd,depth,salt,temp,extc,fa,pc_wet,limitation)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_habitat_seagrass_data_t),INTENT(in) :: data
   AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp,rhtd
!  AED_REAL :: depth,salt,temp,extc,fa,Io,vel,pc_wet
   AED_REAL :: depth,salt,temp,extc,fa,pc_wet
   AED_REAL :: limitation(:,:)
!
!LOCALS
   AED_REAL :: rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess
   AED_REAL :: light
   INTEGER  :: model
!
!-----------------------------------------------------------------------
!BEGIN

   rupp_salt=one_; rupp_temp=one_; rupp_lght=one_;
   rupp_falg=one_; rupp_matz=one_; rupp_dess=one_;

   model = data%simRuppiaHabitat ! select either Gen 0 or Gen II model

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
      rupp_salt = seagrass_salinity(salt, "adult")
      rupp_temp = seagrass_temp    (temp, "adult")
      rupp_lght = seagrass_light   (light,"adult")
      rupp_falg = seagrass_filalgae(fa,   "adult")
      rupp_dess = seagrass_depth   (depth,"adult")
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
     rupp_salt = seagrass_salinity(salt, "flower")
     rupp_temp = seagrass_temp    (temp, "flower")
     rupp_lght = seagrass_light   (light,"flower")
     rupp_falg = seagrass_filalgae(fa,   "flower")
     rupp_dess = seagrass_depth   (depth,"flower")
   ENDIF
   ! Habitat suitability for flowering
   rhfl = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
   limitation(2,1) = rupp_salt
   limitation(2,2) = rupp_temp
   limitation(2,3) = rupp_lght
   limitation(2,4) = rupp_falg
   limitation(2,5) = rupp_dess
   limitation(2,6) = rupp_matz

   !-- Third do seed germination
   IF( pc_wet < 0.1 ) THEN
     ! Dry cell - set dessication factor

     rupp_dess = 0.5    ! maybe need a time counter here.

   ELSE
     ! Wet cell
     rupp_salt = seagrass_salinity(salt, "seed")
     rupp_temp = seagrass_temp    (temp, "seed")
     rupp_lght = one_
     rupp_falg = one_
     rupp_dess = seagrass_depth   (depth,"seed")

   ENDIF
   ! Habitat suitability for seed germination
   rhsd = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
   limitation(3,1) = rupp_salt
   limitation(3,2) = rupp_temp
   limitation(3,3) = rupp_lght
   limitation(3,4) = rupp_falg
   limitation(3,5) = rupp_dess
   limitation(3,6) = rupp_matz


   !-- Fourth do sediment suitability for turion sprouting
   IF( pc_wet < 0.1 ) THEN
     ! Dry cell - set dessication factor

     rupp_dess = zero_    ! maybe need a time counter here.

   ELSE
     ! Wet cell
     rupp_salt = seagrass_salinity(salt, "turion")
     rupp_temp = seagrass_temp    (temp, "turion")
     rupp_lght = seagrass_light   (light,"turion")
     rupp_falg = one_
     rupp_dess = seagrass_depth   (depth,"turion")

   ENDIF
   !IF( .NOT. _STATE_VAR_S_(data%id_E_matz)==in_zone_set ) rupp_matz = zero_
   ! overall habitat suitability for seed germination
   rhtr = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
   limitation(4,1) = rupp_salt
   limitation(4,2) = rupp_temp
   limitation(4,3) = rupp_lght
   limitation(4,4) = rupp_falg
   limitation(4,5) = rupp_dess
   limitation(4,6) = rupp_matz

   !-- Fifth do sprout tolerance
   IF( pc_wet < 0.1 ) THEN
      ! Dry cell - set dessication factor

      rupp_dess = zero_    ! maybe need a time counter here.

   ELSE
      ! Wet cell
      rupp_salt = seagrass_salinity(salt, "sprout")
      rupp_temp = seagrass_temp    (temp, "sprout")
      rupp_lght = seagrass_light   (light,"sprout")
      rupp_falg = seagrass_filalgae(fa,   "sprout")
      rupp_dess = seagrass_depth   (depth,"sprout")

   ENDIF
   ! overall habitat suitability for sprouting
   rhsp = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
   limitation(5,1) = rupp_salt
   limitation(5,2) = rupp_temp
   limitation(5,3) = rupp_lght
   limitation(5,4) = rupp_falg
   limitation(5,5) = rupp_dess
   limitation(5,6) = rupp_matz

   !-- Sixth do turion viability during dormancy
   IF( pc_wet < 0.1 ) THEN
      ! Dry cell - set dessication factor

      rupp_dess = one_    ! maybe need a time counter here.
   ELSE
      ! Wet cell
      rupp_salt = seagrass_salinity(salt, "dormant")
      rupp_temp = one_
      rupp_lght = one_
      rupp_falg = one_
      rupp_dess = one_
   ENDIF

   ! overall habitat suitability for turion viability
   rhtd = MIN(rupp_salt,rupp_temp,rupp_lght,rupp_falg,rupp_matz,rupp_dess)
   limitation(6,1) = rupp_salt
   limitation(6,2) = rupp_temp
   limitation(6,3) = rupp_lght
   limitation(6,4) = rupp_falg
   limitation(6,5) = rupp_dess
   limitation(6,6) = rupp_matz


  !---------------------------------------------------------------------
  CONTAINS

  !#############################################################################
  AED_REAL FUNCTION seagrass_salinity(salt,stage)
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

     seagrass_salinity = one_

     IF ( TRIM(stage)=="seed" ) THEN

       IF ( model==1 ) THEN
         ! COORONG GENERATION 0

         ! <10 unsuitable
         !  10-30 suboptimal
         !  40-60 optimal
         !  60-85 suboptimal
         ! >85 unsuitable
         IF( salt<=10. ) THEN
           seagrass_salinity = zero_
         ELSE IF ( salt>10. .AND. salt<=30.  ) THEN
           seagrass_salinity = 0. + ( (salt-10.)/(30.-10.) )
         ELSE IF ( salt>30. .AND. salt<=60. ) THEN
           seagrass_salinity = one_
         ELSE IF ( salt>60. .AND. salt<=85. ) THEN
           seagrass_salinity = 1. - ( (salt-60.)/(85.-60.) )
         ELSE IF ( salt>85. ) THEN
           seagrass_salinity = zero_
         ENDIF
       ELSEIF( model==2 ) THEN
         ! COORONG GENERATION II

         ! <1 unsuitable
         !  1-5 suboptimal
         !  5-40 optimal
         !  40-85 suboptimal
         ! >85 unsuitable
         IF( salt<=1. ) THEN
           seagrass_salinity = zero_
         ELSE IF ( salt>1. .AND. salt<=5.  ) THEN
           seagrass_salinity = 0. + ( (salt-1.)/(5.-1.) )
         ELSE IF ( salt>5. .AND. salt<=40. ) THEN
           seagrass_salinity = one_
         ELSE IF ( salt>40. .AND. salt<=85. ) THEN
           seagrass_salinity = 1. - ( (salt-40.)/(85.-40.) )
         ELSE IF ( salt>85. ) THEN
           seagrass_salinity = zero_
         ENDIF
       ENDIF

     ELSEIF( TRIM(stage)=="sprout" ) THEN

       optsal = 130. ; IF(model==2) optsal = 125.
       ! <20 suboptimal
       !  20-75 optimal
       !  75 - 130 suboptimal
       ! >130 unsuitable
       IF( salt<=0.1 ) THEN
         seagrass_salinity = zero_
       ELSE IF ( salt>0.1 .AND. salt<=20.  ) THEN
         seagrass_salinity = 0. + ( (salt-0.1)/(20.-0.1) )
       ELSE IF ( salt>20. .AND. salt<=75. ) THEN
         seagrass_salinity = one_
       ELSE IF ( salt>75. .AND. salt<=optsal ) THEN
         seagrass_salinity = 1. - ( (salt-75.)/(optsal-75.) )
       ELSE IF ( salt>optsal ) THEN
         seagrass_salinity = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="adult" ) THEN

      IF( model==1 ) THEN
        ! COORONG GENERATION 0

        !    0 - 10 unsuitable
        !   10 - 71 suboptimal
        !   72 - 123 optimal
        !  123 - 230 suboptimal
        ! >230 unsuitable
        IF( salt<=10. ) THEN
         seagrass_salinity = zero_
        ELSE IF ( salt>10. .AND. salt<=31.  ) THEN
         seagrass_salinity = 0. + ( (salt-10.)/(31.-10.) )
        ELSE IF ( salt>31. .AND. salt<=123. ) THEN
         seagrass_salinity = one_
        ELSE IF ( salt>123. .AND. salt<=230. ) THEN
         seagrass_salinity = 1. - ( (salt-123.)/(230.-123.) )
        ELSE IF ( salt>230. ) THEN
         seagrass_salinity = zero_
        ENDIF
      ELSEIF( model==2 ) THEN
        ! COORONG GENERATION II

        !    0 - 10 unsuitable
        !   10 - 19 suboptimal
        !   19 - 124 optimal
        !  124 - 230 suboptimal
        ! >230 unsuitable
        IF( salt<=10. ) THEN
          seagrass_salinity = zero_
        ELSE IF ( salt>10. .AND. salt<=19.  ) THEN
          seagrass_salinity = 0. + ( (salt-10.)/(19.-10.) )
        ELSE IF ( salt>19. .AND. salt<=124. ) THEN
          seagrass_salinity = one_
        ELSE IF ( salt>124. .AND. salt<=230. ) THEN
          seagrass_salinity = 1. - ( (salt-124.)/(230.-124.) )
        ELSE IF ( salt>230. ) THEN
          seagrass_salinity = zero_
        ENDIF
      ENDIF

     ELSEIF( TRIM(stage)=="flower" ) THEN

       IF( model==1 ) THEN
         ! COORONG GENERATION 0

         !  <10 unsuitable
         !   10 - 35 suboptimal
         !   35 - 62 optimal
         !   62 - 100 suboptimal
         ! >100 unsuitable
         IF( salt<=10. ) THEN
          seagrass_salinity = zero_
         ELSE IF ( salt>10. .AND. salt<=35.  ) THEN
          seagrass_salinity = 0. + ( (salt-10.)/(35.-10.) )
         ELSE IF ( salt>35. .AND. salt<=62. ) THEN
          seagrass_salinity = one_
         ELSE IF ( salt>62. .AND. salt<=100. ) THEN
          seagrass_salinity = 1. - ( (salt-62.)/(100.-62.) )
         ELSE IF ( salt>100. ) THEN
          seagrass_salinity = zero_
         ENDIF

       ELSEIF( model==2 ) THEN
         ! COORONG GENERATION II

         !  <12 unsuitable
         !   12 - 47 suboptimal
         !   47 - 62 optimal
         !   62 - 100 suboptimal
         ! >100 unsuitable
         IF( salt<=10. ) THEN
          seagrass_salinity = zero_
         ELSE IF ( salt>10. .AND. salt<=47.  ) THEN
          seagrass_salinity = 0. + ( (salt-10.)/(47.-10.) )
         ELSE IF ( salt>47. .AND. salt<=62. ) THEN
          seagrass_salinity = one_
         ELSE IF ( salt>62. .AND. salt<=100. ) THEN
          seagrass_salinity = 1. - ( (salt-62.)/(100.-62.) )
         ELSE IF ( salt>100. ) THEN
          seagrass_salinity = zero_
         ENDIF
       ENDIF

     ELSEIF( TRIM(stage)=="turion" ) THEN

       IF( model==1 ) THEN
         ! COORONG GENERATION 0

         ! <70 unsuitable
         !  70 - 124 suboptimal
         !  124 - 160 optimal
         !  160 - 230 suboptimal
         ! >230 unsuitable
         IF( salt<=70. ) THEN
           seagrass_salinity = zero_
         ELSE IF ( salt>70. .AND. salt<=124.  ) THEN
           seagrass_salinity = 0. + ( (salt-70.)/(124.-70.) )
         ELSE IF ( salt>124. .AND. salt<=160. ) THEN
           seagrass_salinity = one_
         ELSE IF ( salt>160. .AND. salt<=230. ) THEN
           seagrass_salinity = 1. - ( (salt-160.)/(230.-160.) )
         ELSE IF ( salt>230. ) THEN
           seagrass_salinity = zero_
         ENDIF

       ELSEIF( model==2 ) THEN
         ! COORONG GENERATION II

         ! <40 unsuitable
         !  40 - 70 suboptimal
         !  70 - 160 optimal
         !  160 - 230 suboptimal
         ! >230 unsuitable
         IF( salt<=40. ) THEN
           seagrass_salinity = zero_
         ELSE IF ( salt>40. .AND. salt<=70.  ) THEN
           seagrass_salinity = 0. + ( (salt-40.)/(70.-40.) )
         ELSE IF ( salt>70. .AND. salt<=160. ) THEN
           seagrass_salinity = one_
         ELSE IF ( salt>160. .AND. salt<=230. ) THEN
           seagrass_salinity = 1. - ( (salt-160.)/(230.-160.) )
         ELSE IF ( salt>230. ) THEN
           seagrass_salinity = zero_
         ENDIF

       ENDIF

     ELSEIF( TRIM(stage)=="dormant" ) THEN

       IF (model==2 ) THEN
         IF ( salt<=135. ) THEN
           seagrass_salinity = one_
         ELSE IF ( salt>135. .AND. salt<=165. ) THEN
           seagrass_salinity = 1. - ( (salt-135.)/(165.-135.) )
         ELSE IF ( salt>165. ) THEN
           seagrass_salinity = zero_
         ENDIF
       ENDIF

     ENDIF

  END FUNCTION seagrass_salinity
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  AED_REAL FUNCTION seagrass_temp(temp,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: temp
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN
     seagrass_temp = one_

     IF( model==1 ) THEN
       ! COORONG GENERATION 0

       IF( TRIM(stage)=="seed"   .OR. &
           TRIM(stage)=="sprout" .OR. &
           TRIM(stage)=="adult"  .OR. &
           TRIM(stage)=="flower" .OR. &
           TRIM(stage)=="turion"      ) THEN

         !  <4 unsuitable
         !   4 - 10 suboptimal
         !  10 - 20 optimal
         !  20-30 suboptimal
         ! >30 unsuitable
         IF( temp<=4. ) THEN
          seagrass_temp = zero_
         ELSE IF ( temp>4. .AND. temp<=10.  ) THEN
          seagrass_temp = 0. + ( (temp-4.)/(10.-4.) )
         ELSE IF ( temp>10. .AND. temp<=20. ) THEN
          seagrass_temp = one_
         ELSE IF ( temp>20. .AND. temp<=30. ) THEN
          seagrass_temp = 1. - ( (temp-20.)/(30.-20.) )
         ELSE IF ( temp>30. ) THEN
          seagrass_temp = zero_
         ENDIF
       ENDIF

     ELSEIF( model==2 ) THEN
         ! COORONG GENERATION II

         IF( TRIM(stage)=="adult"  .OR. &
             TRIM(stage)=="seed"   .OR. &
             TRIM(stage)=="sprout" .OR. &
             TRIM(stage)=="flower" .OR. &
             TRIM(stage)=="turion"      ) THEN

           !  <4 unsuitable
           !   4 - 20 suboptimal
           !  20 - 23 optimal
           !  23-30 suboptimal
           ! >30 unsuitable
           IF( temp<=4. ) THEN
            seagrass_temp = zero_
          ELSE IF ( temp>4. .AND. temp<=12.  ) THEN
            seagrass_temp = 0. + ( (temp-4.)/(12.-4.) )
          ELSE IF ( temp>12. .AND. temp<=23. ) THEN
            seagrass_temp = one_
          ELSE IF ( temp>23. .AND. temp<=30. ) THEN
            seagrass_temp = 1. - ( (temp-23.)/(30.-23.) )
           ELSE IF ( temp>30. ) THEN
            seagrass_temp = zero_
           ENDIF
         ENDIF
     ENDIF

  END FUNCTION seagrass_temp
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !#############################################################################
  AED_REAL FUNCTION seagrass_light(light,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: light
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

    seagrass_light = one_

    IF( model==1 ) THEN
     ! COORONG GENERATION 0

     IF( TRIM(stage)=="sprout" .OR. &
         TRIM(stage)=="adult"  ) THEN
       !  0 - 7.5 unsuitable
       !  7.5 - 24 suboptimal
       ! >24 optimal
       IF( light<=7.5 ) THEN
         seagrass_light = zero_
       ELSE IF ( light>7.5 .AND. light<=24.  ) THEN
         seagrass_light = 0. + ( (light-7.5)/(24.-7.5) )
       ELSE IF ( light>24. ) THEN
         seagrass_light = one_
       ENDIF
     ENDIF

    ELSEIF( model==2 ) THEN
     ! COORONG GENERATION II

     IF( TRIM(stage)=="sprout" .OR. &
         TRIM(stage)=="turion" .OR. &
         TRIM(stage)=="flower" .OR. &
         TRIM(stage)=="adult"  ) THEN
       !   0 - 5 unsuitable
       !   5 - 36 suboptimal
       ! >36 optimal
       IF( light<=5.0 ) THEN
         seagrass_light = zero_
       ELSE IF ( light>5.0 .AND. light<=36.  ) THEN
         seagrass_light = 0. + ( (light-5.0)/(36.-5.0) )
       ELSE IF ( light>36. ) THEN
         seagrass_light = one_
       ENDIF
     ENDIF
    ENDIF

  END FUNCTION seagrass_light
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  AED_REAL FUNCTION seagrass_filalgae(fa,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: fa
    CHARACTER(len=*),INTENT(in) :: stage
  !
  !---------------------------------------------------------------------
  !BEGIN

    seagrass_filalgae = one_

    IF( model==1 ) THEN
     ! COORONG GENERATION 0

     IF( TRIM(stage)=="adult"  .OR. &
         TRIM(stage)=="flower" ) THEN

       IF ( fa>=0. .AND. fa<=25. ) THEN
         seagrass_filalgae = one_
       ELSE IF ( fa>25. .AND. fa<=100. ) THEN
         seagrass_filalgae = 1. - ( (fa-25.)/(100.-25.) )
       ELSE IF ( fa>100. ) THEN
         seagrass_filalgae = zero_
       ENDIF
     ENDIF

    ELSEIF( model==2 ) THEN
     ! COORONG GENERATION II

     IF( TRIM(stage)=="adult" ) THEN

       IF ( fa>=0. .AND. fa<=100. ) THEN
         seagrass_filalgae = one_
       ELSE IF ( fa>100. .AND. fa<=368. ) THEN
         seagrass_filalgae = 1. - ( (fa-100.)/(368.-100.) )
       ELSE IF ( fa>368. ) THEN
         seagrass_filalgae = zero_
       ENDIF

     ELSEIF( TRIM(stage)=="flower" ) THEN

       IF ( fa>=0. .AND. fa<=100. ) THEN
         seagrass_filalgae = one_
       ELSE IF ( fa>100. .AND. fa<=184. ) THEN
         seagrass_filalgae = 1. - ( (fa-100.)/(184.-100.) )
       ELSE IF ( fa>184. ) THEN
         seagrass_filalgae = zero_
       ENDIF
     ENDIF

    ENDIF

  END FUNCTION seagrass_filalgae
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !#############################################################################
  AED_REAL FUNCTION seagrass_depth(depth,stage)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    AED_REAL,INTENT(in) :: depth
    CHARACTER(len=*),INTENT(in) :: stage
    AED_REAL :: maxdep
  !
  !---------------------------------------------------------------------
  !BEGIN

     seagrass_depth = one_

     IF( TRIM(stage)=="adult" ) THEN

       IF ( depth<=0.1 ) THEN
         seagrass_depth = zero_
       ELSE IF ( depth>0.1 .AND. depth<=0.2 ) THEN
         seagrass_depth = 1. - ( (depth-0.1)/(0.2-0.1) )
       ELSE IF ( depth>0.2 ) THEN
         seagrass_depth = one_
       ENDIF

     ELSEIF( TRIM(stage)=="sprout" ) THEN

       IF ( depth<=0.01 ) THEN
         seagrass_depth = zero_
       ELSE IF ( depth>0.01 .AND. depth<=0.2 ) THEN
         seagrass_depth = 1. - ( (depth-0.01)/(0.2-0.01) )
       ELSE IF ( depth>0.2 ) THEN
         seagrass_depth = one_
       ENDIF

     ELSEIF( TRIM(stage)=="flower" ) THEN

       maxdep = 1.0 ; IF(model==2) maxdep=0.9

       IF ( depth<=0.01 ) THEN
         seagrass_depth = zero_
       ELSE IF ( depth>0.01 .AND. depth<=0.1 ) THEN
         seagrass_depth = 1. - ( (depth-0.01)/(0.1-0.01) )
       ELSE IF ( depth>0.1 .AND. depth<=0.4 ) THEN
         seagrass_depth = one_
       ELSE IF ( depth>0.4 .AND. depth<=maxdep ) THEN
         seagrass_depth = 1. - ( (depth-0.4)/(maxdep-0.4) )
       ELSE IF ( depth>maxdep ) THEN
         seagrass_depth = zero_
       ENDIF

     ENDIF

  END FUNCTION seagrass_depth
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE seagrass_habitat_suitability
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_habitat_seagrass
