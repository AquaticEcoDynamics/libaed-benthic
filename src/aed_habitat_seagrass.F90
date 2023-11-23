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
      INTEGER :: id_d_env
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
!  INTEGER,ALLOCATABLE :: id_d_rupfs(:), id_d_rupft(:), id_d_rupfl(:), id_d_rupfa(:), id_d_rupfd(:)
      INTEGER :: id_rhsi, id_rhfl, id_rhsd, id_rhtr, id_rhsp, id_rhtd
      INTEGER :: id_wettime, id_drytime
!     INTEGER, ALLOCATABLE :: id_rhpl(:)

      !# Dependencies
      INTEGER :: id_l_salg, id_l_falg, id_d_turb, id_l_ncs1, id_l_ncs2, id_l_tau0

      !# Environment variables
      INTEGER :: id_E_temp, id_E_salt, id_E_bathy, id_E_matz, id_E_depth
      INTEGER :: id_E_nearlevel, id_E_extc, id_E_Io

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

! print*,'############################################################'
   !# Do one pass to count the number of entries
   start = .true.
   DO WHILE ( aed_rcsv_read_row(unit, values) > 0 )
      IF (values(1)%length > 0) THEN
         CALL copy_name(values(1), name)
         IF ( start ) THEN
            print *,"Counting: '", TRIM(name), "'"
            n_grs = n_grs + 1
! else
! print *,"Ignoring: '", TRIM(name), "' length = ", values(1)%length
         ENDIF
         start = .false.
      ELSE
! print*,"zero length"
         start = .true.
      ENDIF
   ENDDO
   print*,"load_csv found ",n_grs," entries"
! print*,'############################################################'

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
            IF ( name(1:1) == '!' ) THEN        !# yes, it was a comment
               n = aed_rcsv_read_row(unit, values)
            ENDIF
            CALL copy_name(values(1), seagrass_param(g_no)%name)
            seagrass_param(g_no)%n_funcs = (n / 2)
            ALLOCATE(seagrass_param(g_no)%funcs(seagrass_param(g_no)%n_funcs))
            DO i = 1, seagrass_param(g_no)%n_funcs
               CALL copy_name(values((i*2)-1), seagrass_param(g_no)%funcs(i)%name)
               seagrass_param(g_no)%funcs(i)%id_env = 0
               seagrass_param(g_no)%funcs(i)%n_env = extract_integer(values(i*2))
               ALLOCATE(seagrass_param(g_no)%funcs(i)%envs(seagrass_param(g_no)%funcs(i)%n_env))
            ENDDO
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
   INTEGER,INTENT(in)          :: count   ! Number of seagrass groups
   CHARACTER(30),INTENT(in)    :: list(*) ! List of seagrass groups to simulate
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
       data%grasses(i)%name = ''
       data%grasses(i)%n_funcs = 0
!print*,"looking for ", list(i)
       DO j=1,dbsize
          IF ( list(i) == seagrass_param(j)%name ) THEN
             data%grasses(i) = seagrass_param(j)
             seagrass_param(j)%name = ''
             DO k=1,seagrass_param(j)%n_funcs
                data%grasses(i)%funcs(k)%id_env = &
                              aed_locate_global(data%grasses(i)%funcs(k)%name)
             ENDDO
!ELSE
!print*,"Its not '",TRIM(seagrass_param(j)%name),"'"
          ENDIF
       ENDDO
    ENDDO

    !# Now cleanup the unused database entries
    DO j=1,dbsize
       IF ( seagrass_param(j)%name /= '' ) THEN
          DO K=1,seagrass_param(j)%n_funcs
             IF ( ALLOCATED(seagrass_param(j)%funcs(k)%envs) ) &
                 DEALLOCATE(seagrass_param(j)%funcs(k)%envs)
          ENDDO
          IF ( ALLOCATED(seagrass_param(j)%funcs) ) &
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
   INTEGER :: g, f, status

!  %% NAMELIST   %%  /aed_habitat_seagrass/

   CHARACTER(len=128) :: dbase = 'seagrass_threshold.csv'

   INTEGER :: n_grasses
   CHARACTER(len=30)  :: seagrass_list(10)

! %% From Module Globals
!  INTEGER :: diag_level = 10
!  %% END NAMELIST   %%  /aed_habitat_seagrass/

   NAMELIST /aed_habitat_seagrass/ &
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
!   data%simBenthicProd   = simBenthicProd   ; IF(simBenthicProd) data%num_habitats=data%num_habitats+1
!   data%simMetalTox      = simMetalTox      ; IF(simMetalTox) data%num_habitats=data%num_habitats+1
!   data%simCyanoRisk     = simCyanoRisk     ; IF(simCyanoRisk) data%num_habitats=data%num_habitats+1
!   data%simRuppiaHabitat = simRuppiaHabitat ; IF(simRuppiaHabitat>0) data%num_habitats=data%num_habitats+1

   print *,"          ... # habitat templates simulated: ",data%num_habitats

   !----------------------------------------------------------------------------
   ! Define variables and dependencies

   !-- SEAGRASS : diagnostic outputs of each species/group
!  data%id_rhsi =  aed_define_sheet_diag_variable('seagrass_hsi','-', 'Ruppia Habitat Suitability Index')

!  IF (rhsi_falg_link .EQ. "") THEN
!      STOP 'need to set rhsi_falg_link and rhsi_salg_link'
!  ENDIF
!
!  data%id_l_salg  = aed_locate_variable(TRIM(rhsi_salg_link))
!  data%id_l_falg  = aed_locate_sheet_variable(TRIM(rhsi_falg_link))


   ! Register variables and store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed_seagrass_load_params to values per second.
   data%n_grasses = 0
   CALL aed_seagrass_load_params(data, dbase, n_grasses, seagrass_list)
   data%n_grasses = n_grasses

   IF (diag_level>1) THEN
!     ALLOCATE(data%id_rhpl(n_grasses))
      DO g=1,n_grasses
         DO f=1,data%grasses(g)%n_funcs
            data%grasses(g)%funcs(f)%id_d_env = aed_define_sheet_diag_variable(&
                         data%grasses(g)%name//data%grasses(g)%funcs(f)%name,  &
                                 '-', 'Seagrass Habitat Suitability')
         ENDDO
      ENDDO
   ENDIF

   !-- GENERAL
!  data%id_wettime = aed_define_sheet_diag_variable('wettime','d','time cell has been innundated')
!  data%id_drytime = aed_define_sheet_diag_variable('drytime','d','time cell has been exposed')

   ! Register environmental dependencies
   data%id_E_salt  = aed_locate_global('salinity')
   data%id_E_extc  = aed_locate_global('extc_coef')
   data%id_E_temp  = aed_locate_global('temperature')
   data%id_E_depth = aed_locate_global('layer_ht')
   data%id_E_bathy = aed_locate_sheet_global('bathy')
   data%id_E_matz  = aed_locate_sheet_global('material')
   data%id_E_Io    = aed_locate_sheet_global('par_sf')
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
   AED_REAL :: temp, salt, extc, Io, vel

   ! State
   AED_REAL :: depth

   ! Temporary variables
   INTEGER  :: i, j

   ! Parameters
!  AED_REAL, PARAMETER :: crit_soil_acidity = 100.000
!  AED_REAL, PARAMETER :: crit_water_ph = 6.0
!  AED_REAL, PARAMETER :: crit_soil_type = 5.5 ! (matz 6&7 is sand)
!  AED_REAL, PARAMETER :: crit_rveg_depth = 0.6
!  AED_REAL, PARAMETER :: sav_height = 0.1 !assume plant are 10cm to middle
!  AED_REAL, PARAMETER :: crit_salinity = 50.0
!  AED_REAL, PARAMETER :: crit_leg_depth = 0.12
!  AED_REAL, PARAMETER :: crit_hab_conc = 500.

!  AED_REAL :: rhpl,rhfl,rhsd,rhtr,rhsp =0.,falg,rhtd=1.0
!  AED_REAL :: limitation(6,6)

!-------------------------------------------------------------------------------
!BEGIN
    depth = _STATE_VAR_(data%id_E_depth)  ! metres
    salt  = _STATE_VAR_(data%id_E_salt)   ! salinity g/L
    temp  = _STATE_VAR_(data%id_E_temp)   ! degC
    extc  = _STATE_VAR_(data%id_E_extc)   ! /m
    Io    = _STATE_VAR_S_(data%id_E_Io)   ! W/m2
!   falg  =(_STATE_VAR_S_(data%id_l_falg) +   &
!           _STATE_VAR_(data%id_l_salg)*depth) * 12. * 1e-3/0.5  ! convert mmolC/m2 to gDW/m2
    vel   = 0. !

    ! This is not right. CAB
    DO i=1,data%n_grasses
       DO j=1,data%grasses(i)%n_funcs
    !     _DIAG_VAR_S_(data%id_rhpl(i)) = &
          _DIAG_VAR_S_(data%grasses(i)%funcs(j)%id_d_env) = &
                                    seagrass_limitation(data%grasses(i)%funcs(j))
       ENDDO
    ENDDO

    ! Inundation time counter and wetness checker
    IF( pc_wet < 0.1 ) THEN
      _DIAG_VAR_S_(data%id_drytime) = _DIAG_VAR_S_(data%id_drytime) + DDT
      IF( _DIAG_VAR_S_(data%id_drytime)>2. ) _DIAG_VAR_S_(data%id_wettime) = zero_
    ELSE
      _DIAG_VAR_S_(data%id_wettime) = _DIAG_VAR_S_(data%id_wettime) + DDT
      IF( _DIAG_VAR_S_(data%id_wettime)>2. )_DIAG_VAR_S_(data%id_drytime) = zero_
    ENDIF

!-------------------------------------------------------------------------------
  CONTAINS

  !#############################################################################
  AED_REAL FUNCTION seagrass_limitation(func)
  !-----------------------------------------------------------------------------
  ! Salinity function
  !-----------------------------------------------------------------------------
  !ARGUMENTS
    TYPE(env_var_t),INTENT(in) :: func
  !
    AED_REAL :: val
    INTEGER :: i
  !
  !---------------------------------------------------------------------
  !BEGIN

     val = _STATE_VAR_(func%id_env)
     seagrass_limitation = zero_

     i = 1
     DO WHILE ( i < func%n_env )
        IF ( val > func%envs(i)%env .AND. val <= func%envs(i+1)%env ) THEN
           IF ( func%envs(i)%value == func%envs(i+1)%value ) THEN
              seagrass_limitation = val
           ELSE
              seagrass_limitation = (val - func%envs(i)%value) / (func%envs(i+1)%value - func%envs(i)%value)

              IF ( func%envs(i)%value > func%envs(i+1)%value ) THEN
                 seagrass_limitation = val - seagrass_limitation
              ENDIF
           ENDIF
           EXIT  ! Drop out - we're done
        ENDIF
        i = i + 1
     ENDDO

  END FUNCTION seagrass_limitation
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE aed_calculate_riparian_habitat_seagrass
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_habitat_seagrass
