!###############################################################################
!#                                                                             #
!# aed_bivalve.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  In collaboration with :                                                    #
!#     Cornell University, Department of Natural Resources                     #
!#                                                                             #
!#  Copyright 2015 - 2022 -  The University of Western Australia               #
!#                                                                             #
!#   AED2 is free software: you can redistribute it and/or modify              #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED2 is distributed in the hope that it will be useful,                   #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#   For use in a commercial seeting please contact the authors.               #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Originally created January 2015 by Matthew Hipsey & Lars Rudstam            #
!# Follow updates @ https://github.com/AquaticEcoDynamics/libaed2-plus         #
!#                                                                             #
!###############################################################################

#include "aed.h"

#define _MODPFX_ 3
#define _PHYMOD_ 'PHY'
#define _OGMPOC_ 'OGM'

MODULE aed_bivalve
!-------------------------------------------------------------------------------
!  aed_bivalve --- multi-group bivalve biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util
!  USE aed_zoop_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_bivalve_data_t
!
!  %% NAMELIST   %%  type : bivalve_prey_t
   TYPE bivalve_prey_t
      !State variable name for bivalvelankton prey
      CHARACTER(64) :: bivalve_prey
      !Preference factors for bivalvelankton predators grazing on prey
      AED_REAL      :: Pbiv_prey
   END TYPE bivalve_prey_t
!  %% END NAMELIST   %%  type : bivalve_prey_t


!  %% NAMELIST   %% type : bivalve_params_t
   TYPE bivalve_params_t
      ! General Attributes
      CHARACTER(64) :: name
      AED_REAL :: initial_conc, min
      INTEGER  :: Length
      AED_REAL :: X_ci

      ! Nutrient parameters
       AED_REAL :: INC, IPC

      ! Growth rate parameters
      AED_REAL :: Rgrz
      AED_REAL :: Ing, WaI, WbI
      AED_REAL :: fassim

      ! Minumum prey concentration parameters
      AED_REAL :: Cmin_grz, Kgrz
      AED_REAL :: minT, Tmin, Tmax, maxT, Dmax, maxD, SSmax, maxSS

      ! Respiration, mortaility and excretion parameters
      AED_REAL :: Rexcr, Regst, gegst, Rresp

      ! Salinity parameters
      INTEGER  :: saltfunc
      AED_REAL :: minS, Smin, Smax, maxS

      AED_REAL :: fR20, War, Wbr, fR

      AED_REAL :: theta_resp
      AED_REAL :: TmaxR, maxTR, Qresp

      AED_REAL :: SDA, Rmort, Rpred
      AED_REAL :: fDO, K_BDO, KDO

      ! The prey
      INTEGER  :: num_prey
      TYPE(bivalve_prey_t) :: prey(MAX_ZOOP_PREY)
   END TYPE
!  %% END NAMELIST   %% type : bivalve_params_t

   TYPE,extends(bivalve_params_t) :: bivalve_data_t
      INTEGER  :: id_prey(MAX_ZOOP_PREY)
      INTEGER  :: id_phyIN(MAX_ZOOP_PREY), id_phyIP(MAX_ZOOP_PREY)
   END TYPE


   TYPE,extends(aed_model_data_t) :: aed_bivalve_data_t
      !# Variable identifiers
      INTEGER  :: id_biv(MAX_ZOOP_TYPES)
      INTEGER  :: id_Nexctarget,id_Nmorttarget
      INTEGER  :: id_Pexctarget,id_Pmorttarget
      INTEGER  :: id_Cexctarget,id_Cmorttarget
      INTEGER  :: id_DOupttarget,id_SSupttarget
      INTEGER  :: id_tem,id_sal,id_sed_zone
      INTEGER  :: id_grz,id_resp,id_mort,id_excr,id_egst
      INTEGER  :: id_excr_n,id_excr_p,id_egst_n,id_egst_p,id_grz_n,id_grz_p
      INTEGER  :: id_tbiv,id_nmp,id_fT,id_fD,id_fG
      INTEGER  :: id_bivtr,id_pf,id_FR
      INTEGER  :: id_3d_grz,id_3d_exn,id_3d_exp
      INTEGER  :: id_bnum,id_x_c

      !# Model parameters
      INTEGER  :: num_biv
      TYPE(bivalve_data_t),DIMENSION(:),ALLOCATABLE :: bivalves
      LOGICAL  :: simDNexcr, simDPexcr, simDCexcr
      LOGICAL  :: simPNexcr, simPPexcr, simPCexcr
      LOGICAL  :: simSSlim, simBivFeedback, simFixedEnv, simStaticBiomass
      LOGICAL  :: initFromDensity
      INTEGER  :: n_zones
      AED_REAL,ALLOCATABLE :: active_zones(:)
      AED_REAL :: fixed_temp, fixed_sal, fixed_oxy, fixed_food
      AED_REAL :: bt_renewal

     CONTAINS
         PROCEDURE :: define             => aed_define_bivalve
         PROCEDURE :: initialize_benthic => aed_initialize_benthic_bivalve
         PROCEDURE :: calculate_benthic  => aed_calculate_benthic_bivalve
        !PROCEDURE :: mobility           => aed_mobility_bivalve
        !PROCEDURE :: light_extinction   => aed_light_extinction_bivalve
        !PROCEDURE :: delete             => aed_delete_bivalve

   END TYPE

   ! MODULE GLOBALS
   LOGICAL :: extra_diag = .false.
   INTEGER :: diag_level = 10                 ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

CONTAINS
!===============================================================================


!###############################################################################
INTEGER FUNCTION load_csv(dbase,bivalve_param)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(bivalve_params_t),INTENT(out) :: bivalve_param(MAX_ZOOP_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0
   LOGICAL :: meh
   INTEGER :: ret = 0
!
!BEGIN
!-------------------------------------------------------------------------------
   unit = aed_csv_read_header(dbase, csvnames, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   ALLOCATE(values(nccols))

   DO WHILE ( aed_csv_read_row(unit, values) )
      DO ccol=2,nccols
         bivalve_param(ccol)%name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('initial_conc')         ; bivalve_param(ccol)%initial_conc = extract_double(values(ccol))
            CASE ('min')                  ; bivalve_param(ccol)%min          = extract_double(values(ccol))
            CASE ('length')               ; bivalve_param(ccol)%length       = extract_integer(values(ccol))
            CASE ('INC')                  ; bivalve_param(ccol)%INC          = extract_double(values(ccol))
            CASE ('IPC')                  ; bivalve_param(ccol)%IPC          = extract_double(values(ccol))
            CASE ('Rgrz')                 ; bivalve_param(ccol)%Rgrz         = extract_double(values(ccol))
            CASE ('Ing')                  ; bivalve_param(ccol)%Ing          = extract_integer(values(ccol))
            CASE ('WaI')                  ; bivalve_param(ccol)%WaI          = extract_double(values(ccol))
            CASE ('WbI')                  ; bivalve_param(ccol)%WbI          = extract_double(values(ccol))
            CASE ('fassim')               ; bivalve_param(ccol)%fassim       = extract_double(values(ccol))
            CASE ('Cmin_grz')             ; bivalve_param(ccol)%Cmin_grz     = extract_double(values(ccol))
            CASE ('Kgrz')                 ; bivalve_param(ccol)%Kgrz         = extract_double(values(ccol))
            CASE ('minT')                 ; bivalve_param(ccol)%minT         = extract_double(values(ccol))
            CASE ('Tmin')                 ; bivalve_param(ccol)%Tmin         = extract_double(values(ccol))
            CASE ('Tmax')                 ; bivalve_param(ccol)%Tmax         = extract_double(values(ccol))
            CASE ('maxT')                 ; bivalve_param(ccol)%maxT         = extract_double(values(ccol))
            CASE ('Dmax')                 ; bivalve_param(ccol)%Dmax         = extract_double(values(ccol))
            CASE ('maxD')                 ; bivalve_param(ccol)%maxD         = extract_double(values(ccol))
            CASE ('SSmax')                ; bivalve_param(ccol)%SSmax        = extract_double(values(ccol))
            CASE ('maxSS')                ; bivalve_param(ccol)%maxSS        = extract_double(values(ccol))
            CASE ('Rexcr')                ; bivalve_param(ccol)%Rexcr        = extract_double(values(ccol))
            CASE ('Regst')                ; bivalve_param(ccol)%Regst        = extract_double(values(ccol))
            CASE ('gegst')                ; bivalve_param(ccol)%gegst        = extract_double(values(ccol))
            CASE ('Rresp')                ; bivalve_param(ccol)%Rresp        = extract_double(values(ccol))
            CASE ('saltfunc')             ; bivalve_param(ccol)%saltfunc     = extract_integer(values(ccol))
            CASE ('minS')                 ; bivalve_param(ccol)%minS         = extract_double(values(ccol))
            CASE ('Smin')                 ; bivalve_param(ccol)%Smin         = extract_double(values(ccol))
            CASE ('Smax')                 ; bivalve_param(ccol)%Smax         = extract_double(values(ccol))
            CASE ('maxS')                 ; bivalve_param(ccol)%maxS         = extract_double(values(ccol))
            CASE ('fR20')                 ; bivalve_param(ccol)%fR20         = extract_integer(values(ccol))
            CASE ('War')                  ; bivalve_param(ccol)%War          = extract_double(values(ccol))
            CASE ('Wbr')                  ; bivalve_param(ccol)%Wbr          = extract_double(values(ccol))
            CASE ('fR')                   ; bivalve_param(ccol)%fR           = extract_integer(values(ccol))
            CASE ('theta_resp')           ; bivalve_param(ccol)%theta_resp   = extract_double(values(ccol))
            CASE ('TmaxR')                ; bivalve_param(ccol)%TmaxR        = extract_double(values(ccol))
            CASE ('maxTR')                ; bivalve_param(ccol)%maxTR        = extract_double(values(ccol))
            CASE ('Qresp')                ; bivalve_param(ccol)%Qresp        = extract_double(values(ccol))
            CASE ('SDA')                  ; bivalve_param(ccol)%SDA          = extract_double(values(ccol))
            CASE ('Rmort')                ; bivalve_param(ccol)%Rmort        = extract_double(values(ccol))
            CASE ('Rpred')                ; bivalve_param(ccol)%Rpred        = extract_double(values(ccol))
            CASE ('fDO')                  ; bivalve_param(ccol)%fDO          = extract_integer(values(ccol))
            CASE ('K_BDO')                ; bivalve_param(ccol)%K_BDO        = extract_double(values(ccol))
            CASE ('KDO')                  ; bivalve_param(ccol)%KDO          = extract_double(values(ccol))
            CASE ('num_prey')             ; bivalve_param(ccol)%num_prey     = extract_integer(values(ccol))
            CASE ('prey(1)%bivalve_prey') ; CALL copy_name(values(ccol), bivalve_param(ccol)%prey(1)%bivalve_prey)
            CASE ('prey(1)%Pbiv_prey')    ; bivalve_param(ccol)%prey(1)%Pbiv_prey = extract_double(values(ccol))
            CASE ('prey(2)%bivalve_prey') ; CALL copy_name(values(ccol), bivalve_param(ccol)%prey(3)%bivalve_prey)
            CASE ('prey(2)%Pbiv_prey')    ; bivalve_param(ccol)%prey(2)%Pbiv_prey = extract_double(values(ccol))
            CASE ('prey(3)%bivalve_prey') ; CALL copy_name(values(ccol), bivalve_param(ccol)%prey(3)%bivalve_prey)
            CASE ('prey(3)%Pbiv_prey')    ; bivalve_param(ccol)%prey(3)%Pbiv_prey = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit)
   !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bivalve_load_params(data, dbase, count, list, X_c)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count   !Number of bivalve groups
   INTEGER,INTENT(in)          :: list(*) !List of bivalve groups to simulate
   AED_REAL,INTENT(in)         :: X_c(*)  !
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   AED_REAL :: Pbiv_prey(MAX_ZOOP_PREY)

   TYPE(bivalve_params_t)  :: bivalve_param(MAX_ZOOP_TYPES)
   NAMELIST /bivalve_params/ bivalve_param   ! %% bivalve_params_t - see above
!-------------------------------------------------------------------------------
!BEGIN
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, bivalve_param)
       CASE (NML_TYPE)
           tfil = find_free_lun()
           open(tfil,file=dbase, status='OLD',iostat=status)
           IF (status /= 0) STOP 'Error opening bivalves_params namelist file'
           read(tfil,nml=bivalve_params,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist bivalves_params'

    data%num_biv = count
    allocate(data%bivalves(count))
    DO i=1,count
       ! General
       data%bivalves(i)%name          = bivalve_param(list(i))%name
       data%bivalves(i)%initial_conc  = bivalve_param(list(i))%initial_conc
       data%bivalves(i)%min           = bivalve_param(list(i))%min
       data%bivalves(i)%Length        = bivalve_param(list(i))%Length
       data%bivalves(i)%X_ci          = X_c(i) ! 1. ! ADD DETAIL
       data%bivalves(i)%INC           = bivalve_param(list(i))%INC
       data%bivalves(i)%IPC           = bivalve_param(list(i))%IPC
       ! Filtration & Ingestion
       data%bivalves(i)%Rgrz          = bivalve_param(list(i))%Rgrz/secs_per_day
       data%bivalves(i)%Ing           = bivalve_param(list(i))%Ing
       data%bivalves(i)%WaI           = bivalve_param(list(i))%WaI/secs_per_day
       data%bivalves(i)%WbI           = bivalve_param(list(i))%WbI
       data%bivalves(i)%fassim        = bivalve_param(list(i))%fassim
       data%bivalves(i)%Cmin_grz      = bivalve_param(list(i))%Cmin_grz
       data%bivalves(i)%Kgrz          = bivalve_param(list(i))%Kgrz
       data%bivalves(i)%minT          = bivalve_param(list(i))%minT
       data%bivalves(i)%Tmin          = bivalve_param(list(i))%Tmin
       data%bivalves(i)%Tmax          = bivalve_param(list(i))%Tmax
       data%bivalves(i)%maxT          = bivalve_param(list(i))%maxT
       data%bivalves(i)%Dmax          = bivalve_param(list(i))%Dmax
       data%bivalves(i)%maxD          = bivalve_param(list(i))%maxD
       data%bivalves(i)%SSmax         = bivalve_param(list(i))%SSmax
       data%bivalves(i)%maxSS         = bivalve_param(list(i))%maxSS
       ! Excretion & Egestion
       data%bivalves(i)%Rexcr         = bivalve_param(list(i))%Rexcr
       data%bivalves(i)%Regst         = bivalve_param(list(i))%Regst
       data%bivalves(i)%gegst         = bivalve_param(list(i))%gegst
       ! Respiration
       data%bivalves(i)%Rresp         = bivalve_param(list(i))%Rresp/secs_per_day
       data%bivalves(i)%saltfunc      = bivalve_param(list(i))%saltfunc
       data%bivalves(i)%minS          = bivalve_param(list(i))%minS
       data%bivalves(i)%Smin          = bivalve_param(list(i))%Smin
       data%bivalves(i)%Smax          = bivalve_param(list(i))%Smax
       data%bivalves(i)%maxS          = bivalve_param(list(i))%maxS
       data%bivalves(i)%fR20          = bivalve_param(list(i))%fR20
       data%bivalves(i)%War           = bivalve_param(list(i))%War
       data%bivalves(i)%Wbr           = bivalve_param(list(i))%Wbr
       data%bivalves(i)%fR            = bivalve_param(list(i))%fR
       data%bivalves(i)%theta_resp    = bivalve_param(list(i))%theta_resp
       data%bivalves(i)%TmaxR         = bivalve_param(list(i))%TmaxR
       data%bivalves(i)%maxTR         = bivalve_param(list(i))%maxTR
       data%bivalves(i)%Qresp         = bivalve_param(list(i))%Qresp
       data%bivalves(i)%SDA           = bivalve_param(list(i))%SDA
       ! Mortality
       data%bivalves(i)%Rmort         = bivalve_param(list(i))%Rmort/secs_per_day
       data%bivalves(i)%Rpred         = bivalve_param(list(i))%Rpred/secs_per_day
       data%bivalves(i)%fDO           = bivalve_param(list(i))%fDO
       data%bivalves(i)%K_BDO         = bivalve_param(list(i))%K_BDO
       data%bivalves(i)%KDO           = bivalve_param(list(i))%KDO

       data%bivalves(i)%num_prey      = bivalve_param(list(i))%num_prey

       !Loop through prey variables assigning a target variable and preference factor
       !First sort in descending order of food preferences
       DO j=1,data%bivalves(i)%num_prey
          sort_i(j) = j
          Pbiv_prey(j) = bivalve_param(list(i))%prey(j)%Pbiv_prey
       ENDDO
       CALL qsort(Pbiv_prey,sort_i,1,data%bivalves(i)%num_prey)
       DO j=1,data%bivalves(i)%num_prey
          data%bivalves(i)%prey(j)%bivalve_prey = &
                   bivalve_param(list(i))%prey(sort_i(data%bivalves(i)%num_prey-j+1))%bivalve_prey
          data%bivalves(i)%prey(j)%Pbiv_prey = &
                   bivalve_param(list(i))%prey(sort_i(data%bivalves(i)%num_prey-j+1))%Pbiv_prey
       ENDDO

       ! Register group as a state variable
       data%id_biv(i) = aed_define_sheet_variable(                 &
                              bivalve_param(list(i))%name,          &
                              'mmolC/m**2', 'bivalve',              &
                              bivalve_param(list(i))%initial_conc,  &
                              minimum=bivalve_param(list(i))%min)
    ENDDO
!
END SUBROUTINE aed_bivalve_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_bivalve(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the bivalve biogeochemical model
!
!  Here, the aed_bivalve namelist is read and te variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status, i

!  %% NAMELIST   %%  /aed_bivalve/
!  %% Last Checked 20/08/2021
   INTEGER  :: num_biv = 0
   INTEGER  :: the_biv(MAX_ZOOP_TYPES) = 0
   INTEGER  :: n_zones = 0
   INTEGER  :: active_zones(MAX_ZONES)
   LOGICAL  :: initFromDensity = .false.  !ability to read in map of mussel #'s
   LOGICAL  :: simBivTracer = .false.     !include filtration rate tracer
   LOGICAL  :: simBivFeedback = .false.   !allow module to change prey/nutrient concs
   LOGICAL  :: simStaticBiomass = .false. !keep biomass fixed over time
   LOGICAL  :: simFixedEnv = .false.      !special case to overwrite env and food factors with constants
   AED_REAL :: X_c(MAX_ZOOP_TYPES) = 1.   !mmolC/organism
   AED_REAL :: bt_renewal = zero_         !biv-tracer "renewal" timescale (days)
   AED_REAL :: fixed_temp
   AED_REAL :: fixed_sal
   AED_REAL :: fixed_oxy
   AED_REAL :: fixed_food

   CHARACTER(len=64)  :: dn_target_variable = '' !dissolved nitrogen target variable
   CHARACTER(len=64)  :: pn_target_variable = '' !particulate nitrogen target variable
   CHARACTER(len=64)  :: dp_target_variable = '' !dissolved phosphorus target variable
   CHARACTER(len=64)  :: pp_target_variable = '' !particulate phosphorus target variable
   CHARACTER(len=64)  :: dc_target_variable = '' !dissolved carbon target variable
   CHARACTER(len=64)  :: pc_target_variable = '' !particulate carbon target variable
   CHARACTER(len=64)  :: do_uptake_variable = '' !oxy uptake variable
   CHARACTER(len=64)  :: ss_uptake_variable = '' !sus. solids uptake variable
   CHARACTER(len=128) :: dbase = 'aed_bivalve_pars.nml'

! %% From Module Globals
!  LOGICAL  :: extra_diag = .false.      !## Obsolete Use diag_level = 10
!  INTEGER :: diag_level = 10                 ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_bivalve/

   INTEGER  :: biv_i, prey_i, phy_i

   NAMELIST /aed_bivalve/ num_biv, the_biv, n_zones, active_zones, extra_diag, &
                    simBivTracer, simBivFeedback, simStaticBiomass,            &
                    dn_target_variable, pn_target_variable, dp_target_variable,&
                    pp_target_variable, dc_target_variable, pc_target_variable,&
                    do_uptake_variable, ss_uptake_variable, dbase,             &
                    simFixedEnv, fixed_temp, fixed_sal, fixed_oxy, fixed_food, &
                    initFromDensity, diag_level, X_c, bt_renewal
!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_bivalve initialization"

   ! Read the namelist
   read(namlst,nml=aed_bivalve,iostat=status)
   IF (status /= 0) STOP 'ERROR reading namelist aed_bivalve'

   IF( extra_diag )   diag_level = 10           ! legacy use of extra_debug

   ! Assign module level variables
   data%simFixedEnv = simFixedEnv
   IF (simFixedEnv) THEN
     PRINT *,'    NOTE - bivalves driven by fixed T,DO,A values'
     PRINT *,'         - biv feedback disabled'
     data%fixed_temp = fixed_temp
     data%fixed_sal = fixed_sal
     data%fixed_oxy = fixed_oxy
     data%fixed_food = fixed_food
     simBivFeedback = .false.
   END IF

   data%simBivFeedback = simBivFeedback
   PRINT *,'    NOTE - bivalve feedbacks on water column properties: ',simBivFeedback

   data%bt_renewal = bt_renewal
   PRINT *,'    NOTE - bivalve filtration tracer variable: ',simBivTracer

   data%initFromDensity = initFromDensity
   PRINT *,'    NOTE - bivalve mass to be initialised from initial density map: ',initFromDensity

   data%simStaticBiomass = simStaticBiomass
   PRINT *,'    NOTE - bivalve mass to be held static: ',simStaticBiomass

   ! Configure which sediment zones to run within
   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
      PRINT *,'    NOTE - bivalves growing in these zones: ',data%active_zones
   ELSE
      PRINT *,'    NOTE - bivalves growing in all zones'
   ENDIF

   ! Register variables and store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed_bivalve_load_params to values per second.
   data%num_biv = 0
   CALL aed_bivalve_load_params(data, dbase, num_biv, the_biv, X_c)


   !Register link to prey state variables
   DO biv_i = 1,num_biv
      phy_i = 0
      DO prey_i = 1,data%bivalves(biv_i)%num_prey
          data%bivalves(biv_i)%id_prey(prey_i) = aed_locate_variable( &
                                       data%bivalves(biv_i)%prey(prey_i)%bivalve_prey)
          !If the prey is phytoplankton then also register state dependency on
          !internal nitrogen and phosphorus
          IF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_MODPFX_).EQ. _PHYMOD_) THEN
              phy_i = phy_i + 1
              data%bivalves(biv_i)%id_phyIN(phy_i) = aed_locate_variable( &
                                       TRIM(data%bivalves(biv_i)%prey(prey_i)%bivalve_prey)//'_IN')
              data%bivalves(biv_i)%id_phyIP(phy_i) = aed_locate_variable( &
                                       TRIM(data%bivalves(biv_i)%prey(prey_i)%bivalve_prey)//'_IP')
          ENDIF
      ENDDO
   ENDDO

   ! Register link to nutrient pools, if variable names are provided in namelist
   !  and feedbacks are on
   data%simDNexcr = dn_target_variable .NE. ''
   IF (data%simDNexcr) THEN
     data%id_Nexctarget = aed_locate_variable(dn_target_variable)
   ENDIF
   data%simDPexcr = dp_target_variable .NE. ''
   IF (data%simDPexcr) THEN
     data%id_Pexctarget = aed_locate_variable(dp_target_variable)
   ENDIF
   data%simDCexcr = dc_target_variable .NE. ''
   IF (data%simDCexcr) THEN
     data%id_Cexctarget = aed_locate_variable(dc_target_variable)
   ENDIF

   data%simPNexcr = pn_target_variable .NE. ''
   IF (data%simPNexcr) THEN
     data%id_Nmorttarget = aed_locate_variable(pn_target_variable)
   ENDIF
   data%simPPexcr = pp_target_variable .NE. ''
   IF (data%simPPexcr) THEN
     data%id_Pmorttarget = aed_locate_variable(pp_target_variable)
   ENDIF
   data%simPCexcr = pc_target_variable .NE. ''
   IF (data%simPCexcr) THEN
     data%id_Cmorttarget = aed_locate_variable(pc_target_variable)
   ENDIF

   IF (do_uptake_variable .EQ. '') STOP 'bivalve needs DO uptake variable'
   data%id_DOupttarget = aed_locate_variable(do_uptake_variable)

   data%simSSlim = ss_uptake_variable .NE. ''
   IF (data%simSSlim) THEN
     data%id_SSupttarget = aed_locate_variable(ss_uptake_variable)
   ENDIF

   IF (simBivTracer) THEN
       ! Register group as a state variable
       data%id_bivtr = aed_define_variable(                                   &
                              'filtfrac',                                      &
                              '-', 'fraction of water filtered by bivalves',   &
                              ZERO_,  &
                              minimum=ZERO_)
   ENDIF

   ! Register diagnostic variables
   IF (initFromDensity) THEN
     data%id_bnum = aed_define_sheet_diag_variable('num' ,'#/m**2','mussel density')
     IF (diag_level>0) data%id_x_c = aed_define_sheet_diag_variable('X_c','mmolC/organism','mean bivalve mass / organism')
   ENDIF
   IF (diag_level>0) data%id_tbiv = aed_define_sheet_diag_variable('tbiv','mmolC/m**2','total bivalve mass')
   IF (diag_level>1) data%id_nmp  = aed_define_sheet_diag_variable('nmp' ,'mmolC/m**2/d','net mussel production')
   IF (diag_level>1) data%id_resp = aed_define_sheet_diag_variable('resp','/d','bivalve respiration')
   IF (diag_level>2) THEN
     data%id_3d_grz = aed_define_diag_variable('net_c','mmolC/m**3/day','water column change due to grazing/excretion')
     data%id_3d_exn = aed_define_diag_variable('net_n','mmolN/m**3/day','water column change due to grazing/excretion')
     data%id_3d_exp = aed_define_diag_variable('net_p','mmolP/m**3/day','water column change due to grazing/excretion')
   ENDIF
   IF (diag_level>9) THEN
     data%id_FR = aed_define_sheet_diag_variable('FR','m3/mmolC/m2/day','filtration rate')
     data%id_fT = aed_define_sheet_diag_variable('fT','-','temp limitation')
     data%id_fD = aed_define_sheet_diag_variable('fD','-','density limitation')
     data%id_fG = aed_define_sheet_diag_variable('fG','-','grazing limitation')
     data%id_pf = aed_define_sheet_diag_variable('pf','/d','pseudofeaces production rate')
     data%id_grz  = aed_define_sheet_diag_variable('grz' ,'/d', 'bivalve grazing')
     data%id_mort = aed_define_sheet_diag_variable('mort','/d','bivalve mortality')
     data%id_egst = aed_define_sheet_diag_variable('egst','/d','bivalve egestion')
     data%id_excr = aed_define_sheet_diag_variable('excr','/d','bivalve excretion')
     data%id_excr_n = aed_define_sheet_diag_variable('excr_n','/d','bivalve excretion')
     data%id_excr_p = aed_define_sheet_diag_variable('excr_p','/d','bivalve excretion')
     data%id_egst_n = aed_define_sheet_diag_variable('egst_n','/d','bivalve excretion')
     data%id_egst_p = aed_define_sheet_diag_variable('egst_p','/d','bivalve excretion')
     data%id_grz_n = aed_define_sheet_diag_variable('grz_n','/d','bivalve excretion')
     data%id_grz_p = aed_define_sheet_diag_variable('grz_p','/d','bivalve excretion')
   ENDIF

   ! Register environmental dependencies
   data%id_tem = aed_locate_global('temperature')
   data%id_sal = aed_locate_global('salinity')
   data%id_sed_zone = aed_locate_sheet_global('sed_zone')
!
END SUBROUTINE aed_define_bivalve
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic_bivalve(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to initialize bivalve mass from spatial density map
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: biv_i
   AED_REAL :: matz, biv_density
!-------------------------------------------------------------------------------
!BEGIN
   ! Check user wishes to initialise from spatial map of biv density
   IF ( .NOT. data%initFromDensity ) RETURN

   ! Check to ensure this zone is colonisable
   matz = _STATE_VAR_S_(data%id_sed_zone)
   IF ( .NOT. in_zone_set(matz, data%active_zones) ) RETURN

   ! Get the bivalve density (orgs/m2)
   biv_density = _DIAG_VAR_S_(data%id_bnum)

   ! Distribute the read in biv density over the simulated groups
   DO biv_i=1,data%num_biv
     _STATE_VAR_S_(data%id_biv(biv_i)) = &
                            (biv_density/data%num_biv)*data%bivalves(biv_i)%X_ci
   ENDDO
END SUBROUTINE aed_initialize_benthic_bivalve
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_bivalve(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of bivalve biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: biv, temp, salinity, oxy, ss ,matz !State variables
   AED_REAL :: prey(MAX_ZOOP_PREY), grazing_prey(MAX_ZOOP_PREY) !Prey state variables
   AED_REAL :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for phytos
   AED_REAL :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL :: pon, pop, poc !Mortaility and fecal pellet state variables
   AED_REAL :: FGrazing_Limitation, f_Temp, f_Salinity, f_SS, I_max
   AED_REAL :: pref_factor, Ctotal_prey !total concentration of available prey
   AED_REAL :: food, grazing, respiration, mortality !Growth & decay functions
   AED_REAL :: grazing_n, grazing_p, gn, gp !Grazing on nutrients
   AED_REAL :: pon_excr, pop_excr, poc_excr !POM excretion rates
   AED_REAL :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   AED_REAL :: f_Dens, W, Imax, psuedofaeces, ingestion, excretion, egestion, iteg, R20
   AED_REAL :: bt, fr, Rbt !BivTracer vars
   INTEGER  :: biv_i,prey_i,prey_j,phy_i !Counters
!
!-------------------------------------------------------------------------------
!BEGIN
   bt = 0.0 ; poc = 0.0 ; pon = 0.0 ; pop = 0.0 ; Ctotal_prey = 0.0   !## CAB [-Wmaybe-uninitialized]

   ! Check to ensure this zone is colonisable, and numbers are adequate
   matz = _STATE_VAR_S_(data%id_sed_zone)

   IF ( .NOT. in_zone_set(matz, data%active_zones) ) RETURN
   IF ( .NOT. data%initFromDensity ) THEN
      IF ( _DIAG_VAR_S_(data%id_bnum) <1e-3 ) RETURN
   ENDIF


   ! Retrieve current environmental conditions.
   IF(data%simFixedEnv) THEN
     temp = data%fixed_temp + matz             ! user provided temp
     salinity = data%fixed_sal                 ! user provided salinity
     oxy = data%fixed_oxy                      ! user provided oxygen
     Ctotal_prey = data%fixed_food             ! user provided food
     ss = zero_
   ELSE
     ! Normal operation
     temp = _STATE_VAR_(data%id_tem)           ! local temperature
     salinity = _STATE_VAR_(data%id_sal)       ! local salinity
     oxy = _STATE_VAR_(data%id_DOupttarget)    ! local oxygen
     IF (data%simSSlim) THEN
        ss = _STATE_VAR_(data%id_SSupttarget)  ! local suspended solids (inorganic)
     ELSE
        ss = zero_
     ENDIF
   ENDIF

   ! Retrieve current (local) state variable values.
   IF (data%simDNexcr)  dn_excr = _STATE_VAR_(data%id_Nexctarget)
   IF (data%simDPexcr)  dp_excr = _STATE_VAR_(data%id_Pexctarget)
   IF (data%simDCexcr)  dc_excr = _STATE_VAR_(data%id_Cexctarget)

   IF (data%simPNexcr)  pon = _STATE_VAR_(data%id_Nmorttarget)
   IF (data%simPPexcr)  pop = _STATE_VAR_(data%id_Pmorttarget)
   IF (data%simPCexcr)  poc = _STATE_VAR_(data%id_Cmorttarget)

   ! Zero cumulative diagnostics
   IF(diag_level>0) _DIAG_VAR_S_(data%id_tbiv) = zero_
   IF(diag_level>2) _DIAG_VAR_(data%id_3d_grz) = zero_
   IF(diag_level>2) _DIAG_VAR_(data%id_3d_exn) = zero_
   IF(diag_level>2) _DIAG_VAR_(data%id_3d_exp) = zero_
   IF(diag_level>9) _DIAG_VAR_S_(data%id_grz_n)= zero_
   IF(diag_level>9) _DIAG_VAR_S_(data%id_grz_p)= zero_
   IF(diag_level>9) _DIAG_VAR_S_(data%id_excr_n)= zero_
   IF(diag_level>9) _DIAG_VAR_S_(data%id_excr_p)= zero_
   IF(diag_level>9) _DIAG_VAR_S_(data%id_egst_n)= zero_
   IF(diag_level>9) _DIAG_VAR_S_(data%id_egst_p)= zero_


   ! Setup for bivalve filtration "tracer"
   Rbt = data%bt_renewal/secs_per_day       ! biv_tracer decay (%/day)
   fr = zero_                               ! filtration rate (m3/mmolC/day)
   IF (data%id_bivtr>0) THEN
     bt = _STATE_VAR_(data%id_bivtr)
   ENDIF
   Rbt = Rbt*bt

   ! Loop through number of configured groups
   DO biv_i=1,data%num_biv

      ! Retrieve this bivalve group biomass (mmolC/m2)
      biv = MAX(_STATE_VAR_S_(data%id_biv(biv_i)),zero_)
      !_DIAG_VAR_S_(data%id_pbiv) = _DIAG_VAR_S_(data%id_pbiv) + biv

      grazing       = zero_
      respiration   = zero_
      mortality     = zero_
      don_excr      = zero_
      dop_excr      = zero_
      doc_excr      = zero_


      ! Retrieve prey groups (mmolC/m3)
      IF(.NOT.data%simFixedEnv) THEN
        Ctotal_prey   = zero_
        DO prey_i=1,data%bivalves(biv_i)%num_prey
           prey(prey_i) = _STATE_VAR_(data%bivalves(biv_i)%id_prey(prey_i))
           Ctotal_prey = Ctotal_prey + prey(prey_i)
        ENDDO
      ENDIF

      ! Get the grazing limitation function
      !fGrazing_Limitation = fPrey_Limitation(data%bivalves,biv_i,Ctotal_prey)
      fGrazing_Limitation = min(Ctotal_prey/data%bivalves(biv_i)%Kgrz,one_)


      ! Get the temperature function impacting ingestion
      f_Temp = fTemp_function_biv(data,biv_i, temp)

      ! Get the suspended solids function
      f_SS = fSS_function(data,biv_i,ss)

      ! Get the density limitation function
      f_Dens = fD_function(data,biv_i,biv)

      ! Now compute the final grazing rate (/s). This is the
      ! amount grazed in units of mass consumed/mass bivalve/unit time
      IF(data%bivalves(biv_i)%Ing==1) THEN
        W = (0.071/1000.) * data%bivalves(biv_i)%Length**2.8
        Imax = data%bivalves(biv_i)%WaI * W** data%bivalves(biv_i)%WbI
      ELSE
        Imax = data%bivalves(biv_i)%Rgrz
      END IF
      grazing = Imax * fGrazing_Limitation * f_Temp * f_Dens * f_SS

      IF(Ctotal_prey<data%bivalves(biv_i)%Kgrz) THEN
        FR = grazing / data%bivalves(biv_i)%Kgrz
      ELSE
        FR = grazing / Ctotal_prey
      END IF

      ! Now determine available prey and limit grazing amount to availability of prey
      ! food is total amount of food in units of mass/unit volume/unit time (mmolC/m2/s)
      food = grazing * biv
      IF (Ctotal_prey < data%bivalves(biv_i)%num_prey * data%bivalves(biv_i)%Cmin_grz ) THEN
         food = zero_
         grazing = zero_
      ELSEIF (food > Ctotal_prey - data%bivalves(biv_i)%num_prey * data%bivalves(biv_i)%Cmin_grz ) THEN
         food = Ctotal_prey - data%bivalves(biv_i)%num_prey * data%bivalves(biv_i)%Cmin_grz
         grazing = food / MAX(biv,1e-2)  !/s
      ENDIF

      ! Now determine prey composition based on preference factors and availability of prey

      ! Prey has been ordered in grazing preference;  So take food in the
      ! order of preference up to availability minus value of minimum residual
      ! note: grazing_prey is in units of mass consumed/unit volume/unit time
      IF(.NOT.data%simFixedEnv) THEN
       DO prey_i = 1,data%bivalves(biv_i)%num_prey
          ! Add up preferences for remaining prey
          pref_factor = zero_
          DO prey_j = prey_i,data%bivalves(biv_i)%num_prey
             pref_factor = pref_factor + data%bivalves(biv_i)%prey(biv_i)%Pbiv_prey
          ENDDO
          IF (food * data%bivalves(biv_i)%prey(prey_i)%Pbiv_prey / pref_factor <= &
                              prey(prey_i) - data%bivalves(biv_i)%Cmin_grz) THEN
             ! Take fraction of left over food based on preference factor
             grazing_prey(prey_i) = food * data%bivalves(biv_i)%prey(prey_i)%Pbiv_prey / pref_factor
          ELSEIF (prey(prey_i) > data%bivalves(biv_i)%Cmin_grz) THEN
             grazing_prey(prey_i) = prey(prey_i) - data%bivalves(biv_i)%Cmin_grz
          ELSE
             grazing_prey(prey_i) = zero_
          ENDIF
          ! Food remaining after grazing from current food group (prey)
          food = food - grazing_prey(prey_i)
       ENDDO
       ! Reset food to total amount consumed
       food = SUM(grazing_prey(:))
      ENDIF

      ! Now determine nutrient composition of food based on mix of food types.
      ! note:
      !   grazing_n is in units of mass N consumed/unit volume/unit time
      !   grazing_p is in units of mass P consumed/unit volume/unit time
      !   these are the same units as "food" & "grazing_prey", but not "grazing"
      IF(data%simFixedEnv) THEN
        grazing_n = food * 1./5.5   !grazing * 16./106.
        grazing_p = food * 1./130.  !grazing * 1./106.
      ELSE
        grazing_n = zero_
        grazing_p = zero_
        phy_i = 0
        DO prey_i = 1,data%bivalves(biv_i)%num_prey
         IF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_MODPFX_) .EQ. _OGMPOC_) THEN
            IF (poc > 1e-6) THEN
                grazing_n = grazing_n + grazing_prey(prey_i) *  MIN(pon/poc,2.) !pon/poc
                grazing_p = grazing_p + grazing_prey(prey_i) *  MIN(pop/poc,2.) !pop/poc
            ELSE
                grazing_n = zero_
                grazing_p = zero_
            ENDIF
         ELSEIF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_MODPFX_).EQ. _PHYMOD_) THEN
            phy_i = phy_i + 1
            phy_INcon(phy_i) = _STATE_VAR_(data%bivalves(biv_i)%id_phyIN(phy_i))
            phy_IPcon(phy_i) = _STATE_VAR_(data%bivalves(biv_i)%id_phyIP(phy_i))
            grazing_n = grazing_n + grazing_prey(prey_i) / MAX(prey(prey_i),1e-6) * phy_INcon(phy_i) !/14.0
            grazing_p = grazing_p + grazing_prey(prey_i) / MAX(prey(prey_i),1e-6) * phy_IPcon(phy_i) !/31.0
         ELSEIF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_MODPFX_).EQ.'ZOO') THEN
            grazing_n = grazing_n + grazing_prey(prey_i) * data%bivalves(biv_i)%INC
            grazing_p = grazing_p + grazing_prey(prey_i) * data%bivalves(biv_i)%IPC
         ENDIF
        ENDDO
    !    IF(grazing_n>food)print *,'gn?food',grazing_n,grazing_prey(1),pon,poc,pon/poc
      ENDIF


      ! Now compute the rate of food ingestion (/s)
      psuedofaeces = (one_ - data%bivalves(biv_i)%fassim) * grazing
      ingestion = data%bivalves(biv_i)%fassim * grazing

      IF (biv <= data%bivalves(biv_i)%min) THEN
        ! Don't excrete or die if we are at the min biomass otherwise we have a
        ! mass conservation leak in the C mass balance

        respiration = zero_
        mortality   = zero_
        excretion   = zero_
        egestion    = zero_

      ELSE
        ! Get the egestion (/s) => I * alphaEG * EXP(gammaEG * MIN([A]/KA,1))
        egestion = ingestion * data%bivalves(biv_i)%Regst * &
                   exp(data%bivalves(biv_i)%gegst * min(Ctotal_prey/data%bivalves(biv_i)%Kgrz,one_))

        ! Get the respiration rate (/s)
        iteg = ingestion - egestion
        respiration = aed_bivalve_respiration(data,biv_i,iteg,temp,salinity)

        ! Get the excretion rate (of carbon!) (/s)
        excretion =  data%bivalves(biv_i)%Rexcr * iteg

        ! Get the mortality rate (/s)
        mortality = data%bivalves(biv_i)%Rmort * f_DO(data,biv_i,oxy)

        ! Add predation losses to mortality
        mortality = mortality + data%bivalves(biv_i)%Rpred

      ENDIF

      ! Calculate losses into the particulate organic matter pools - (mmolC/m2/s)
      poc_excr = (psuedofaeces + egestion + mortality)*biv

      pon_excr = (psuedofaeces * (grazing_n / MAX(food,1e-10))   &
                   + (egestion + mortality)*data%bivalves(biv_i)%INC )*biv

!print *,'d',

      pop_excr = (psuedofaeces * (grazing_p / MAX(food,1e-10))   &
                   + (egestion + mortality)*data%bivalves(biv_i)%IPC )*biv

      ! Now we know the rates of carbon consumption and excretion, calculate
      ! the rates of n & p excretion to maintain internal nutrient stores

      ! First, compute rate of change so far of bivalve carbon (mmolC/m2/s)
      delta_C = (ingestion - respiration - egestion - excretion - mortality)*biv

      ! Then calc nutrient excretion require to balance internal nutrient store
      ! Note pon_excr includes loss due to psuedofaeces so no need to include
      ! assimilation fraction on the grazing_n & grazing_p amount

      !don_excr = grazing_n*data%bivalves(biv_i)%fassim - &
      ! ((egestion + mortality)*data%bivalves(biv_i)%INC*biv) - (delta_C * data%bivalves(biv_i)%INC)

      don_excr = grazing_n - pon_excr - (delta_C * data%bivalves(biv_i)%INC)
      dop_excr = grazing_p - pop_excr - (delta_C * data%bivalves(biv_i)%IPC)

!IF(don_excr>0) print *,'ss',grazing_n,pon_excr,(delta_C * data%bivalves(biv_i)%INC), &
!                                  food* data%bivalves(biv_i)%INC, don_excr*secs_per_day
      ! If nutrients are limiting then must excrete doc to maintain balance

      ! IF ( don_excr < zero_ .AND. dop_excr < zero_ ) THEN
      !    !Determine which nutrient is more limiting, and remove more DOC
      !
      !   ! IF((data%bivalves(biv_i)%INC * (grazing_n - pon_excr) - delta_C) .GT. &  !MH I think this should be divide by INC and INP
      !   !    (data%bivalves(biv_i)%IPC * (grazing_p - pop_excr) - delta_C)) THEN
      !
      !    IF( ((grazing_n-pon_excr)/data%bivalves(biv_i)%INC  - delta_C) > &
      !        ((grazing_p-pop_excr)/data%bivalves(biv_i)%IPC  - delta_C) ) THEN
      !        !N excr is less negative than P excr - retain N
      !        don_excr = zero_
      !        doc_excr =  (grazing_n - pon_excr) / data%bivalves(biv_i)%INC - delta_C
      !        delta_C = delta_C - doc_excr
      !        dop_excr = grazing_p - pop_excr - delta_C*data%bivalves(biv_i)%IPC
      !    ELSE
      !        !P excr is less negative than N excr - retain P
      !        dop_excr = zero_
      !        doc_excr = (grazing_p - pop_excr) / data%bivalves(biv_i)%IPC - delta_C
      !        delta_C = delta_C - doc_excr
      !        don_excr = grazing_n - pon_excr - delta_C*data%bivalves(biv_i)%INC
      !    ENDIF
      IF (don_excr < zero_ .OR. dop_excr < zero_ ) THEN
       IF (don_excr/data%bivalves(biv_i)%INC < dop_excr/data%bivalves(biv_i)%IPC) THEN !nitrogen limited
        !  print *,'don_excr1', don_excr,delta_C
          doc_excr = - don_excr / data%bivalves(biv_i)%INC
          delta_C = delta_C - doc_excr    ! reduce the growth rate due to extra C loss
          don_excr = zero_                ! hold onto the N
          dop_excr = grazing_p - pop_excr - delta_C*data%bivalves(biv_i)%IPC
        !  print *,'don_excr2', don_excr,doc_excr,delta_C,dop_excr
        ELSE
      !    print *,'dop_excr1', dop_excr,delta_C
          doc_excr = - dop_excr / data%bivalves(biv_i)%IPC
          delta_C = delta_C - doc_excr    ! reduce the growth rate due to extra C loss
          dop_excr = zero_                ! hold onto the N
          don_excr = grazing_n - pon_excr - delta_C*data%bivalves(biv_i)%INC
      !    print *,'dop_excr2', dop_excr,doc_excr,delta_C,don_excr

        ENDIF
      ELSE
        doc_excr = zero_
      ENDIF
      excretion = excretion+doc_excr/biv

      ! ELSEIF (don_excr < zero_) THEN !nitrogen limited
      !    don_excr = zero_
      !    doc_excr = (grazing_n - pon_excr) / data%bivalves(biv_i)%INC - delta_C
      !    delta_C = delta_C - doc_excr
      !    dop_excr = grazing_p - pop_excr - delta_C*data%bivalves(biv_i)%IPC
      ! ELSEIF (dop_excr < zero_) THEN !phosphorus limited
      !    dop_excr = zero_
      !    doc_excr = (grazing_p - pop_excr) / data%bivalves(biv_i)%IPC - delta_C
      !    delta_C = delta_C - doc_excr
      !    don_excr = grazing_n - pon_excr - delta_C*data%bivalves(biv_i)%INC
      ! ELSE !just excrete nutrients, no need to balance c
      !     doc_excr = zero_
      ! ENDIF

      !write(*,"(4X,'limitations (f_T,f_Salinity): ',2F8.2)")f_T,f_Salinity
      !write(*,"(4X,'sources/sinks (grazing,respiration,mortaility): ',3F8.2)")grazing,excretion,mortality


      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Biv production & losses (mmolC/m2/s)
      IF ( .NOT.data%simStaticBiomass ) THEN
        _FLUX_VAR_B_(data%id_biv(biv_i)) = _FLUX_VAR_B_(data%id_biv(biv_i))    &
          + ( (ingestion - respiration - excretion - egestion - mortality)*biv )
      ENDIF

      ! Effects on prey and nutrients
      IF (data%simBivFeedback) THEN
         ! Now take food grazed by mussels from food pools in mmol/m3/s
         ! Note : grazing_prey is in the units of food, specific for food type
         phy_i = 0
         DO prey_i = 1,data%bivalves(biv_i)%num_prey
            gn=0.;gp=0.
            _FLUX_VAR_(data%bivalves(biv_i)%id_prey(prey_i)) = &
            _FLUX_VAR_(data%bivalves(biv_i)%id_prey(prey_i)) - grazing_prey(prey_i)
            IF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_MODPFX_) .EQ. _OGMPOC_) THEN
               IF (poc > 1e-6) THEN
                  gn = ( grazing_prey(prey_i) * MIN(pon/poc,2.))
                  gp = ( grazing_prey(prey_i) * MIN(pop/poc,2.))
                  _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) - gn
                  _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) - gp
               ENDIF
            ELSEIF (data%bivalves(biv_i)%prey(prey_i)%bivalve_prey(1:_MODPFX_) .EQ. _PHYMOD_) THEN
               phy_i = phy_i + 1
               IF(prey(prey_i)>1e-6) THEN
                gn = grazing_prey(prey_i) * phy_INcon(phy_i)/prey(prey_i)
                gp = grazing_prey(prey_i) * phy_IPcon(phy_i)/prey(prey_i)
               ENDIF
               _FLUX_VAR_(data%bivalves(biv_i)%id_phyIN(phy_i)) = &
                           _FLUX_VAR_(data%bivalves(biv_i)%id_phyIN(phy_i)) - gn
               _FLUX_VAR_(data%bivalves(biv_i)%id_phyIP(phy_i)) = &
                           _FLUX_VAR_(data%bivalves(biv_i)%id_phyIP(phy_i)) - gp
            ENDIF
            IF(diag_level>2) _DIAG_VAR_(data%id_3d_grz) = _DIAG_VAR_(data%id_3d_grz)  &
                                             - grazing_prey(prey_i)*secs_per_day
            IF(diag_level>2) _DIAG_VAR_(data%id_3d_exn)=_DIAG_VAR_(data%id_3d_exn) - gn*secs_per_day
            IF(diag_level>2) _DIAG_VAR_(data%id_3d_exp)=_DIAG_VAR_(data%id_3d_exp) - gp*secs_per_day

            IF(diag_level>9) _DIAG_VAR_S_(data%id_grz_n)=_DIAG_VAR_S_(data%id_grz_n) + gn*secs_per_day
            IF(diag_level>9) _DIAG_VAR_S_(data%id_grz_p)=_DIAG_VAR_S_(data%id_grz_p) + gp*secs_per_day

         ENDDO

        ! Now manage excretion contributions to DOM/DIM pool
        IF (data%simDCexcr) THEN
           _FLUX_VAR_(data%id_Cexctarget) = &
                           _FLUX_VAR_(data%id_Cexctarget) + excretion*biv !+ doc_excr
           IF(diag_level>2) &
              _DIAG_VAR_(data%id_3d_grz) = _DIAG_VAR_(data%id_3d_grz) + excretion*biv*secs_per_day
        ENDIF
        IF (data%simDNexcr) THEN
           _FLUX_VAR_(data%id_Nexctarget) = _FLUX_VAR_(data%id_Nexctarget) + don_excr
           IF(diag_level>2) &
              _DIAG_VAR_(data%id_3d_exn) = _DIAG_VAR_(data%id_3d_exn)   + don_excr*secs_per_day
           IF(diag_level>9) &
              _DIAG_VAR_S_(data%id_excr_n)=_DIAG_VAR_S_(data%id_excr_n) - don_excr*secs_per_day
        ENDIF
        IF (data%simDPexcr) THEN
           _FLUX_VAR_(data%id_Pexctarget) = _FLUX_VAR_(data%id_Pexctarget) + dop_excr
           IF(diag_level>2) &
              _DIAG_VAR_(data%id_3d_exp) = _DIAG_VAR_(data%id_3d_exp)   + dop_excr*secs_per_day
           IF(diag_level>9) &
              _DIAG_VAR_S_(data%id_excr_p)=_DIAG_VAR_S_(data%id_excr_p) - dop_excr*secs_per_day
        ENDIF

        ! Now manage psuedofaeces, egestion and mortality contributions to POM
        IF (data%simPCexcr) THEN
           _FLUX_VAR_(data%id_Cmorttarget) = _FLUX_VAR_(data%id_Cmorttarget) + poc_excr
           IF(diag_level>2) _DIAG_VAR_(data%id_3d_grz) = _DIAG_VAR_(data%id_3d_grz) + poc_excr*secs_per_day
        ENDIF
        IF (data%simPNexcr) THEN
           _FLUX_VAR_(data%id_Nmorttarget) = _FLUX_VAR_(data%id_Nmorttarget) + pon_excr
           IF(diag_level>2) _DIAG_VAR_(data%id_3d_exn) = _DIAG_VAR_(data%id_3d_exn) + pon_excr*secs_per_day
           IF(diag_level>9) _DIAG_VAR_S_(data%id_egst_n)=_DIAG_VAR_S_(data%id_egst_n) - pon_excr*secs_per_day
        ENDIF
        IF (data%simPPexcr) THEN
           _FLUX_VAR_(data%id_Pmorttarget) = _FLUX_VAR_(data%id_Pmorttarget) + pop_excr
           IF(diag_level>2) _DIAG_VAR_(data%id_3d_exp) = _DIAG_VAR_(data%id_3d_exp) + pop_excr*secs_per_day
           IF(diag_level>9) _DIAG_VAR_S_(data%id_egst_p)=_DIAG_VAR_S_(data%id_egst_p) - pop_excr*secs_per_day
      ENDIF

        ! Now effects of respiration
        !IF (data%id_DOupttarget) THEN
           _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) - respiration*biv
        !ENDIF
        !IF (data%id_CO2upttarget) THEN
        !   _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) + respiration*biv
        !ENDIF

      ENDIF

      ! Export diagnostic variables
      IF (diag_level>0)_DIAG_VAR_S_(data%id_tbiv) = _DIAG_VAR_S_(data%id_tbiv) + biv
      IF (diag_level>1)_DIAG_VAR_S_(data%id_nmp)  = (ingestion-respiration-excretion-egestion-mortality)*biv*secs_per_day
      IF (diag_level>1)_DIAG_VAR_S_(data%id_resp) = respiration*secs_per_day
      IF (diag_level>9) THEN
        _DIAG_VAR_S_(data%id_grz)  = grazing*secs_per_day
        _DIAG_VAR_S_(data%id_mort) = mortality*secs_per_day
        _DIAG_VAR_S_(data%id_excr) = excretion*secs_per_day
        _DIAG_VAR_S_(data%id_egst) = egestion*secs_per_day
        _DIAG_VAR_S_(data%id_fT) = f_Temp
        _DIAG_VAR_S_(data%id_fD) = f_Dens
        _DIAG_VAR_S_(data%id_fG) = fGrazing_Limitation
        _DIAG_VAR_S_(data%id_FR) = FR
        _DIAG_VAR_S_(data%id_pf) = psuedofaeces*secs_per_day
      ENDIF

      IF (data%initFromDensity .AND. diag_level>0) THEN
        _DIAG_VAR_S_(data%id_x_c) = _DIAG_VAR_S_(data%id_tbiv)/_DIAG_VAR_S_(data%id_bnum)
      ENDIF

      ! Update biv_tracer
      IF (data%id_bivtr>0) THEN
        _FLUX_VAR_(data%id_bivtr) = _FLUX_VAR_(data%id_bivtr) + (fr*biv)*(1.-MIN(bt,1.)) - Rbt
      END IF
   ENDDO
!
END SUBROUTINE aed_calculate_benthic_bivalve
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_bivalve_respiration(data,biv_i,iteg,temp,sal) RESULT(resp)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER  :: biv_i ! Invertebrate group
   AED_REAL, INTENT(IN)                   :: temp  ! Temp value being used
   AED_REAL, INTENT(IN)                   :: sal   ! Salinity value being used
   AED_REAL, INTENT(IN)                   :: iteg  ! Ingestion-Egestion
!
!LOCALS
   AED_REAL :: W, TmaxR, maxTR, VV,WW,YY,XX,fT
   AED_REAL :: resp, Q, R20
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Get R20 value
   IF(data%bivalves(biv_i)%fR20==1) THEN
     ! Compute respiration rate coefficient from length
     W = (0.071/1000.) * data%bivalves(biv_i)%Length**2.8
     R20 = data%bivalves(biv_i)%War * W** data%bivalves(biv_i)%Wbr
   ELSE
     ! Use constant respiration rate coefficient
     R20 = data%bivalves(biv_i)%Rresp
   END IF

   ! Now compute actual respiration
   IF(data%bivalves(biv_i)%fR==1) THEN
     ! Option 1) Spillman et al 2008; Bocioniv et al 2013
     resp = R20 * data%bivalves(biv_i)%theta_resp**(temp-20.0)

   ELSEIF(data%bivalves(biv_i)%fR==2) THEN
     ! Option 2) Modified Schneider 1992; Bierman 2005; Gudimov et al. 2015

     TmaxR = data%bivalves(biv_i)%TmaxR
     maxTR = data%bivalves(biv_i)%maxTR
     Q     = data%bivalves(biv_i)%Qresp

     VV = ((TmaxR - temp)/(TmaxR - maxTR))
     WW = LOG(Q)*(TmaxR - maxTR)
     YY = LOG(Q)*(TmaxR - maxTR + 2.)
     XX = (WW * (1. + SQRT(1. + (40. / YY))) / 20.)**2
     fT = VV**XX * EXP(XX*(1.-VV))

     resp = R20 * fT + data%bivalves(biv_i)%SDA * iteg

   ELSE
     ! Unknown option
     resp = data%bivalves(biv_i)%Rresp

   END IF

!   ! Add the salinity stress effect.
!   resp = resp * fSalinity_Limitation(data%bivalves,biv_i,sal)

END FUNCTION aed_bivalve_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION f_DO(data,biv,oxy)
!-------------------------------------------------------------------------------
! Dissolved oxygen effect on bivalve mortality
!
! Add f(DO) effect to base mortality, ie. mort is only ever enhanced by low DO
! Two options:
! 1) Abrupt increase in mortality once below a minimum DO threshold is exceeded
! 2) Steady exponential increase as DO decreases, enhancing mortality
!    by K_BDO so that when DO = 0, f(DO) = 1 + K_BDO
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data
   INTEGER,INTENT(in)   :: biv              ! Invertebrate group
   AED_REAL, INTENT(IN) :: oxy              ! Dissolved oxygen
!
!LOCALS
   AED_REAL :: bot
!
!-------------------------------------------------------------------------------
!BEGIN

   IF(data%bivalves(biv)%fDO==0) THEN
     ! Option 0) Abrupt increase below DO threshold
      IF (oxy < data%bivalves(biv)%KDO) THEN
        f_DO = 10.0
      ELSE
        f_DO = 1.0
      END IF

   ELSEIF(data%bivalves(biv)%fDO==1) THEN
     ! Option 1) Exponential increase: steepness determined by KDO and K_BDO
     ! i.e. when DO = 0, final f(DO) = 1 + K_BDO
     bot  = data%bivalves(biv)%KDO + oxy
     f_DO = 1.0 + data%bivalves(biv)%K_BDO * data%bivalves(biv)%KDO/bot

   ELSE
     ! Unknown option
     f_DO = 1.0
   END IF

END FUNCTION f_DO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fTemp_function_biv(data,biv,temp)
!-------------------------------------------------------------------------------
! Temperature growth multiplier for bivalves
!
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                    :: biv   ! Invertebrate group
   AED_REAL, INTENT(IN)                   :: temp  ! Temp value being used
!
!LOCALS
   AED_REAL , PARAMETER :: a = 1.00
   AED_REAL             :: MINt,Tmin,Tmax,MAXt

!
!-------------------------------------------------------------------------------
!BEGIN

   MINt = data%bivalves(biv)%MINt
   Tmin = data%bivalves(biv)%Tmin
   Tmax = data%bivalves(biv)%Tmax
   MAXt = data%bivalves(biv)%MAXt

   ! If temp below extreme, temp fn = 0, ie filtration ceases
   IF(temp <= MINt) &
      fTemp_function_biv = zero_

   ! If below min temp for optimal range, limited as below
   IF (temp>MINt .and. temp<Tmin) THEN
      fTemp_function_biv = a*(1.0/(-((Tmin-MINt)*(Tmin-MINt)/ &
            (Tmin*Tmin))+(2.0*(Tmin-MINt)/ &
            Tmin)))*(-((temp-MINt)*(temp-MINt)/&
            (Tmin*Tmin))+2.0*((temp-MINt) / Tmin))
   END IF

   ! If between Tmin and Tmax, ie in the opt temp range, then not limited
   IF (temp >=Tmin .and. temp<=Tmax) THEN
      fTemp_function_biv = a
   END IF

   ! If above max temp for optimal range, limited as below
   IF (temp>Tmax .and. temp <MAXt) THEN
      fTemp_function_biv = a*((-(temp*temp)+2.0*(Tmax*temp)- &
            2.0*(Tmax*MAXt)+(MAXt*MAXt))/ &
            ((Tmax-MAXt)*(Tmax-MAXt)))
   END IF

   ! If temp above extreme, temp fn = 0, ie filtration ceases
   IF (temp>=MAXt) &
       fTemp_function_biv = zero_

END FUNCTION fTemp_function_biv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fSS_function(data,biv,SS)
!-------------------------------------------------------------------------------
! Suspended solids function for clams
!
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                    :: biv   ! Invertebrate group
   AED_REAL, INTENT(IN)                   :: SS    ! Susp solids value being used
!
!LOCALS
   AED_REAL , PARAMETER :: a = 1.00
   AED_REAL             :: pseudo, maxSS

!
!-------------------------------------------------------------------------------
!BEGIN

   fSS_function = one_

   pseudo =  data%bivalves(biv)%SSmax    ! SS conc where ingestion decreases
   maxSS  =  data%bivalves(biv)%maxSS    ! SS conc where mussel is buggered


   IF(SS <= pseudo) THEN
        fSS_function = one_
   ELSEIF(SS > pseudo .AND. SS <= maxSS) THEN
        ! From Spillman et al 2008
        fSS_function = (-SS*SS + 2.0*pseudo*SS - &
               2.0*pseudo*maxSS + maxSS*maxSS)/  &
               ((pseudo-maxSS) * (pseudo-maxSS))
   ELSE
        fSS_function = zero_
   END IF

   ! Ensure suspended solids function is not negative
   IF (fSS_function <= zero_) &
        fSS_function = zero_

END FUNCTION fSS_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fD_function(data,biv,dens)
!-------------------------------------------------------------------------------
! Density dependance for bivalves
!
!ARGUMENTS
   CLASS (aed_bivalve_data_t),INTENT(in) :: data  ! Module data, with params
   INTEGER, INTENT(in)                    :: biv   ! Invertebrate group
   AED_REAL, INTENT(in)                   :: dens  ! Density value being used
!
!LOCALS
   AED_REAL             :: maxD, Dmax

!
!-------------------------------------------------------------------------------
!BEGIN

   fD_function = one_

   Dmax =  data%bivalves(biv)%Dmax    ! density where ingestion decreases
   maxD =  data%bivalves(biv)%maxD    ! density where mussel is buggered


   IF(dens <= Dmax) THEN
        fD_function = one_
   ELSEIF(dens > Dmax .AND. dens <= maxD) THEN
        fD_function = ( -dens*dens + 2.*Dmax*dens - 2.*Dmax*maxD + maxD*maxD ) &
                    / ( (Dmax-maxD)*(Dmax-maxD) )
   ELSE
        fD_function = zero_
   END IF

   ! Ensure fD is not negative
   IF (fD_function <= zero_) fD_function = zero_

END FUNCTION fD_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if 0
!###############################################################################
FUNCTION clamSalinity(Sbot) RESULT(ret)
!-------------------------------------------------------------------------------
! Salinity tolerance of clams/mussels
   !-- Incoming
   REAL (r_wq) :: Sbot(:)            ! Salinity in the bottom layer
   !-- Returns the salinity function
   REAL (r_wq), DIMENSION(SIZE(Sbot)) :: ret
   REAL (r_wq), DIMENSION(SIZE(Sbot)) :: tmp1,tmp2

   ! BSmin minimum sal and BSmax is max sal set in WQcons
   ! Don't seem to need Bep, Aep or Sop from consts file
   ! Salinity is non-limiting for sals > BSmin and <BSmax
   ! clamLowCut is salinity when shells close up and feeding ceases
   ! CS 020609 Need to make an exception if FW species and BSmin and
   ! clamLowCut are both equal to zero.

   ! Salinity is within the tolerance; no limitation.
   ! If sal in bott cell is > min sal & < max sal set in WQcons, pf=1
   WHERE(Sbot >= BSmin .and. Sbot <= BSmax)
     ret = wq_one
   END WHERE

   ! Salinity is greater than the upper bound
   ! maxS is set in globals at 45psu
   WHERE(Sbot > BSmax)
     ret = (-Sbot*Sbot+2.0*BSmax*Sbot-  &
            2.0*BSmax*maxS+maxS*maxS)/((BSmax-maxS)*(BSmax-maxS))
   END WHERE

   ! Salinity is less than the lower bound but greater than low cut
   ! If sal is < min set in WQcons but > clamLowCut set at clamCons.dat)
   tmp1 = wq_zero
   tmp2 = wq_zero
   WHERE(Sbot < BSmin .AND. Sbot > clamLowCut)
     tmp1 = Sbot-clamLowCut
     tmp2 = BSmin-clamLowCut
     ret =  (2*tmp1/BSmin-(tmp1*tmp1/(BSmin*BSmin)))/ &
            (2*tmp2/BSmin-(tmp2*tmp2/(BSmin*BSmin)))
   END WHERE

   ! Salinity is less than the clamLowCut
   ! If sal < lowest sal (hardwired at start of fn), shells close
   WHERE(Sbot <= clamLowCut)
     ret = wq_zero
   END WHERE

   ! If lower bound and low cut are both zero i.e. Freshwater species
   ! then need to set f(S) to one
   IF(BSmin==wq_zero) THEN
     WHERE(Sbot <= BSmin)
         ret =  wq_one
     END WHERE
   ENDIF

   ! Ensure temp function is not negative
   WHERE(ret <= wq_zero)
     ret = wq_zero
   END WHERE

 END FUNCTION clamSalinity
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


END MODULE aed_bivalve
