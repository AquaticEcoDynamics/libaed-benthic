!###############################################################################
!#                                                                             #
!# aed_macrophyte.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2015 - 2022 -  The University of Western Australia               #
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
!# Created July 2015                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_macrophyte
!-------------------------------------------------------------------------------
!  aed_macrophyte --- multi-group macrophyte (/seagrass) model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util
   USE aed_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_macrophyte_data_t

!  %% NAMELIST    %% type :  macrophyte_params_t
   TYPE :: macrophyte_params_t
       INTEGER       :: growthForm
       CHARACTER(64) :: m_name
       AED_REAL      :: m0
       AED_REAL      :: R_growth
       INTEGER       :: fT_Method
       AED_REAL      :: theta_growth
       AED_REAL      :: T_std
       AED_REAL      :: T_opt
       AED_REAL      :: T_max
       INTEGER       :: lightModel
       AED_REAL      :: I_K
       AED_REAL      :: I_S
       AED_REAL      :: KeMAC
       AED_REAL      :: f_pr
       AED_REAL      :: R_resp
       AED_REAL      :: theta_resp
       INTEGER       :: salTol
       AED_REAL      :: S_bep
       AED_REAL      :: S_maxsp
       AED_REAL      :: S_opt
       AED_REAL      :: K_CD
       AED_REAL      :: f_bg
       AED_REAL      :: k_omega
       AED_REAL      :: Xcc
       AED_REAL      :: K_N
       AED_REAL      :: X_ncon
       AED_REAL      :: K_P
       AED_REAL      :: X_pcon
   END TYPE
!  %% END NAMELIST    %% type :  macrophyte_params_t

   TYPE,extends(aed_model_data_t) :: aed_macrophyte_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_mphy(:)
      INTEGER :: id_par, id_tem, id_sal, id_dz, id_extc, id_I_0
      INTEGER :: id_d_par, id_gpp, id_p2r, id_mac, id_sed_zone
      INTEGER :: id_mac_ag, id_mac_bg, id_lai, id_root_d, id_root_o
      INTEGER :: id_atem,id_theta
      AED_REAL,ALLOCATABLE :: active_zones(:)

      !# Model parameters
      INTEGER  :: num_mphy
      INTEGER  :: n_zones
      LOGICAL  :: simMacFeedback, simStaticBiomass
      TYPE(macrophyte_params_t),DIMENSION(:),ALLOCATABLE :: mphydata
      AED_REAL :: bg_gpp_frac,coef_bm_hgt

     CONTAINS
         PROCEDURE :: define             => aed_define_macrophyte
         PROCEDURE :: initialize_benthic => aed_initialize_benthic_macrophyte
         PROCEDURE :: calculate_benthic  => aed_calculate_benthic_macrophyte
        !PROCEDURE :: calculate_riparian => aed_calculate_riparian_macrophyte
         PROCEDURE :: bio_drag           => aed_bio_drag_macrophyte
         PROCEDURE :: light_extinction   => aed_light_extinction_macrophyte
   END TYPE


   INTEGER, PARAMETER :: SUBMERGED = 1
   INTEGER, PARAMETER :: EMERGENT  = 2
   INTEGER, PARAMETER :: FLOATING  = 3
   INTEGER            :: diag_level = 10     ! 0 = no diagnostic outputs
                                             ! 1 = basic diagnostic outputs
                                             ! 2 = flux rates, and supporitng
                                             ! 3 = other metrics
                                             !10 = all debug & checking outputs

CONTAINS
!===============================================================================




!###############################################################################
INTEGER FUNCTION load_csv(dbase, md)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(macrophyte_params_t) :: md(MAX_PHYTO_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0, ret = 0
   LOGICAL :: meh
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
         md(ccol)%m_name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('m0')           ; md(ccol)%m0           = extract_double(values(ccol))
            CASE ('R_growth')     ; md(ccol)%R_growth     = extract_double(values(ccol))
            CASE ('fT_Method')    ; md(ccol)%fT_Method    = extract_integer(values(ccol))
            CASE ('theta_growth') ; md(ccol)%theta_growth = extract_double(values(ccol))
            CASE ('T_std')        ; md(ccol)%T_std        = extract_double(values(ccol))
            CASE ('T_opt')        ; md(ccol)%T_opt        = extract_double(values(ccol))
            CASE ('T_max')        ; md(ccol)%T_max        = extract_double(values(ccol))
            CASE ('lightModel')   ; md(ccol)%lightModel   = extract_integer(values(ccol))
            CASE ('I_K')          ; md(ccol)%I_K          = extract_double(values(ccol))
            CASE ('I_S')          ; md(ccol)%I_S          = extract_double(values(ccol))
            CASE ('KeMAC')        ; md(ccol)%KeMAC        = extract_double(values(ccol))
            CASE ('f_pr')         ; md(ccol)%f_pr         = extract_double(values(ccol))
            CASE ('R_resp')       ; md(ccol)%R_resp       = extract_double(values(ccol))
            CASE ('theta_resp')   ; md(ccol)%theta_resp   = extract_double(values(ccol))
            CASE ('salTol')       ; md(ccol)%salTol       = extract_integer(values(ccol))
            CASE ('S_bep')        ; md(ccol)%S_bep        = extract_double(values(ccol))
            CASE ('S_maxsp')      ; md(ccol)%S_maxsp      = extract_double(values(ccol))
            CASE ('S_opt')        ; md(ccol)%S_opt        = extract_double(values(ccol))
            CASE ('K_CD')         ; md(ccol)%K_CD         = extract_double(values(ccol))
            CASE ('f_bg')         ; md(ccol)%f_bg         = extract_double(values(ccol))
            CASE ('k_omega')      ; md(ccol)%k_omega      = extract_double(values(ccol))
            CASE ('Xcc')          ; md(ccol)%Xcc          = extract_double(values(ccol))
            CASE ('K_N')          ; md(ccol)%K_N          = extract_double(values(ccol))
            CASE ('X_ncon')       ; md(ccol)%X_ncon       = extract_double(values(ccol))
            CASE ('K_P')          ; md(ccol)%K_P          = extract_double(values(ccol))
            CASE ('X_pcon')       ; md(ccol)%X_pcon       = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit) !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_macrophyte_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(macrophyte_params_t),ALLOCATABLE :: md(:)
   NAMELIST /macrophyte_data/ md  ! %% type : macrophyte_params_t - see above
!-------------------------------------------------------------------------------
!BEGIN
    ALLOCATE(md(MAX_PHYTO_TYPES))
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, md)
       CASE (NML_TYPE)
           tfil = find_free_lun()
           open(tfil,file=dbase, status='OLD', iostat=status)
           IF (status /= 0) STOP 'Cannot open macrophyte_data namelist file for macrophytes'
           read(tfil,nml=macrophyte_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist macrophyte_data for macrophytes'

    data%num_mphy = count
    ALLOCATE(data%mphydata(count))

    ALLOCATE(data%id_mphy(count))

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%mphydata(i)%growthForm   = 1
         IF(list(i)>6) data%mphydata(i)%growthForm   = 2   ! HACK FOR GELERAH
       data%mphydata(i)%m_name       = md(list(i))%m_name
       data%mphydata(i)%m0           = md(list(i))%m0
       data%mphydata(i)%R_growth     = md(list(i))%R_growth/secs_per_day
       data%mphydata(i)%fT_Method    = md(list(i))%fT_Method
       data%mphydata(i)%theta_growth = md(list(i))%theta_growth
       data%mphydata(i)%T_std        = md(list(i))%T_std
       data%mphydata(i)%T_opt        = md(list(i))%T_opt
       data%mphydata(i)%T_max        = md(list(i))%T_max
       data%mphydata(i)%lightModel   = md(list(i))%lightModel
       data%mphydata(i)%I_K          = md(list(i))%I_K
       data%mphydata(i)%I_S          = md(list(i))%I_S
       data%mphydata(i)%KeMAC        = md(list(i))%KeMAC
       data%mphydata(i)%f_pr         = md(list(i))%f_pr
       data%mphydata(i)%R_resp       = md(list(i))%R_resp/secs_per_day
       data%mphydata(i)%theta_resp   = md(list(i))%theta_resp
       data%mphydata(i)%salTol       = md(list(i))%salTol
       data%mphydata(i)%S_bep        = md(list(i))%S_bep
       data%mphydata(i)%S_maxsp      = md(list(i))%S_maxsp
       data%mphydata(i)%S_opt        = md(list(i))%S_opt
       data%mphydata(i)%K_CD         = md(list(i))%K_CD
       data%mphydata(i)%f_bg         = md(list(i))%f_bg
       data%mphydata(i)%k_omega      = md(list(i))%k_omega
       data%mphydata(i)%Xcc          = md(list(i))%Xcc
       data%mphydata(i)%K_N          = md(list(i))%K_N
       data%mphydata(i)%X_ncon       = md(list(i))%X_ncon
       data%mphydata(i)%K_P          = md(list(i))%K_P
       data%mphydata(i)%X_pcon       = md(list(i))%X_pcon

       ! Register group as a state variable
       data%id_mphy(i) = aed_define_sheet_variable(                    &
                              md(list(i))%m_name,                      &
                              'mmol C/m2', 'macrophyte biomass',       &
                              md(list(i))%m0,                          &
                              minimum=zero_)
    ENDDO
    DEALLOCATE(md)
END SUBROUTINE aed_macrophyte_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_macrophyte(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the macrophyte/seagrass  model
!
!  Here, the aed_ namelist is read and the variables exported
!  by the model are registered
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status, i

!  %% NAMELIST   %%  /aed_macrophyte/
!  %% Last Checked 20/08/2021
   INTEGER            :: num_mphy = 0
   INTEGER            :: the_mphy(MAX_PHYTO_TYPES) = 0
   CHARACTER(len=128) :: dbase = 'aed_macrophyte_pars.nml'
   INTEGER            :: n_zones = 0
   INTEGER            :: active_zones(MAX_ZONES)
   LOGICAL            :: simMacFeedback = .FALSE.
   LOGICAL            :: simStaticBiomass = .FALSE.
   AED_REAL           :: bg_gpp_frac = 0.1
   AED_REAL           :: coef_bm_hgt = 1e-5
!  %% END NAMELIST   %%  /aed_macrophyte/

   NAMELIST /aed_macrophyte/ num_mphy, the_mphy, dbase, n_zones, active_zones, &
                              simMacFeedback, simStaticBiomass,                &
                              bg_gpp_frac,coef_bm_hgt

!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_macrophyte configuration"

! now done in the declaration
!  simMacFeedback = .FALSE.
!  simStaticBiomass = .FALSE.

   ! Read the namelist
   read(namlst,nml=aed_macrophyte,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_macrophyte'

   data%simMacFeedback = simMacFeedback
   data%simStaticBiomass = simStaticBiomass
   data%bg_gpp_frac = bg_gpp_frac  ! make this species specific
   data%coef_bm_hgt = coef_bm_hgt  ! make this species specific

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   ! Macrophyte state variable allocated in here
   CALL aed_macrophyte_load_params(data, dbase, num_mphy, the_mphy)

  !CALL aed_bio_temp_function(data%num_phytos,                &
  !                             data%phytos%theta_growth,     &
  !                             data%phytos%T_std,            &
  !                             data%phytos%T_opt,            &
  !                             data%phytos%T_max,            &
  !                             data%phytos%aTn,              &
  !                             data%phytos%bTn,              &
  !                             data%phytos%kTn,              &
  !                             data%phytos%p_name)


   ! Register diagnostic variables
   data%id_d_par  = aed_define_sheet_diag_variable('par','W/m2','benthic light intensity')
   data%id_gpp    = aed_define_sheet_diag_variable('gpp','mmol C/m2/d','benthic plant productivity')
   data%id_p2r    = aed_define_sheet_diag_variable('npp','/d','macrophyte net productivity')
   data%id_mac    = aed_define_sheet_diag_variable('mac_ben','mmol C/m2','total macrophyte biomass')
   data%id_lai    = aed_define_sheet_diag_variable('mac_lai','m2/m2','macrophyte leaf area density')
   data%id_mac_ag = aed_define_sheet_diag_variable('mac_ag','mmol C/m2','total above ground macrophyte biomass')
   data%id_mac_bg = aed_define_sheet_diag_variable('mac_bg','mmol C/m2','total below ground macrophyte biomass')
   data%id_root_d = aed_define_sheet_diag_variable('mac_root_depth','m','mean depth of roots below the sediment surface')
   data%id_root_o = aed_define_sheet_diag_variable('mac_root_o2','mmol O2/m2','mean O2 injection rate of roots into sediment')

   ! Register environmental dependencies
   data%id_tem  = aed_locate_global('temperature')
   data%id_extc = aed_locate_global('extc_coef')
   data%id_sal  = aed_locate_global('salinity')
   data%id_dz   = aed_locate_global('layer_ht')
   data%id_par  = aed_locate_global('par')
   data%id_sed_zone = aed_locate_sheet_global('sed_zone')
   data%id_atem     = aed_locate_sheet_global('air_temp')
   data%id_I_0      = aed_locate_sheet_global('par_sf')

END SUBROUTINE aed_define_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if 0
!# Already in aed_util
!###############################################################################
LOGICAL FUNCTION in_zone_set(matz, active_zones)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: matz
   AED_REAL,INTENT(in) :: active_zones(:)
!
!LOCALS
   INTEGER :: i, l
   LOGICAL :: res
!BEGIN
!-------------------------------------------------------------------------------
   res = .FALSE.
   l = size(active_zones)
   DO i=1,l
      IF ( active_zones(i) == matz ) THEN
         res = .TRUE.
         EXIT
      ENDIF
   ENDDO

   in_zone_set = res
END FUNCTION in_zone_set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif


!###############################################################################
SUBROUTINE aed_initialize_benthic_macrophyte(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to initialize bottom diagnostics, and other checks
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mphy_i
   AED_REAL :: matz, mphy
!-------------------------------------------------------------------------------
!BEGIN

   ! Check to ensure this cell/zone is colonisable
   matz = _STATE_VAR_S_(data%id_sed_zone)
   IF ( .NOT. in_zone_set(matz, data%active_zones) ) THEN
     _DIAG_VAR_S_(data%id_mac) = zero_
     _DIAG_VAR_S_(data%id_mac_ag) = zero_
     _DIAG_VAR_S_(data%id_mac_bg) = zero_
     _DIAG_VAR_S_(data%id_lai) = zero_
     _DIAG_VAR_S_(data%id_gpp) = zero_
     _DIAG_VAR_S_(data%id_root_o) = zero_
     _DIAG_VAR_S_(data%id_root_d) = 0.01
     RETURN
   ENDIF

   ! Set initial diagnostics
   _DIAG_VAR_S_(data%id_mac) = zero_
   DO mphy_i=1,data%num_mphy
       mphy = _STATE_VAR_S_(data%id_mphy(mphy_i))
       _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + mphy
       _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + mphy*(one_-data%mphydata(mphy_i)%f_bg)
       _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + mphy*(data%mphydata(mphy_i)%f_bg)
       _DIAG_VAR_S_(data%id_lai) = _DIAG_VAR_S_(data%id_lai) +       &
                     (one_ - exp(-data%mphydata(mphy_i)%k_omega * mphy*(one_-data%mphydata(mphy_i)%f_bg)))
   ENDDO
   _DIAG_VAR_S_(data%id_gpp) = zero_
   _DIAG_VAR_S_(data%id_root_o) = zero_
   _DIAG_VAR_S_(data%id_root_d) = MAX( MIN(_DIAG_VAR_S_(data%id_mac_ag) * data%coef_bm_hgt,0.25),0.01)

END SUBROUTINE aed_initialize_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_macrophyte(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate submerged macrophyte production and respiration etc
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: mphy        ! State
   INTEGER  :: mphy_i
   AED_REAL :: mphy_flux
   AED_REAL :: fT, fI, fSal, fDO
   AED_REAL :: extc, dz, par, Io, temp, salinity
   AED_REAL :: primprod(data%num_mphy)
   AED_REAL :: respiration(data%num_mphy)
   AED_REAL :: matz
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Retrieve current environmental conditions
   temp = _STATE_VAR_(data%id_tem)      ! local temperature
   salinity = _STATE_VAR_(data%id_sal)  ! local salinity
   par = _STATE_VAR_(data%id_par)       ! local photosynthetically active radiation
   Io = _STATE_VAR_S_(data%id_I_0)      ! surface short wave radiation

   ! Initialise cumulative biomass diagnostics
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_
   _DIAG_VAR_S_(data%id_lai) = zero_
   _DIAG_VAR_S_(data%id_gpp) = zero_

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i))! macrophyte group i

      ! LIGHT
      extc = _STATE_VAR_(data%id_extc)
      dz   = _STATE_VAR_(data%id_dz)     ! dz = 0.5
      fI   = photosynthesis_irradiance(data%mphydata(mphy_i)%lightModel, &
                     data%mphydata(mphy_i)%I_K, data%mphydata(mphy_i)%I_S, par, extc, Io, dz)
      fT   = 1. ! / fTemp_function(temp,     ) !1.
      fSal = fSal_function(salinity,10.,72.,123.,230.)

      ! Primary productivity
      primprod(mphy_i) = data%mphydata(mphy_i)%R_growth * fI * fT * fSal

      ! Respiration and general metabolic loss
      respiration(mphy_i) = bio_respiration(data%mphydata(mphy_i)%R_resp, data%mphydata(mphy_i)%theta_resp, temp)

      ! Salinity stress effect on respiration
      fSal = 1.0 ! phyto_salinity(data%mphydata,mphy_i,salinity)
      fDO  = 1.0 ! phyto_oxygen(data%mphydata,mphy_i,oxygen)

      respiration(mphy_i) = respiration(mphy_i) * fSal * fDO

      IF( .NOT.data%simStaticBiomass ) THEN
        mphy_flux = (primprod(mphy_i) - respiration(mphy_i)) *  mphy
        ! Set bottom fluxes for the pelagic (change per surface area per second)
        _FLUX_VAR_B_(data%id_mphy(mphy_i)) = _FLUX_VAR_B_(data%id_mphy(mphy_i)) + mphy_flux
      ENDIF
      IF( data%simMacFeedback ) THEN
    !    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + mphy_flux
    !    _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - mphy_flux
      ENDIF

     _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + mphy
     _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + mphy*(one_-data%mphydata(mphy_i)%f_bg)
     _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + mphy*(data%mphydata(mphy_i)%f_bg)
     _DIAG_VAR_S_(data%id_lai) = _DIAG_VAR_S_(data%id_lai) +       &
                   (one_ - exp(-data%mphydata(mphy_i)%k_omega * mphy*(one_-data%mphydata(mphy_i)%f_bg)))

     _DIAG_VAR_S_(data%id_gpp) = _DIAG_VAR_S_(data%id_gpp) + primprod(mphy_i)*mphy*secs_per_day

   ENDDO

   _DIAG_VAR_S_(data%id_root_o) = _DIAG_VAR_S_(data%id_gpp) * data%bg_gpp_frac    !mmolO2/m2/d
   _DIAG_VAR_S_(data%id_root_d) = MAX( MIN(_DIAG_VAR_S_(data%id_mac_ag) * data%coef_bm_hgt,0.25),0.01) !m

   ! Export diagnostic variables
   _DIAG_VAR_S_(data%id_d_par)= par
   _DIAG_VAR_S_(data%id_p2r)  = SUM(primprod) - SUM(respiration)
!   IF( SUM(respiration(:)) > 1e-5 ) THEN
!     _DIAG_VAR_S_(data%id_p2r)  = (SUM(primprod(:))/data%num_mphy) / (SUM(respiration(:))/data%num_mphy)
!   ELSE
!     _DIAG_VAR_S_(data%id_p2r)  = 9999.
!   ENDIF


END SUBROUTINE aed_calculate_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_macrophyte(data, column, layer_idx, extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to macrophyte biomass
!
!  WARNING - THIS IS ADDDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: dz,mphy,matz
   INTEGER  :: mphy_i
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i)) *       &
                   (one_-mphy*data%mphydata(mphy_i)%f_bg) ! above ground density of macrophyte group i

      ! Self-shading depending on amount of carbon in water volume
      extinction = extinction + (data%mphydata(mphy_i)%KeMAC * (mphy/dz) )
   ENDDO
END SUBROUTINE aed_light_extinction_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_riparian_macrophyte(data,column,layer_idx,pc_wet)
!-------------------------------------------------------------------------------
! Calculate emergent macrophytes in the riparian zone.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: pc_wet
!
!LOCALS
   AED_REAL :: mphy        ! State
   INTEGER  :: mphy_i
   AED_REAL :: mphy_flux
   AED_REAL :: fT, fI, fSal, fDO, fSM
   AED_REAL :: extc, dz, par, Io, temp, salinity, theta
   AED_REAL :: primprod(data%num_mphy)
   AED_REAL :: respiration(data%num_mphy)
   AED_REAL :: matz
   AED_REAL :: emergent_portion
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Retrieve current environmental conditions
   temp = _STATE_VAR_S_(data%id_atem)      ! local air temperature
   theta = 1.0 !_STATE_VAR_(data%id_theta)    ! local soil moisture
   Io = _STATE_VAR_S_(data%id_I_0)       ! surface short wave radiation

   ! Initialise cumulative biomass diagnostics
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_
   _DIAG_VAR_S_(data%id_lai) = zero_

   DO mphy_i=1,data%num_mphy

     ! Check for emergent species ( can grow in dry and wet )
     IF( data%mphydata(mphy_i)%GrowthForm == EMERGENT ) THEN

        ! Retrieve current (local) state variable values
        mphy = _STATE_VAR_S_(data%id_mphy(mphy_i))! macrophyte group i

        ! LIGHT
        IF( pc_wet >0.8 ) THEN
          par   = _STATE_VAR_(data%id_par)
          emergent_portion = 1.0 ! will depend on pc_wet and depth of water
        ELSE
          par   = 0.45*Io
          emergent_portion = 1.0
        ENDIF
        fI   = photosynthesis_irradiance(data%mphydata(mphy_i)%lightModel, &
                       data%mphydata(mphy_i)%I_K, data%mphydata(mphy_i)%I_S, par, extc, Io, dz)
        fI   = fI * emergent_portion

        ! TEMPERATURE
        fT   = 1.

        ! SOIL WATER AVAILABILITY
        IF( pc_wet >0.8 ) THEN
          fSM   = 1.
        ELSE
          fSM   = 0.1
        ENDIF

        primprod(mphy_i) = data%mphydata(mphy_i)%R_growth * fI * fT * fSM

        ! Respiration and general metabolic loss
        respiration(mphy_i) = bio_respiration(data%mphydata(mphy_i)%R_resp, data%mphydata(mphy_i)%theta_resp, temp)

        ! Salinity stress effect on respiration
        fSal = 1.0 ! phyto_salinity(data%mphydata,mphy_i,salinity)
        fDO  = 1.0 ! phyto_oxygen(data%mphydata,mphy_i,oxygen)

        respiration(mphy_i) = respiration(mphy_i) * fSal * fDO

        IF( .NOT.data%simStaticBiomass ) THEN
          mphy_flux = (primprod(mphy_i) - respiration(mphy_i)) *  mphy

          ! Set bottom fluxes for the pelagic (change per surface area per second)
          _FLUX_VAR_B_(data%id_mphy(mphy_i)) = _FLUX_VAR_B_(data%id_mphy(mphy_i)) + mphy_flux
        ENDIF
        IF( data%simMacFeedback ) THEN
        !    _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + mphy_flux
        !    _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - mphy_flux
        ENDIF

        _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + mphy
        _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + mphy*(one_-data%mphydata(mphy_i)%f_bg)
        _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + mphy*(data%mphydata(mphy_i)%f_bg)
        _DIAG_VAR_S_(data%id_lai) = _DIAG_VAR_S_(data%id_lai) +    &
                         (one_ - exp(-data%mphydata(mphy_i)%k_omega * mphy*(one_-data%mphydata(mphy_i)%f_bg)))

      ELSEIF( data%mphydata(mphy_i)%GrowthForm == SUBMERGED ) THEN
        ! Dry cell
        IF( pc_wet <0.5 ) THEN
          ! Retrieve current (local) state variable values
          mphy = _STATE_VAR_S_(data%id_mphy(mphy_i))! macrophyte group i

          IF( .NOT.data%simStaticBiomass ) THEN
             mphy_flux = -(0.5/86400) *  mphy

            ! Set bottom fluxes for the pelagic (change per surface area per second)
            _FLUX_VAR_B_(data%id_mphy(mphy_i)) = _FLUX_VAR_B_(data%id_mphy(mphy_i)) + mphy_flux
          ENDIF
          _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + mphy
          _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + mphy*(one_-data%mphydata(mphy_i)%f_bg)
          _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + mphy*(data%mphydata(mphy_i)%f_bg)
          _DIAG_VAR_S_(data%id_lai) = _DIAG_VAR_S_(data%id_lai) +    &
                       (one_ - exp(-data%mphydata(mphy_i)%k_omega * mphy*(one_-data%mphydata(mphy_i)%f_bg)))

      ENDIF
    ENDIF
   ENDDO

   ! Export diagnostic variables
   IF( pc_wet <1.0 ) _DIAG_VAR_S_(data%id_d_par)= Io
   _DIAG_VAR_S_(data%id_gpp) = zero_ !_DIAG_VAR_S_(data%id_gpp) + SUM(primprod)*secs_per_day


END SUBROUTINE aed_calculate_riparian_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_drag_macrophyte(data, column, layer_idx, drag)
!-------------------------------------------------------------------------------
! Get the effect of macrophyte biomass on benthic drag
!
!  WARNING - THIS IS ADDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   AED_REAL :: dz,mphy,matz
   INTEGER  :: mphy_i
!-------------------------------------------------------------------------------
!BEGIN
   mphy = zero_  !## CAB [-Wmaybe-uninitialized]

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   DO mphy_i=1,data%num_mphy
      ! Retrieve current (local) state variable values
      dz   = _STATE_VAR_(data%id_dz)  ! dz = 0.5

      !# above ground density of macrophyte group i
      mphy = _STATE_VAR_S_(data%id_mphy(mphy_i)) * (one_-mphy*data%mphydata(mphy_i)%f_bg)

      !# additional drag due to biomass
      drag = drag + (data%mphydata(mphy_i)%K_CD * (mphy/dz) )
   ENDDO
END SUBROUTINE aed_bio_drag_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_macrophyte
