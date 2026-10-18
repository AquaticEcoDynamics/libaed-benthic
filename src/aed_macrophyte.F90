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
!#  Copyright 2015-2026 : The University of Western Australia                  #
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
!# Created July 2024                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_macrophyte
!===============================================================================
!  aed_macrophyte --- multi-group macrophyte (/seagrass) model, with options
!
!  ORDER OF CALLS, each time step (every host calls them in this order)
!
!   1. calculate_column   aed_calculate_column_macrophyte    epi_model 2 only
!        Works up from the bed through the layers the canopy occupies, with the
!        canopy geometry from the current leaf C (mac_canopy). In each layer:
!        canopy share, light, temperature, nutrients -> epiphyte net production,
!        stored as diagnostics (no fluxes are set here):
!          epi_npp  each layer (mmol C/m3/d);  epi_net  column total (mmol C/m2/d)
!          epi_col = 1  (the column result is ready)
!   2. calculate_benthic  aed_calculate_benthic_macrophyte   always
!        a. plants: growth, respiration, mortality, translocation, fruiting, and
!           their exchange with the bottom water layer and the sediment
!        b. canopy geometry diagnostics (height, shoot density, biovolume)
!        c. epiphytes: the state changes by net production less losses (R_epib,
!           and sloughing with the leaves), the losses going to the mortality links
!             epi_model 1, or 2 with epi_col = 0 (host made no column call):
!               bulk calculation, exchanging with the bottom water layer here
!             epi_model 2 with epi_col = 1: takes epi_net, and sets epi_col = 2
!   3. calculate          aed_calculate_macrophyte           epi_model 2 only
!        In each layer, if epi_col = 2: applies that layer's epi_npp to its O2,
!        DIC, N and P (or DON/DOP), matching the state change of step 2c.
!
!  Benthic fluxes (step 2) are areal (mmol/m2/s), pelagic ones (step 3) are per
!  m3 of the layer. No fluxes are set in step 1, as hosts treat them differently.
!  Hosts zero diagnostics at each step, so nothing is carried between steps in
!  them: the canopy geometry is recomputed from the state where needed.
!  With sediment zones (the module flagged for zone averaging, eg GLM
!  benthic_mode 2, 3), step 2 runs on each zone but step 3 on the whole column,
!  so their results can't be matched: epi_model 2 is then set to 1 at startup.
!===============================================================================
   USE aed_core
   USE aed_util
  !USE aed_maths
   USE aed_bio_utils
  !USE kobe_sav_posture

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_macrophyte_data_t

!-------------------------------------------------------------------------------
!  MODULE SPECIES PARAMETER TYPE DEFINITION

!  %% NAMELIST    %% type :  macrophyte_params_t
   TYPE :: macrophyte_params_t

       CHARACTER(64) :: m_name
       INTEGER       :: growth_form
       INTEGER       :: zone_lock
       AED_REAL      :: m_initial
       AED_REAL      :: m_initial_minb
       AED_REAL      :: m_initial_mind
       AED_REAL      :: m_initial_d0 = 1.2 ! mac_initial 1: depth above which biomass is m_initial
       AED_REAL      :: m0

       ! Canopy/physical characteristics
       INTEGER       :: drag_model
       AED_REAL      :: K_CD
       AED_REAL      :: shoot_diameter   ! shoot/leaf diameter or blade width
       AED_REAL      :: shoot_length     ! shoot/leaf typical length
       AED_REAL      :: tissue_density   ! shoot/leaf density (kg/m3 of tissue)
       AED_REAL      :: elastic_modulus  ! shoot/leaf elasticity
       INTEGER       :: shading_model
       AED_REAL      :: KeMAC
       AED_REAL      :: kA               ! shelf shading effect
       AED_REAL      :: k_omega          ! Omega_MAC?
       AED_REAL      :: sine_blade       ! sine blade shape
       AED_REAL      :: X_ldw            ! leaf length to biomass
       AED_REAL      :: K_ldw            ! half-biomass of leaf to biomass curve
       AED_REAL      :: X_sdw            ! shoot to biomass ratio

       ! Growth rate options/parameters
       AED_REAL      :: R_growth
       INTEGER       :: temp_method
       AED_REAL      :: theta_growth
       AED_REAL      :: T_std
       AED_REAL      :: T_opt
       AED_REAL      :: T_max
       AED_REAL      :: kTn, aTn, bTn
       INTEGER       :: light_model         ! Light_model
       AED_REAL      :: I_K
       AED_REAL      :: I_S
       AED_REAL      :: f_pr
       AED_REAL      :: K_epi

       ! Loss rate options/parameters
       AED_REAL      :: R_resp             ! maximum AG respiration rate, /day
       AED_REAL      :: theta_resp
       AED_REAL      :: E_comp             ! compensation scalar PAR irradiance, mol photon/m2/d;
       AED_REAL      :: R_mort             ! maximum AG mortality rate, /day
       INTEGER       :: sal_method
       AED_REAL      :: S_bep
       AED_REAL      :: S_maxsp
       AED_REAL      :: S_opt

       ! Belowground options/parameters
       AED_REAL      :: m_initial_bg
       AED_REAL      :: f_bg               ! f_below
       AED_REAL      :: tau_tran           ! AG/BG translocation rate, /day
       AED_REAL      :: R_resp_bg          ! maximum BG respiration rate, /day
       AED_REAL      :: R_mort_bg          ! maximum BG mortality rate, /day

       ! Stoichiometry and limitation parameters
       AED_REAL      :: X_cchl
       AED_REAL      :: X_cdw
       AED_REAL      :: K_N
       AED_REAL      :: X_ncon
       AED_REAL      :: K_P
       AED_REAL      :: X_pcon
       AED_REAL      :: K_S
       AED_REAL      :: K_NPP

       ! Fruiting/recruitment options/parameters
       INTEGER       :: fruit_model         ! Fruiting
       AED_REAL      :: m_initial_f
       AED_REAL      :: tau_tran_fruit  ! translocation rate from leaves to fruit, /day
       AED_REAL      :: r_release         ! fruit releasing rate, /day (assume all fruit released in 10 days)
       AED_REAL      :: f_seed          ! seed biomass/(AG+seed biomass) fraction
       AED_REAL      :: t_start_g       ! calenda day to grow, days of the year
       AED_REAL      :: t_dur_g          ! duration to reach maximum growth, days
       AED_REAL      :: t_start_r       ! calenda day to release fruits, days of the year
       AED_REAL      :: t_dur_r           ! duration to reach maximum release, days

   END TYPE
!  %% END NAMELIST    %% type :  macrophyte_params_t
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  MODULE DATA STRUCTURE DEFINITION

   TYPE,extends(aed_model_data_t) :: aed_macrophyte_data_t
      !# Variable identifiers (local and linked)
      INTEGER,ALLOCATABLE :: id_mphya(:),id_mphyb(:),id_mphyf(:)
      INTEGER :: id_epi, id_epib, id_epig, id_epir, id_epil
      INTEGER :: id_epi_npp, id_epi_net, id_epi_col   ! epi_model 2, column -> benthic/pelagic
      INTEGER :: id_par, id_tem, id_sal, id_extc, id_I_0, id_yearday, id_vel
      INTEGER :: id_depth, id_dz
      INTEGER :: id_atem,id_theta
      INTEGER :: id_mac_ag, id_mac_bg, id_root_d, id_root_o, id_mac_tr
      INTEGER :: id_d_par, id_gpp, id_p2r, id_mac, id_sed_zone
      INTEGER :: id_shdiam(MAX_PHYTO_TYPES)
      INTEGER :: id_shdens(MAX_PHYTO_TYPES)
      INTEGER :: id_shhght(MAX_PHYTO_TYPES)
      INTEGER :: id_shabve(MAX_PHYTO_TYPES)
      INTEGER :: id_shbvol(MAX_PHYTO_TYPES)
      INTEGER :: id_canopy_stem_diam, id_canopy_stem_dens, id_canopy_height, id_canopy_biovolume
      INTEGER :: id_canopy_sh_dens, id_canopy_sh_diam, id_canopy_frarea
      INTEGER :: id_canopy_velocity, id_canopy_drag, id_canopy_blockage, id_canopy_lai
      INTEGER :: id_parc, id_aeff, id_kemac, id_mac_fr, id_mac_ft, id_mac_gt
      INTEGER,ALLOCATABLE :: id_fT(:), id_fSal(:)
      INTEGER :: id_oxy, id_dic, id_nox, id_nh4, id_po4
      INTEGER :: id_doc, id_don, id_dop
      INTEGER :: id_poc, id_pon, id_pop, id_fruit
      LOGICAL :: poc_sheet, pon_sheet, pop_sheet, fruit_sheet
      INTEGER :: id_sed_nup, id_sed_pup, id_sed_dic
      INTEGER :: id_sed_doc, id_sed_don, id_sed_dop
      INTEGER :: id_sed_poc, id_sed_pon, id_sed_pop
      AED_REAL :: sc_nup, sc_pup, sc_dic, sc_doc, sc_don, sc_dop  ! areal flux -> target units
      AED_REAL :: sc_poc, sc_pon, sc_pop
      LOGICAL  :: sed_nup_no3      ! sediment N uptake pool is nitrate (O2 from its reduction)

      !# Macrophyte simulation settings
      INTEGER  :: num_mphy
      INTEGER  :: n_zones
      LOGICAL  :: simMacFeedback, simStaticBiomass
      LOGICAL  :: simEpiphytes, simFruiting, simMacDrag
      INTEGER  :: drag_model, mac_pattern, mac_mode, mac_initial, resp_mode
      AED_REAL :: water_nutrient_frac, water_excr_frac
      AED_REAL :: active_layer_limit, f_mort_dom, sed_layer_dz
      AED_REAL,ALLOCATABLE :: active_zones(:)

      !# Macrophyte module species parameters
      TYPE(macrophyte_params_t),DIMENSION(:),ALLOCATABLE :: mpars
      TYPE(phyto_data_t),DIMENSION(:),ALLOCATABLE :: salpars  ! sal_method etc, for phyto_salinity

      !# Epiphyte sub-module option & parameters
      INTEGER  :: epi_model
      AED_REAL :: R_epig, R_epir, R_epib, I_Kepi, epi_max
      AED_REAL :: theta_epi_growth, theta_epi_resp
      AED_REAL :: epi_initial, epi_Xpc, epi_Xnc, epi_K_N, epi_K_P
      LOGICAL  :: do_Cuptake, do_Nuptake, do_Puptake

   CONTAINS      ! Selected AED methods to activate:
         PROCEDURE :: define             => aed_define_macrophyte
         PROCEDURE :: initialize_benthic => aed_initialize_benthic_macrophyte
         PROCEDURE :: calculate_column   => aed_calculate_column_macrophyte  ! step 1
         PROCEDURE :: calculate_benthic  => aed_calculate_benthic_macrophyte ! step 2
         PROCEDURE :: calculate          => aed_calculate_macrophyte         ! step 3
       ! PROCEDURE :: light_extinction   => aed_light_extinction_macrophyte
       ! PROCEDURE :: bio_drag           => aed_bio_drag_macrophyte
   END TYPE
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  MODULE WORKING VARIABLES

   INTEGER :: i
   AED_REAL, parameter :: h = 6.626e-34        ! Planck constant, Js
   AED_REAL, parameter :: c = 2.998e8          ! light speed, m/s
   AED_REAL, parameter :: Av= 6.02e23          ! Avagadro number, /mol, convert light capture to units of mol photon/m2/s
   AED_REAL, parameter :: WL_range(501) = (/ (I, I = 300, 800) /) ![300 800] ! wavelength range
   AED_REAL, parameter :: WL_mean = 550        ! mean wavelength, nm
   AED_REAL, ALLOCATABLE :: WLint(:), landaiint(:), ALlint(:)
   AED_REAL :: landaiSUM, ALl

   ! clear-sky wave lengths and irradiance at a particular wave-band W/m2/nm, from EMS model
   AED_REAL, parameter :: wavei(150) = (/ &
       140.00, 150.00, 160.00, 170.00, &
       180.00, 190.00, 200.00, 205.00, 210.00, 215.00, 220.00, &
       225.00, 230.00, 235.00, 240.00, 245.00, 250.00, 255.00, &
       260.00, 265.00, 270.00, 275.00, 280.00, 285.00, 290.00, &
       295.00, 300.00, 305.00, 310.00, 315.00, 320.00, 325.00, &
       330.00, 335.00, 340.00, 345.00, 350.00, 355.00, 360.00, &
       365.00, 370.00, 375.00, 380.00, 385.00, 390.00, 395.00, &
       400.00, 405.00, 410.00, 415.00, 420.00, 425.00, 430.00, &
       435.00, 440.00, 445.00, 450.00, 455.00, 460.00, 465.00, &
       470.00, 475.00, 480.00, 485.00, 490.00, 495.00, 500.00, &
       505.00, 510.00, 515.00, 520.00, 525.00, 530.00, 535.00, &
       540.00, 545.00, 550.00, 555.00, 560.00, 565.00, 570.00, &
       575.00, 580.00, 585.00, 590.00, 595.00, 600.00, 610.00, &
       620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, &
       690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, &
       800.00, 850.00, 900.00, 950.00, 1000.00, 1100.00, 1200.00, &
       1300.00, 1400.00, 1500.00, 1600.00, 1700.00, 1800.00, 1900.00, &
       2000.00, 2100.00, 2200.00, 2300.00, 2400.00, 2500.00, 2600.00, &
       2700.00, 2800.00, 2900.00, 3000.00, 3100.00, 3200.00, 3300.00, &
       3400.00, 3500.00, 3600.00, 3700.00, 3800.00, 3900.00, 4000.00, &
       4100.00, 4200.00, 4300.00, 4400.00, 4500.00, 4600.00, 4700.00, &
       4800.00, 4900.00, 5000.00, 6000.00, 7000.00, 8000.00   /)

   AED_REAL, parameter :: landai(150) = (/ &
       0.0000, 0.0001, 0.0000, 0.0004, &
       0.0009,0.0017,0.0030,0.0050,0.0100,0.0180,0.0300,0.0420,0.0520, &
       0.0540,0.0580,0.0640,0.0640,0.1000,0.1300,0.2000,0.2500,0.2200, &
       0.2400,0.3400,0.5200,0.6300,0.6100,0.6700,0.7600,0.8200,0.8500, &
       1.0200,1.1500,1.1100,1.1100,1.1700,1.1800,1.1600,1.1600,1.2900, &
       1.3300,1.3200,1.2300,1.1500,1.1200,1.2000,1.5400,1.8800,1.9400, &
       1.9200,1.9200,1.8900,1.7800,1.8200,2.0300,2.1500,2.2000,2.1900, &
       2.1600,2.1500,2.1700,2.2000,2.1600,2.0300,1.9900,2.0400,1.9800, &
       1.9700,1.9600,1.8900,1.8700,1.9200,1.9500,1.9700,1.9800,1.9800, &
       1.9500,1.9200,1.9000,1.8900,1.8700,1.8700,1.8700,1.8500,1.8400, &
       1.8300,1.8100,1.7700,1.7400,1.7000,1.6600,1.6200,1.5900,1.5500, &
       1.5100,1.4800,1.4400,1.4100,1.3700,1.3400,1.3000,1.2700,1.1270, &
       1.0030,0.8950,0.8030,0.7250,0.6060,0.5010,0.4060,0.3280,0.2670, &
       0.2200,0.1820,0.1520,0.1274,0.1079,0.0917,0.0785,0.0676,0.0585, &
       0.0509,0.0445,0.0390,0.0343,0.0303,0.0268,0.0230,0.0214,0.0191, &
       0.0171,0.0153,0.0139,0.0125,0.0114,0.0103,0.0095,0.0087,0.0080, &
       0.0073,0.0067,0.0061,0.0056,0.0051,0.0048,0.0044,0.0042,0.0021, &
       0.0012,0.0006 /)

   INTEGER, PARAMETER :: SUBMERGED = 1
   INTEGER, PARAMETER :: EMERGENT  = 2
   INTEGER, PARAMETER :: FLOATING  = 3
   INTEGER            :: diag_level = 10     ! 0 = no diagnostic outputs
                                             ! 1 = basic diagnostic outputs
                                             ! 2 = flux rates, and supporitng
                                             ! 3 = other metrics
                                             !10 = all debug & checking outputs
!-------------------------------------------------------------------------------

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
   INTEGER :: unit, nccols, ccol, dcol
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
         dcol = ccol-1
         md(dcol)%m_name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('growth_form')    ; md(dcol)%growth_form  = extract_integer(values(ccol))
            CASE ('zone_lock')      ; md(dcol)%zone_lock    = extract_integer(values(ccol))
            CASE ('m_initial')      ; md(dcol)%m_initial    = extract_double(values(ccol))
            CASE ('m_initial_minb') ; md(dcol)%m_initial_minb = extract_double(values(ccol))
            CASE ('m_initial_mind') ; md(dcol)%m_initial_mind = extract_double(values(ccol))
            CASE ('m_initial_d0')   ; md(dcol)%m_initial_d0 = extract_double(values(ccol))
            CASE ('m0')             ; md(dcol)%m0           = extract_double(values(ccol))

            CASE ('drag_model')     ; md(dcol)%drag_model   = extract_integer(values(ccol))
            CASE ('K_CD')           ; md(dcol)%K_CD         = extract_double(values(ccol))
            CASE ('shoot_diameter') ; md(dcol)%shoot_diameter = extract_double(values(ccol))
            CASE ('shoot_length')   ; md(dcol)%shoot_length = extract_double(values(ccol))
            CASE ('tissue_density') ; md(dcol)%tissue_density = extract_double(values(ccol))
            CASE ('elastic_modulus'); md(dcol)%elastic_modulus = extract_double(values(ccol))

            CASE ('shading_model')  ; md(dcol)%shading_model= extract_integer(values(ccol))
            CASE ('KeMAC')          ; md(dcol)%KeMAC        = extract_double(values(ccol))
            CASE ('kA')             ; md(dcol)%kA           = extract_double(values(ccol))
            CASE ('k_omega')        ; md(dcol)%k_omega      = extract_double(values(ccol))
            CASE ('sine_blade')     ; md(dcol)%sine_blade   = extract_double(values(ccol))
            CASE ('X_ldw')          ; md(dcol)%X_ldw        = extract_double(values(ccol))
            CASE ('K_ldw')          ; md(dcol)%K_ldw        = extract_double(values(ccol))
            CASE ('X_sdw')          ; md(dcol)%X_sdw        = extract_double(values(ccol))

            CASE ('R_growth')       ; md(dcol)%R_growth     = extract_double(values(ccol))
            CASE ('temp_method')    ; md(dcol)%temp_method  = extract_integer(values(ccol))
            CASE ('theta_growth')   ; md(dcol)%theta_growth = extract_double(values(ccol))
            CASE ('T_std')          ; md(dcol)%T_std        = extract_double(values(ccol))
            CASE ('T_opt')          ; md(dcol)%T_opt        = extract_double(values(ccol))
            CASE ('T_max')          ; md(dcol)%T_max        = extract_double(values(ccol))
            CASE ('light_model')    ; md(dcol)%light_model  = extract_integer(values(ccol))
            CASE ('I_K')            ; md(dcol)%I_K          = extract_double(values(ccol))
            CASE ('I_S')            ; md(dcol)%I_S          = extract_double(values(ccol))
            CASE ('f_pr')           ; md(dcol)%f_pr         = extract_double(values(ccol))
            CASE ('K_epi')          ; md(dcol)%K_epi        = extract_double(values(ccol))

            CASE ('R_resp')         ; md(dcol)%R_resp       = extract_double(values(ccol))
            CASE ('theta_resp')     ; md(dcol)%theta_resp   = extract_double(values(ccol))
            CASE ('E_comp')         ; md(dcol)%E_comp       = extract_double(values(ccol))
            CASE ('R_mort')         ; md(dcol)%R_mort       = extract_double(values(ccol))
            CASE ('sal_method')     ; md(dcol)%sal_method   = extract_integer(values(ccol))
            CASE ('S_bep')          ; md(dcol)%S_bep        = extract_double(values(ccol))
            CASE ('S_maxsp')        ; md(dcol)%S_maxsp      = extract_double(values(ccol))
            CASE ('S_opt')          ; md(dcol)%S_opt        = extract_double(values(ccol))

            CASE ('m_initial_bg')   ; md(dcol)%m_initial_bg = extract_double(values(ccol))
            CASE ('f_bg')           ; md(dcol)%f_bg         = extract_double(values(ccol))
            CASE ('R_resp_bg')      ; md(dcol)%R_resp_bg    = extract_double(values(ccol))
            CASE ('R_mort_bg')      ; md(dcol)%R_mort_bg    = extract_double(values(ccol))
            CASE ('tau_tran')       ; md(dcol)%tau_tran     = extract_double(values(ccol))
            CASE ('X_cchl')         ; md(dcol)%X_cchl       = extract_double(values(ccol))
            CASE ('X_cdw')          ; md(dcol)%X_cdw        = extract_double(values(ccol))
            CASE ('K_N')            ; md(dcol)%K_N          = extract_double(values(ccol))
            CASE ('X_ncon')         ; md(dcol)%X_ncon       = extract_double(values(ccol))
            CASE ('K_P')            ; md(dcol)%K_P          = extract_double(values(ccol))
            CASE ('X_pcon')         ; md(dcol)%X_pcon       = extract_double(values(ccol))
            CASE ('K_S')            ; md(dcol)%K_S          = extract_double(values(ccol))
            CASE ('K_NPP')          ; md(dcol)%K_NPP        = extract_double(values(ccol))

            CASE ('fruit_model')    ; md(dcol)%fruit_model  = extract_integer(values(ccol))
            CASE ('m_initial_f')    ; md(dcol)%m_initial_f  = extract_double(values(ccol))
            CASE ('tau_tran_fruit') ; md(dcol)%tau_tran_fruit = extract_double(values(ccol))
            CASE ('r_release')      ; md(dcol)%r_release    = extract_double(values(ccol))
            CASE ('f_seed')         ; md(dcol)%f_seed       = extract_double(values(ccol))
            CASE ('t_start_g')      ; md(dcol)%t_start_g    = extract_double(values(ccol))
            CASE ('t_dur_g')        ; md(dcol)%t_dur_g      = extract_double(values(ccol))
            CASE ('t_start_r')      ; md(dcol)%t_start_r    = extract_double(values(ccol))
            CASE ('t_dur_r')        ; md(dcol)%t_dur_r      = extract_double(values(ccol))

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
!
!-------------------------------------------------------------------------------
!BEGIN
    ALLOCATE(md(MAX_PHYTO_TYPES))
    md(:)%R_mort = zero_ ; md(:)%R_mort_bg = zero_   ! optional in the parameter file
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, md)
       CASE (NML_TYPE)
           open(NEWUNIT=tfil,file=dbase, status='OLD', iostat=status)
           IF (status /= 0) STOP 'Cannot open macrophyte_data namelist file for macrophytes'
           read(tfil,nml=macrophyte_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist macrophyte_data for macrophytes'

    data%num_mphy = count
    ALLOCATE(data%mpars(count))

    ALLOCATE(data%id_mphya(count))
    ALLOCATE(data%id_mphyb(count))
    ALLOCATE(data%id_mphyf(count))

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%mpars(i)%growth_form    = md(list(i))%growth_form
       data%mpars(i)%zone_lock      = md(list(i))%zone_lock
       data%mpars(i)%m_name         = md(list(i))%m_name
       data%mpars(i)%m_initial      = md(list(i))%m_initial
       data%mpars(i)%m_initial_minb = md(list(i))%m_initial_minb
       data%mpars(i)%m_initial_mind = md(list(i))%m_initial_mind
       data%mpars(i)%m_initial_d0   = md(list(i))%m_initial_d0
       data%mpars(i)%m0             = md(list(i))%m0

       data%mpars(i)%drag_model     = md(list(i))%drag_model
       data%mpars(i)%K_CD           = md(list(i))%K_CD
       data%mpars(i)%shoot_diameter = md(list(i))%shoot_diameter
       data%mpars(i)%shoot_length   = md(list(i))%shoot_length
       data%mpars(i)%tissue_density = md(list(i))%tissue_density
       data%mpars(i)%elastic_modulus= md(list(i))%elastic_modulus

       data%mpars(i)%shading_model  = md(list(i))%shading_model
       data%mpars(i)%KeMAC          = md(list(i))%KeMAC
       data%mpars(i)%kA             = md(list(i))%kA
       data%mpars(i)%k_omega        = md(list(i))%k_omega
       data%mpars(i)%sine_blade     = md(list(i))%sine_blade
       data%mpars(i)%X_ldw          = md(list(i))%X_ldw
       data%mpars(i)%K_ldw          = md(list(i))%K_ldw
       data%mpars(i)%X_sdw          = md(list(i))%X_sdw

       data%mpars(i)%R_growth       = md(list(i))%R_growth/secs_per_day
       data%mpars(i)%temp_method    = md(list(i))%temp_method
       data%mpars(i)%theta_growth   = md(list(i))%theta_growth
       data%mpars(i)%T_std          = md(list(i))%T_std
       data%mpars(i)%T_opt          = md(list(i))%T_opt
       data%mpars(i)%T_max          = md(list(i))%T_max
       data%mpars(i)%light_model    = md(list(i))%light_model
       data%mpars(i)%I_K            = md(list(i))%I_K
       data%mpars(i)%I_S            = md(list(i))%I_S
       data%mpars(i)%f_pr           = md(list(i))%f_pr
       data%mpars(i)%K_epi          = md(list(i))%K_epi

       data%mpars(i)%R_resp         = md(list(i))%R_resp/secs_per_day
       data%mpars(i)%theta_resp     = md(list(i))%theta_resp
       data%mpars(i)%E_comp         = md(list(i))%E_comp
       data%mpars(i)%R_mort         = md(list(i))%R_mort/secs_per_day
       data%mpars(i)%sal_method     = md(list(i))%sal_method
       data%mpars(i)%S_bep          = md(list(i))%S_bep
       data%mpars(i)%S_maxsp        = md(list(i))%S_maxsp
       data%mpars(i)%S_opt          = md(list(i))%S_opt

       data%mpars(i)%m_initial_bg   = md(list(i))%m_initial_bg
       data%mpars(i)%f_bg           = md(list(i))%f_bg
       data%mpars(i)%R_resp_bg      = md(list(i))%R_resp_bg/secs_per_day
       data%mpars(i)%R_mort_bg      = md(list(i))%R_mort_bg/secs_per_day
       data%mpars(i)%tau_tran       = md(list(i))%tau_tran/secs_per_day
       data%mpars(i)%X_cchl         = md(list(i))%X_cchl
       data%mpars(i)%X_cdw          = md(list(i))%X_cdw
       data%mpars(i)%K_N            = md(list(i))%K_N
       data%mpars(i)%X_ncon         = md(list(i))%X_ncon
       data%mpars(i)%K_P            = md(list(i))%K_P
       data%mpars(i)%X_pcon         = md(list(i))%X_pcon
       data%mpars(i)%K_S            = md(list(i))%K_S
       data%mpars(i)%K_NPP          = md(list(i))%K_NPP

       data%mpars(i)%fruit_model    = md(list(i))%fruit_model
       data%mpars(i)%m_initial_f    = md(list(i))%m_initial_f
       data%mpars(i)%tau_tran_fruit = md(list(i))%tau_tran_fruit/secs_per_day
       data%mpars(i)%r_release      = md(list(i))%r_release/secs_per_day
       data%mpars(i)%f_seed         = md(list(i))%f_seed
       data%mpars(i)%t_start_g      = md(list(i))%t_start_g
       data%mpars(i)%t_dur_g        = md(list(i))%t_dur_g
       data%mpars(i)%t_start_r      = md(list(i))%t_start_r
       data%mpars(i)%t_dur_r        = md(list(i))%t_dur_r
       IF (.NOT. data%simFruiting)  data%mpars(i)%fruit_model = 0


       ! Register group's above ground biomass as a state variable
       data%id_mphya(i) = aed_define_sheet_variable(                         &
                              TRIM(md(list(i))%m_name)//'_ag', 'mmol C/m2',  &
                               'macrophyte group above-ground biomass',      &
                              md(list(i))%m_initial,                         &
                              minimum=md(list(i))%m0)
       ! Register group's below ground biomass as a state variable
       data%id_mphyb(i) = aed_define_sheet_variable(                         &
                              TRIM(md(list(i))%m_name)//'_bg', 'mmol C/m2',  &
                              'macrophyte group below-ground biomass',       &
                              md(list(i))%m_initial_bg,                      &
                              minimum=md(list(i))%m0)
       ! Register group's fruit biomass as a state variable
       IF( data%mpars(i)%fruit_model >0 ) THEN
          data%id_mphyf(i) = aed_define_sheet_variable(                      &
                              TRIM(md(list(i))%m_name)//'_fr', 'mmol C/m2',  &
                               'macrophyte group seed/fruit biomass',        &
                              md(list(i))%m_initial_f,                       &
                              minimum=md(list(i))%m0/2.)
       ENDIF
!       ! Register group's shoot number as diagnostic
!       IF( data%mpars(i)%drag_model >0 ) THEN
!         data%id_shdens(i) = aed_define_sheet_diag_variable(                 &
!                             TRIM(md(list(i))%m_name)//'_shdens', '# sh/m2', &
!                              'macrophyte group shoot density')
!         data%id_shhght(i) = aed_define_sheet_diag_variable(                 &
!                             TRIM(md(list(i))%m_name)//'_shhght', '# m', &
!                              'macrophyte group shoot/leaf length')
!       ENDIF

   ENDDO
   DEALLOCATE(md)
END SUBROUTINE aed_macrophyte_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION mac_locate_link(name, is_sheet) RESULT(ret)
!-------------------------------------------------------------------------------
! Locate a link that may be a water column or a sediment (sheet) variable.
!  NB: a sediment variable is only recognised if its module is listed before
!      aed_macrophyte; otherwise it is taken to be a water column variable.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
   LOGICAL,INTENT(out) :: is_sheet
!
!LOCALS
   INTEGER :: ret
!
!-------------------------------------------------------------------------------
!BEGIN
   is_sheet = .FALSE. ; ret = 0
   IF ( name == '' ) RETURN

   ret = aed_locate_sheet_variable(name)
   is_sheet = ret > 0
   IF ( .NOT. is_sheet ) ret = aed_locate_variable(name)
END FUNCTION mac_locate_link
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION mac_locate_sed_link(name, dz, scale) RESULT(ret)
!-------------------------------------------------------------------------------
! Locate a sediment (sheet) variable link; it must already be defined, ie its
! module is listed before aed_macrophyte. Targets held as a concentration
! (units per m3, eg porewater nutrients) take the areal flux / dz; areal targets
! (per m2) take it as is.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
   AED_REAL,INTENT(in)  :: dz
   AED_REAL,INTENT(out) :: scale
!
!LOCALS
   INTEGER :: ret
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   ret = 0 ; scale = one_
   IF ( name == '' ) RETURN

   ret = aed_locate_sheet_variable(name)
   IF ( ret == 0 ) THEN
      print *,"          ERROR: sediment variable '",TRIM(name),"' not found;"
      print *,"          its module must be listed before aed_macrophyte"
      STOP 'aed_macrophyte: sediment link not found'
   ENDIF
   IF ( aed_get_var(ret, tvar) ) THEN
      IF ( INDEX(tvar%units,'/m3') > 0 .OR. INDEX(tvar%units,'m-3') > 0 ) scale = one_/dz
   ENDIF
END FUNCTION mac_locate_sed_link
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE mac_add_flux(column, layer_idx, id, is_sheet, flux)
!-------------------------------------------------------------------------------
! Add an areal flux (mmol/m2/s) to a linked water column or sediment variable
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx, id
   LOGICAL,INTENT(in) :: is_sheet
   AED_REAL,INTENT(in) :: flux
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( id <= 0 ) RETURN
   IF ( is_sheet ) THEN
      _FLUX_VAR_B_(id) = _FLUX_VAR_B_(id) + flux
   ELSE
      _FLUX_VAR_(id) = _FLUX_VAR_(id) + flux
   ENDIF
END SUBROUTINE mac_add_flux
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE mac_nutrient_supply(wnf, wat_linked, c_w, sed_linked, c_s, K, s_w, a_s, f)
!-------------------------------------------------------------------------------
! Supply of one nutrient (N or P) for growth, as fractions of the demand:
!  s_w : met by the water - its share wnf, less as the water runs low
!  a_s : of the remainder, met by the roots - from a linked sediment pool, or
!        an assumed ample supply if no pool is linked; none if wnf = 1
!  f   : availability, s_w + (1-s_w)*a_s, ie the growth limitation
! With K = 0 there is no limitation: the fixed shares wnf, 1-wnf.
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in)  :: wnf, c_w, c_s, K
   LOGICAL,INTENT(in)   :: wat_linked, sed_linked
   AED_REAL,INTENT(out) :: s_w, a_s, f
!
!-------------------------------------------------------------------------------
!BEGIN
   s_w = wnf
   IF ( wat_linked .AND. K > zero_ ) s_w = wnf * MAX(c_w,zero_) / (MAX(c_w,zero_) + K)

   IF ( wnf >= one_ ) THEN
      a_s = zero_
   ELSEIF ( sed_linked .AND. K > zero_ ) THEN
      a_s = MAX(c_s,zero_) / (MAX(c_s,zero_) + K)
   ELSE
      a_s = one_
   ENDIF

   f = s_w + (one_-s_w) * a_s
END SUBROUTINE mac_nutrient_supply
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION mac_is_nitrate(name)
!-------------------------------------------------------------------------------
! Whether a linked N pool is nitrate, from its name, eg NIT_nit, SFM_no3_2
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: name
!
!LOCALS
   CHARACTER(len=LEN(name)) :: lc
   INTEGER :: i, u
!
!-------------------------------------------------------------------------------
!BEGIN
   lc = name
   DO i=1,LEN_TRIM(lc)
      IF ( lc(i:i) >= 'A' .AND. lc(i:i) <= 'Z' ) lc(i:i) = ACHAR(IACHAR(lc(i:i))+32)
   ENDDO
   u = INDEX(lc,'_')                       ! skip the module prefix
   lc = lc(u+1:)
   mac_is_nitrate = lc(1:3) == 'no3' .OR. lc(1:3) == 'nox' .OR. TRIM(lc) == 'nit'
END FUNCTION mac_is_nitrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION mac_in_zone(data, matz)
!-------------------------------------------------------------------------------
! Whether plants can grow in sediment zone matz: n_zones = 0 means all zones
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   AED_REAL,INTENT(in) :: matz
!
!-------------------------------------------------------------------------------
!BEGIN
   mac_in_zone = .TRUE.
   IF ( data%n_zones > 0 ) mac_in_zone = in_zone_set(matz, data%active_zones)
END FUNCTION mac_in_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE mac_shoots(data, mi, MAC_A, n_shoot, sh_height, sh_diameter, biovolume)
!-------------------------------------------------------------------------------
! Shoot density, height, diameter and biovolume of group mi from its leaf C
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   INTEGER,INTENT(in)   :: mi
   AED_REAL,INTENT(in)  :: MAC_A                          ! mmol C/m2
   AED_REAL,INTENT(out) :: n_shoot, sh_height, sh_diameter, biovolume
!
!LOCALS
   AED_REAL :: dw
   AED_REAL,PARAMETER :: rho_v = 150e3  !!CHECK, unit in g DW/m3
!
!-------------------------------------------------------------------------------
!BEGIN
   dw = MAX(MAC_A, zero_) / data%mpars(mi)%X_cdw
   n_shoot = dw * data%mpars(mi)%X_sdw
   sh_height =  max(0.05, data%mpars(mi)%X_ldw * (dw/(dw+120.)))   ! set min height 5cm to avoid unrealistic large diameter
   sh_diameter = zero_
   IF ( n_shoot > zero_ ) sh_diameter = 2.* sqrt( (dw/rho_v/(3.142*sh_height*n_shoot)) )
   biovolume = n_shoot* sh_height * 3.142*(sh_diameter/2.)**2.
END SUBROUTINE mac_shoots
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE mac_canopy(data, column, layer_idx, matz, sh_dens, sh_diam, height, biovolume)
!-------------------------------------------------------------------------------
! Canopy geometry from the current leaf C of all groups, as step 2b sets the
! canopy diagnostics. Hosts zero diagnostics each time step, so steps 1 and 2a
! must compute it rather than read last step's diagnostics.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in)   :: layer_idx
   AED_REAL,INTENT(in)  :: matz
   AED_REAL,INTENT(out) :: sh_dens, sh_diam, height, biovolume
!
!LOCALS
   INTEGER  :: mi
   AED_REAL :: n_shoot, sh_height, sh_diameter, bv
!
!-------------------------------------------------------------------------------
!BEGIN
   sh_dens = zero_ ; sh_diam = zero_ ; height = zero_ ; biovolume = zero_
   DO mi=1,data%num_mphy
      IF( data%mpars(mi)%zone_lock >0 .and. data%mpars(mi)%zone_lock .ne. INT(matz) ) CYCLE
      CALL mac_shoots(data, mi, _STATE_VAR_S_(data%id_mphya(mi)), n_shoot, sh_height, sh_diameter, bv)
      sh_dens = sh_dens + n_shoot
      sh_diam = sh_diam + sh_diameter
      height  = height + sh_height
      biovolume = biovolume + bv
   ENDDO
   sh_diam = sh_diam / REAL(data%num_mphy)
   height  = height / REAL(data%num_mphy)
END SUBROUTINE mac_canopy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION mac_epi_prod(data, column, layer_idx, par, extc, Io, dz, temp, epi, epi_max) RESULT(prod)
!-------------------------------------------------------------------------------
! Epiphyte gross production rate (/s), with the light, temperature and
! nutrients of a water layer (an unlinked nutrient is taken as not limiting).
! epi, epi_max are for the same area (benthos), so the crowding is the same in
! each layer of the canopy.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: par, extc, Io, dz, temp, epi, epi_max
   AED_REAL :: prod
!
!LOCALS
   AED_REAL :: din, dip, fN, fP, fI
!
!-------------------------------------------------------------------------------
!BEGIN
   prod = zero_
   IF ( epi_max <= zero_ ) RETURN       ! no leaves to grow on

   fN = one_ ; fP = one_
   IF ( data%do_Nuptake .AND. data%epi_K_N > zero_ ) THEN
      din = zero_
      IF (data%id_nox>0) din = din + MAX(_STATE_VAR_(data%id_nox), zero_)
      IF (data%id_nh4>0) din = din + MAX(_STATE_VAR_(data%id_nh4), zero_)
      fN = din/(din+data%epi_K_N)
   ENDIF
   IF ( data%do_Puptake .AND. data%epi_K_P > zero_ ) THEN
      dip = MAX(_STATE_VAR_(data%id_po4), zero_)
      fP = dip/(dip+data%epi_K_P)
   ENDIF

   fI = photosynthesis_irradiance(10,data%I_Kepi,data%I_Kepi,par,extc,Io,dz)
   prod = data%R_epig*fI*(data%theta_epi_growth**(temp-20.))*MIN(fN,fP) &
        * (one_-(MIN(epi,epi_max)/epi_max))
END FUNCTION mac_epi_prod
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE mac_epi_exchange(data, column, layer_idx, net)
!-------------------------------------------------------------------------------
! Epiphyte metabolism exchange with a water layer, for net production net
! (mmol C/m2/s from the benthic step, mmol C/m3/s from the pelagic one).
! Net growth takes DIC, N and P (ammonium preferred by concentration, 2 O2 per
! nitrate N); net respiration returns DIC and DON/DOP (no mineralisation).
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(in) :: net
!
!LOCALS
   AED_REAL :: nh4, no3, pNH4
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (data%id_oxy>0) _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + net
   IF (data%id_dic>0) _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - net
   IF ( net > zero_ ) THEN
      ! uptake, N split by concentration
      nh4 = zero_ ; no3 = zero_
      IF (data%id_nh4>0) nh4 = MAX(_STATE_VAR_(data%id_nh4), zero_)
      IF (data%id_nox>0) no3 = MAX(_STATE_VAR_(data%id_nox), zero_)
      pNH4 = 0.5
      IF (nh4+no3 > zero_) pNH4 = nh4/(nh4+no3)
      IF (data%id_nox==0) pNH4 = one_
      IF (data%id_nh4==0) pNH4 = zero_
      IF (data%id_nox>0) _FLUX_VAR_(data%id_nox) = &
               _FLUX_VAR_(data%id_nox) - net * data%epi_Xnc * (one_-pNH4)
      IF (data%id_nh4>0) _FLUX_VAR_(data%id_nh4) = &
               _FLUX_VAR_(data%id_nh4) - net * data%epi_Xnc * pNH4
      IF (data%id_nox>0 .AND. data%id_oxy>0) _FLUX_VAR_(data%id_oxy) = &
               _FLUX_VAR_(data%id_oxy) + 2. * net * data%epi_Xnc * (one_-pNH4)
      IF (data%id_po4>0) _FLUX_VAR_(data%id_po4) = &
               _FLUX_VAR_(data%id_po4) - net * data%epi_Xpc
   ELSE
      ! net respiration, N and P released as DON/DOP
      IF (data%id_don>0) _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) - net * data%epi_Xnc
      IF (data%id_dop>0) _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) - net * data%epi_Xpc
   ENDIF
END SUBROUTINE mac_epi_exchange
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
   INTEGER  :: status, i, ww
   CHARACTER(len=20)  :: sal_effect
   AED_REAL :: landaimid

   AED_REAL :: Wli(8), ALli(8)

!  %% NAMELIST   %%  /aed_macrophyte/
!  %% Last Checked 20/08/2021
   INTEGER            :: num_mphy = 0
   INTEGER            :: the_mphy(MAX_PHYTO_TYPES) = 0
   CHARACTER(len=128) :: dbase = 'aed_macrophyte_pars.nml'

   INTEGER            :: n_zones = 0
   INTEGER            :: active_zones(MAX_ZONES)

   LOGICAL            :: simStaticBiomass = .FALSE.
   LOGICAL            :: simMacFeedback = .FALSE.
   LOGICAL            :: simEpiphytes = .FALSE.
   LOGICAL            :: simFruiting = .FALSE.
   LOGICAL            :: simMacDrag = .FALSE.

   INTEGER            :: mac_initial = 0
   INTEGER            :: resp_mode = 1   ! light_model 2,3 respiration: 0 = daylight only
                                         ! (Baird et al. 2016, EMS), 1 = over 24 h

   AED_REAL           :: water_nutrient_frac = zero_
   AED_REAL           :: water_excr_frac = zero_
   AED_REAL           :: active_layer_limit = 0.05

   CHARACTER(len=64)  :: do_uptake_target_variable = 'OXY_oxy'
   CHARACTER(len=64)  :: c_uptake_target_variable = ''
   CHARACTER(len=64)  :: n1_uptake_target_variable = 'NIT_nit'
   CHARACTER(len=64)  :: n2_uptake_target_variable = 'NIT_amm'
   CHARACTER(len=64)  :: p1_uptake_target_variable = 'PHS_frp'
   CHARACTER(len=64)  :: n_excretion_target_variable = 'OGM_don'
   CHARACTER(len=64)  :: p_excretion_target_variable = 'OGM_dop'
   CHARACTER(len=64)  :: c_excretion_target_variable = ''
   CHARACTER(len=64)  :: c_mortality_target_variable = ''
   CHARACTER(len=64)  :: n_mortality_target_variable = ''
   CHARACTER(len=64)  :: p_mortality_target_variable = ''
   CHARACTER(len=64)  :: fruit_target_variable = ''
   AED_REAL           :: f_mort_dom = zero_

   ! sediment pools (sheet variables), receiving the parts not exchanged with
   ! the water: (1-water_nutrient_frac) of uptake, (1-water_excr_frac) of release
   CHARACTER(len=64)  :: n_sed_uptake_target_variable = ''
   CHARACTER(len=64)  :: p_sed_uptake_target_variable = ''
   CHARACTER(len=64)  :: c_sed_uptake_target_variable = ''
   CHARACTER(len=64)  :: c_sed_excretion_target_variable = ''
   CHARACTER(len=64)  :: n_sed_excretion_target_variable = ''
   CHARACTER(len=64)  :: p_sed_excretion_target_variable = ''
   CHARACTER(len=64)  :: c_sed_mortality_target_variable = ''
   CHARACTER(len=64)  :: n_sed_mortality_target_variable = ''
   CHARACTER(len=64)  :: p_sed_mortality_target_variable = ''
   AED_REAL           :: sed_layer_dz = one_   ! (m) sediment targets with units per m3
                                               ! get the areal flux / sed_layer_dz

   !AED_REAL           :: bg_gpp_frac = 0.1
   !AED_REAL           :: coef_bm_hgt = 1e-5

   INTEGER            :: epi_model = 1
   AED_REAL           :: epi_initial = 0.1
   AED_REAL           :: R_epig = zero_
   AED_REAL           :: R_epir = zero_
   AED_REAL           :: R_epib = zero_
   AED_REAL           :: I_Kepi = 100.
   AED_REAL           :: epi_Xnc = 16./106.
   AED_REAL           :: epi_Xpc = 1./106.
   AED_REAL           :: epi_max = 10000.
   AED_REAL           :: theta_epi_growth = 1.06
   AED_REAL           :: theta_epi_resp = 1.06
   AED_REAL           :: epi_K_N = 10.
   AED_REAL           :: epi_K_P = 1.0

!  %% END NAMELIST   %%  /aed_macrophyte/

   NAMELIST /aed_macrophyte/ num_mphy, the_mphy, dbase,               &
                             n_zones, active_zones,                   &
                             simMacFeedback, simStaticBiomass,        &
                             simEpiphytes, simFruiting, simMacDrag,   &
                             water_nutrient_frac,water_excr_frac,     &
                             active_layer_limit,                      &
                             do_uptake_target_variable,               &
                             c_uptake_target_variable,                &
                             n1_uptake_target_variable,               &
                             n2_uptake_target_variable,               &
                             p1_uptake_target_variable,               &
                             n_excretion_target_variable,             &
                             p_excretion_target_variable,             &
                             c_excretion_target_variable,             &
                             c_mortality_target_variable,             &
                             n_mortality_target_variable,             &
                             p_mortality_target_variable,             &
                             fruit_target_variable, f_mort_dom,       &
                             n_sed_uptake_target_variable,            &
                             p_sed_uptake_target_variable,            &
                             c_sed_uptake_target_variable,            &
                             c_sed_excretion_target_variable,         &
                             n_sed_excretion_target_variable,         &
                             p_sed_excretion_target_variable,         &
                             c_sed_mortality_target_variable,         &
                             n_sed_mortality_target_variable,         &
                             p_sed_mortality_target_variable,         &
                             sed_layer_dz,                            &
                             epi_model,R_epig,R_epib,R_epir,I_Kepi,epi_max,  &
                             epi_Xnc,epi_Xpc,epi_K_N,epi_K_P,         &
                             theta_epi_growth,theta_epi_resp,         &
                             mac_initial, epi_initial, diag_level,    &
                             resp_mode

!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_macrophyte configuration"

   ! Read the namelist
   read(namlst,nml=aed_macrophyte,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_macrophyte'

   data%simStaticBiomass = simStaticBiomass
   data%simMacFeedback = simMacFeedback
   data%simEpiphytes = simEpiphytes
   data%simFruiting = simFruiting
   data%simMacDrag = simMacDrag
   data%water_nutrient_frac = water_nutrient_frac
   data%water_excr_frac = water_excr_frac
   data%active_layer_limit = active_layer_limit
   data%f_mort_dom = MIN(MAX(f_mort_dom,zero_),one_)
   IF ( sed_layer_dz <= zero_ ) STOP 'aed_macrophyte: sed_layer_dz must be > 0'
   data%sed_layer_dz = sed_layer_dz

   ! mac_initial: -1 keep the host's values, 0 the dbase values, 1 from depth
   IF ( mac_initial == 4 ) THEN
      print *,"          NOTE: mac_initial 4 is now mac_initial -1 (keep the host's initial values)"
      mac_initial = -1
   ENDIF
   IF ( mac_initial == 2 .OR. mac_initial == 3 ) &
      STOP 'aed_macrophyte: mac_initial 2 (LAI) and 3 (canopy inputs) are not supported; use -1, 0 or 1'
   IF ( mac_initial < -1 .OR. mac_initial > 1 ) &
      STOP 'aed_macrophyte: mac_initial must be -1 (host), 0 (dbase) or 1 (depth)'
   data%mac_initial = mac_initial
   data%resp_mode = resp_mode
   data%epi_initial = epi_initial

   IF( .NOT. simEpiphytes) epi_model = 0; data%epi_model = epi_model
   IF( epi_model < 0 .OR. epi_model > 2 ) &
      STOP 'aed_macrophyte: epi_model must be 1 (bulk) or 2 (resolved over the canopy layers)'
   IF( epi_model == 2 .AND. data%aed_model_zone_avg ) THEN
      print *,"          NOTE: epi_model 2 can't be used with sediment zones (zone averaging) - using epi_model 1"
      epi_model = 1 ; data%epi_model = 1
   ENDIF
   data%theta_epi_growth = theta_epi_growth ; data%theta_epi_resp = theta_epi_resp
   data%epi_Xpc = epi_Xpc ; data%epi_K_P = epi_K_P
   data%epi_K_N = epi_K_N ; data%epi_Xnc = epi_Xnc
   data%I_Kepi = I_Kepi ; data%epi_max = epi_max
   data%R_epig = R_epig /secs_per_day
   data%R_epib = R_epib /secs_per_day
   data%R_epir = R_epir /secs_per_day

  ! !remove
  ! data%bg_gpp_frac = bg_gpp_frac  ! make this species specific
  ! data%coef_bm_hgt = coef_bm_hgt  ! make this species specific

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   ! Store parameter values in modules own data structure (macrophyte_params_t)
   !   NB: all rates are provided in values per day,
   !       and are converted here to values per second for internal AED use.
   !  The requested macrophyte state variables are registered in this function
   CALL aed_macrophyte_load_params(data, dbase, num_mphy, the_mphy)

   CALL aed_bio_temp_function(  data%num_mphy,               &
                                data%mpars%theta_growth,     &
                                data%mpars%T_std,            &
                                data%mpars%T_opt,            &
                                data%mpars%T_max,            &
                                data%mpars%aTn,              &
                                data%mpars%bTn,              &
                                data%mpars%kTn,              &
                                data%mpars%m_name)

   ! Temperature and salinity limitation, as for PHY and MAG: the salinity
   ! settings are copied to the phytoplankton form used by phyto_salinity
   ALLOCATE(data%salpars(data%num_mphy))
   ALLOCATE(data%id_fT(data%num_mphy))   ; data%id_fT(:) = 0
   ALLOCATE(data%id_fSal(data%num_mphy)) ; data%id_fSal(:) = 0
   DO i=1,data%num_mphy
      data%salpars(i)%salTol  = data%mpars(i)%sal_method
      data%salpars(i)%S_bep   = data%mpars(i)%S_bep
      data%salpars(i)%S_maxsp = data%mpars(i)%S_maxsp
      data%salpars(i)%S_opt   = data%mpars(i)%S_opt
      IF ( data%mpars(i)%sal_method < 0 .AND. data%mpars(i)%S_bep > one_ ) THEN
         print *,"          ERROR: ",TRIM(data%mpars(i)%m_name),": sal_method < 0 limits growth, so S_bep must be <= 1"
         STOP
      ENDIF
      IF ( data%mac_initial == 1 .AND. data%mpars(i)%m_initial_mind <= data%mpars(i)%m_initial_d0 ) THEN
         print *,"          ERROR: ",TRIM(data%mpars(i)%m_name),": mac_initial 1 needs m_initial_mind > m_initial_d0"
         STOP
      ENDIF
      sal_effect = ' (none)'
      IF ( data%mpars(i)%sal_method < 0 ) sal_effect = ' (on growth)'
      IF ( data%mpars(i)%sal_method > 0 ) sal_effect = ' (on respiration)'
      print "(11X,A,': temp_method ',I2,', sal_method ',I2,A)", TRIM(data%mpars(i)%m_name), &
            data%mpars(i)%temp_method, data%mpars(i)%sal_method, TRIM(sal_effect)
      IF ( diag_level > 1 ) THEN
         data%id_fT(i)   = aed_define_sheet_diag_variable(TRIM(data%mpars(i)%m_name)//'_fT', '-', &
                                                          'macrophyte temperature function')
         data%id_fSal(i) = aed_define_sheet_diag_variable(TRIM(data%mpars(i)%m_name)//'_fSal', '-', &
                                                          'macrophyte salinity function')
      ENDIF
   ENDDO

   ! fruit_model is zeroed when simFruiting is off; equally, there is no
   ! fruiting if no species has a fruit_model
   data%simFruiting = ANY(data%mpars(:)%fruit_model > 0)
   IF ( simFruiting .AND. .NOT. data%simFruiting ) &
      print *,"          WARNING: simFruiting is on but no species has fruit_model > 0"

   ! Variables related to epiphytes
   IF ( simEpiphytes ) THEN
      print *,"          epiphyte simulation is included: epi_model = ",epi_model
      data%id_epi = aed_define_sheet_variable( 'epiphyte', 'mmol C/m2',       &
                                               'epiphyte biomass',            &
                                                epi_initial, 0.001            )
     !data%id_epib = aed_define_sheet_diag_variable('epi_ben','mmol C/m2','total epiphyte biomass')
      IF ( diag_level > 0) THEN
        data%id_epib = aed_define_sheet_diag_variable('epi_ben','mmol C/m2','total epiphyte biomass')
        IF ( diag_level > 1) THEN
         data%id_epig = aed_define_sheet_diag_variable('epi_gpp','mmol C/m2/d','benthic gross productivity')
         data%id_epir = aed_define_sheet_diag_variable('epi_rsp','mmol C/m2/d','benthic phyto respiration')
         data%id_epil = aed_define_sheet_diag_variable('epi_lss','mmol C/m2/d','epiphyte loss (mortality, sloughing)')
        ENDIF
       ENDIF
      ! epi_model 2: set in step 1 (calculate_column), used in steps 2c and 3
      IF ( epi_model == 2 ) THEN
        data%id_epi_npp = aed_define_diag_variable('epi_npp','mmol C/m3/d','epiphyte net production, in each layer')
        data%id_epi_net = aed_define_sheet_diag_variable('epi_net','mmol C/m2/d','epiphyte net production, column total')
        data%id_epi_col = aed_define_sheet_diag_variable('epi_col','-','epiphyte column step: 0 none, 1 done, 2 used')
        print *,"          NOTE: epi_model 2 needs a host that calls calculate_column (TUFLOW-FV, SCHISM,"
        print *,"                GLM, zeroD); otherwise each cell runs as epi_model 1"
      ENDIF
   ENDIF


   ! Register diagnostic variables, for macrophyte the community canopy
   data%id_d_par  = aed_define_sheet_diag_variable('par','W/m2','light intensity at canopy top')
   data%id_gpp    = aed_define_sheet_diag_variable('npp','mmol C/m2/d','macrophyte net productivity')
  !data%id_gpp    = aed_define_sheet_diag_variable('gpp','mmol C/m2/d','benthic plant productivity')
   data%id_mac    = aed_define_sheet_diag_variable('mac_ben','mmol C/m2','total macrophyte biomass')
   data%id_canopy_lai = aed_define_sheet_diag_variable('mac_lai','m2/m2','macrophyte leaf area density')
   data%id_mac_ag = aed_define_sheet_diag_variable('mac_ag','mmol C/m2','total above ground macrophyte biomass')
   data%id_mac_bg = aed_define_sheet_diag_variable('mac_bg','mmol C/m2','total below ground macrophyte biomass')
   data%id_root_d = aed_define_sheet_diag_variable('mac_root_depth','m','mean depth of roots below the sediment surface')
   data%id_root_o = aed_define_sheet_diag_variable('mac_root_o2','mmol O2/m2','mean O2 injection rate of roots into sediment')
   data%id_mac_tr = aed_define_sheet_diag_variable('trn_bg','mmol C/m2/day', 'macrophyte biomass translocation')

   ! Variables related to fruiting
   IF ( data%simFruiting ) THEN
      print *,"          fruiting simulation is included"
      data%id_mac_fr = aed_define_sheet_diag_variable('mac_fr','mmol C/m2', 'macrophyte community fruit biomass')
      data%id_mac_ft = aed_define_sheet_diag_variable('trn_fr','mmol C/m2/day', 'macrophyte fruit biomass translocation')
      data%id_mac_gt = aed_define_sheet_diag_variable('tri_fr','-', 'macrophyte fruit growth trigger', rezero=.FALSE.)
   ENDIF

   ! Variables for macrophyte canopy light attenuation
   data%id_kemac = aed_define_diag_variable('ke_mac','/m', 'macrophyte light attenuation')
   IF (ANY(data%mpars(:)%light_model == 2) .or. ANY(data%mpars(:)%light_model == 3)) THEN
      IF ( resp_mode == 0 ) THEN
         print *,"          resp_mode = 0: light_model 2,3 respiration in daylight only"
      ELSE
         print *,"          resp_mode = 1: light_model 2,3 respiration over 24 h"
      ENDIF
      ! For spectral light capture
      data%id_parc = aed_define_sheet_diag_variable('mac_ag_parc','', 'macrophyte leaf photon capture')
      data%id_aeff = aed_define_sheet_diag_variable('mac_ag_aeff','', 'macrophyte leaf effective area')
   ENDIF

   ! Variables for drag feedback, computed based on macrophyte community canopy properties
   data%id_canopy_height  = aed_define_sheet_diag_variable('canopy_height','m', '')
   data%id_canopy_sh_dens = aed_define_sheet_diag_variable('canopy_shoot_dens','shoot/m2', '')
   data%id_canopy_sh_diam = aed_define_sheet_diag_variable('canopy_shoot_diam','m/shoot', '')
!  IF (data%drag_model > 1) THEN
      data%id_canopy_drag      = aed_define_sheet_diag_variable('canopy_drag','-', '')       !2D becoming 3D
      data%id_canopy_velocity= aed_define_sheet_diag_variable('canopy_velocity','m/s', '')
      data%id_canopy_biovolume = aed_define_sheet_diag_variable('canopy_biovolume','m3 leaf/m2', '')
      data%id_canopy_blockage  = aed_define_diag_variable('canopy_blockage','m3/m2', '')
      data%id_canopy_frarea    = aed_define_diag_variable('canopy_frarea','m2/m2?', '')
!  ENDIF

   ! Link to variables from other modules
   !  (set a target to '' in the namelist to switch that exchange off)
   data%id_oxy = 0 ; data%id_dic = 0
   data%id_nox = 0 ; data%id_nh4 = 0 ; data%id_po4 = 0
   data%id_doc = 0 ; data%id_don = 0 ; data%id_dop = 0
   data%id_poc = 0 ; data%id_pon = 0 ; data%id_pop = 0 ; data%id_fruit = 0
   data%poc_sheet = .FALSE. ; data%pon_sheet = .FALSE.
   data%pop_sheet = .FALSE. ; data%fruit_sheet = .FALSE.
   IF ( simMacFeedback .or. simEpiphytes ) THEN
      IF (do_uptake_target_variable .NE. '') &
         data%id_oxy = aed_locate_variable(do_uptake_target_variable)
      IF (c_uptake_target_variable .NE. '') &
         data%id_dic = aed_locate_variable(c_uptake_target_variable)
      IF (n1_uptake_target_variable .NE. '') &
         data%id_nox = aed_locate_variable(n1_uptake_target_variable)
      IF (n2_uptake_target_variable .NE. '') &
         data%id_nh4 = aed_locate_variable(n2_uptake_target_variable)
      IF (p1_uptake_target_variable .NE. '') &
         data%id_po4 = aed_locate_variable(p1_uptake_target_variable)
      IF (n_excretion_target_variable .NE. '') &
         data%id_don = aed_locate_variable(n_excretion_target_variable)
      IF (p_excretion_target_variable .NE. '') &
         data%id_dop = aed_locate_variable(p_excretion_target_variable)
   ENDIF
   IF ( simMacFeedback .AND. .NOT. simStaticBiomass ) THEN
      ! plant losses (mortality, fruit release) only leave a dynamic plant
      IF (c_excretion_target_variable .NE. '') &
         data%id_doc = aed_locate_variable(c_excretion_target_variable)
      ! leaf mortality and fruit may go to a water (eg OGM_poc) or a
      ! sediment (eg SFM_pom_c_g1) pool
      data%id_poc = mac_locate_link(c_mortality_target_variable, data%poc_sheet)
      data%id_pon = mac_locate_link(n_mortality_target_variable, data%pon_sheet)
      data%id_pop = mac_locate_link(p_mortality_target_variable, data%pop_sheet)
      IF (data%simFruiting) &
         data%id_fruit = mac_locate_link(fruit_target_variable, data%fruit_sheet)
   ENDIF
   data%id_sed_nup = 0 ; data%id_sed_pup = 0 ; data%id_sed_dic = 0
   data%id_sed_doc = 0 ; data%id_sed_don = 0 ; data%id_sed_dop = 0
   data%id_sed_poc = 0 ; data%id_sed_pon = 0 ; data%id_sed_pop = 0
   data%sc_nup = one_ ; data%sc_pup = one_ ; data%sc_dic = one_
   data%sc_doc = one_ ; data%sc_don = one_ ; data%sc_dop = one_
   data%sc_poc = one_ ; data%sc_pon = one_ ; data%sc_pop = one_
   data%sed_nup_no3 = .FALSE.
   IF ( simMacFeedback ) THEN
      ! sediment pools (eg SFM_nh4_2, SFM_pom_c_g1), as sheet variables
      data%id_sed_nup = mac_locate_sed_link(n_sed_uptake_target_variable, data%sed_layer_dz, data%sc_nup)
      data%sed_nup_no3 = data%id_sed_nup > 0 .AND. mac_is_nitrate(n_sed_uptake_target_variable)
      data%id_sed_pup = mac_locate_sed_link(p_sed_uptake_target_variable, data%sed_layer_dz, data%sc_pup)
      data%id_sed_dic = mac_locate_sed_link(c_sed_uptake_target_variable, data%sed_layer_dz, data%sc_dic)
      data%id_sed_doc = mac_locate_sed_link(c_sed_excretion_target_variable, data%sed_layer_dz, data%sc_doc)
      data%id_sed_don = mac_locate_sed_link(n_sed_excretion_target_variable, data%sed_layer_dz, data%sc_don)
      data%id_sed_dop = mac_locate_sed_link(p_sed_excretion_target_variable, data%sed_layer_dz, data%sc_dop)
      data%id_sed_poc = mac_locate_sed_link(c_sed_mortality_target_variable, data%sed_layer_dz, data%sc_poc)
      data%id_sed_pon = mac_locate_sed_link(n_sed_mortality_target_variable, data%sed_layer_dz, data%sc_pon)
      data%id_sed_pop = mac_locate_sed_link(p_sed_mortality_target_variable, data%sed_layer_dz, data%sc_pop)
   ENDIF
   data%do_Cuptake = data%id_dic > 0
   data%do_Nuptake = data%id_nox > 0 .OR. data%id_nh4 > 0
   data%do_Puptake = data%id_po4 > 0

   ! Report where the mass budget is open, given the options and links
   IF ( simMacFeedback ) THEN
      IF ( ANY(data%mpars(:)%K_N > zero_) .OR. ANY(data%mpars(:)%K_P > zero_) ) THEN
         IF ( water_nutrient_frac >= one_ ) THEN
            print *,"          NOTE: water_nutrient_frac = 1 (no root uptake) - growth is limited by water N,P"
         ELSEIF ( data%id_sed_nup == 0 .AND. data%id_sed_pup == 0 ) THEN
            print *,"          NOTE: no sediment N,P pools linked - roots have an ample supply, so growth is"
            print *,"                not nutrient limited; water uptake eases off as the water runs low"
         ELSE
            print *,"          NOTE: growth is nutrient limited by the water and the linked sediment pools"
         ENDIF
      ENDIF
      IF ( simStaticBiomass ) &
         print *,"          NOTE: simStaticBiomass - water fluxes are applied but the plant does not change"
      IF ( data%id_dic == 0 ) &
         print *,"          NOTE: no c_uptake_target_variable - plant carbon exchange is not tracked"
      IF ( water_nutrient_frac < one_ .AND. (data%id_sed_nup == 0 .OR. data%id_sed_pup == 0) ) &
         print *,"          NOTE: water_nutrient_frac < 1 - sediment N,P uptake is not (fully) linked"
      IF ( water_excr_frac < one_ .AND. (data%id_sed_don == 0 .OR. data%id_sed_dop == 0) ) &
         print *,"          NOTE: water_excr_frac < 1 - sediment N,P release is not (fully) linked"
      IF ( water_excr_frac < one_ .AND. .NOT. simStaticBiomass ) THEN
         IF ( ANY(data%mpars(:)%light_model /= 10 .AND. data%mpars(:)%R_resp_bg > zero_) &
              .AND. data%id_sed_dic == 0 ) &
            print *,"          NOTE: no c_sed_uptake_target_variable - root respiration C is lost"
         IF ( ANY(data%mpars(:)%R_mort_bg > zero_) ) THEN
            IF ( data%f_mort_dom < one_ .AND. data%id_sed_poc == 0 ) &
               print *,"          NOTE: no c_sed_mortality_target_variable - root POM is lost"
            IF ( data%f_mort_dom > zero_ .AND. data%id_sed_doc == 0 ) &
               print *,"          NOTE: no c_sed_excretion_target_variable - root DOM is lost"
         ENDIF
      ENDIF
      IF ( .NOT. simStaticBiomass ) THEN
         IF ( ANY(data%mpars(:)%R_mort > zero_) .OR. ANY(data%mpars(:)%R_mort_bg > zero_) &
              .OR. data%epi_model > 0 ) THEN
            IF ( data%f_mort_dom < one_ .AND. data%id_poc == 0 ) &
               print *,"          NOTE: no c_mortality_target_variable - POM from mortality is lost"
            IF ( data%f_mort_dom > zero_ .AND. data%id_doc == 0 ) &
               print *,"          NOTE: no c_excretion_target_variable - DOM from mortality is lost"
         ENDIF
         IF ( data%simFruiting .AND. data%id_fruit == 0 ) &
            print *,"          NOTE: no fruit_target_variable - released fruit is lost"
      ENDIF
   ENDIF
   IF ( ANY(data%mpars(:)%light_model == 10 .AND. (data%mpars(:)%R_mort_bg > zero_ .OR. &
                                                   data%mpars(:)%R_resp_bg > zero_)) ) &
      print *,"          NOTE: light_model 10 has static below-ground biomass; R_resp_bg, R_mort_bg are not used"
   IF ( ANY((data%mpars(:)%light_model == 2 .OR. data%mpars(:)%light_model == 3) .AND. &
            data%mpars(:)%R_resp_bg > zero_) ) THEN
      print *,"          NOTE: light_model 2/3 with R_resp_bg > 0 - roots respire separately, so E_comp"
      print *,"                should be the leaf-only compensation irradiance (not the whole plant)"
   ENDIF

   ! Register environmental dependenciess
   data%id_tem      = aed_locate_global('temperature')
   data%id_extc     = aed_locate_global('extc_coef')
   data%id_sal      = aed_locate_global('salinity')
   data%id_dz       = aed_locate_global('layer_ht')
   data%id_par      = aed_locate_global('par')
   data%id_sed_zone = aed_locate_sheet_global('sed_zone')
   data%id_atem     = aed_locate_sheet_global('air_temp')
   data%id_yearday  = aed_locate_sheet_global('yearday')
   data%id_I_0      = aed_locate_sheet_global('par_sf')
   data%id_depth    = aed_locate_sheet_global('col_depth')


   ! Light spectra pre-processing tasks
   ! calculate the integrated irradiance over the 300-800 nm wave lengths
   WLint=wavei(27:103)             ! 300-800 nm
   landaiint=landai(27:103)        ! corresponding irradiance at 300-800 nm
   ! total irradiance integrate over the 300-800 nm
   landaiSUM=0
   DO ww=1,SIZE(WLint)-1
     landaimid = (landaiint(ww)+landaiint(ww+1))/2
     !WLmid = (WLint(ww)+WLint(ww+1))/2 !NOT USED
     landaiSUM=landaiSUM+landaimid*(WLint(ww+1)-WLint(ww))
   ENDDO
   ! leaf light absorbance at selected wavelengths; from Baird et al (2016)
   Wli = (/ 300,  390,  500,  530,  640,  680,  705,  800 /)
   ALli= (/ 0.72, 0.72, 0.68, 0.38, 0.38, 0.60, 0.04, 0.0 /)
   ! interpolate into fine-resolution wave lengths for further processing
   ALLOCATE (ALlint(501))
   CALL interp_0d(SIZE(Wli),Wli,ALli,SIZE(WLint),WLint,ALlint) ! ALlint=interp1(Wli,ALli,WLint)
   ALl=0.5412  !SUM(ALlint)/SIZE(ALlint)  ! overall (mean) leaf absorbance, fixed to 0.5412 as SUM(ALlint)/SIZE(ALlint) gives 0.083
END SUBROUTINE aed_define_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
   INTEGER  :: mi
   AED_REAL :: matz

   AED_REAL :: depth, biomass, n_shoots
   AED_REAL :: minBiomass, maxBiomass, minDepth, maxDepth, scale, offset, a
   AED_REAL :: MAC_A, MAC_B, MAC_F
!
!-------------------------------------------------------------------------------
!BEGIN

   !---------------------------------------------------------------------------
   ! Check to ensure this cell/zone is colonisable
   matz = _STATE_VAR_S_(data%id_sed_zone)

   IF ( .NOT. mac_in_zone(data, matz) ) THEN
     _DIAG_VAR_S_(data%id_mac) = zero_
     _DIAG_VAR_S_(data%id_mac_ag) = zero_
     _DIAG_VAR_S_(data%id_mac_bg) = zero_
     _DIAG_VAR_S_(data%id_canopy_lai) = zero_
     _DIAG_VAR_S_(data%id_gpp) = zero_
     _DIAG_VAR_S_(data%id_root_o) = zero_
     _DIAG_VAR_S_(data%id_root_d) = zero_
     _DIAG_VAR_S_(data%id_canopy_sh_dens) = zero_

     DO mi=1,data%num_mphy
      _STATE_VAR_S_(data%id_mphya(mi))  =  zero_
      _STATE_VAR_S_(data%id_mphyb(mi))  =  zero_
      IF( data%mpars(mi)%fruit_model >0 ) THEN
       _STATE_VAR_S_(data%id_mphyf(mi))  =  zero_
      ENDIF
     ENDDO
     RETURN
   ENDIF

   IF ( ANY(data%mpars(:)%fruit_model >0) ) _DIAG_VAR_S_(data%id_mac_gt) = one_

   !---------------------------------------------------------------------------
   ! (Re)set the biomass of each group
   !  -1 : keep the values the host has set (registered m_initial etc, or the
   !       host's initial conditions / restart)
   !   0 : m_initial, m_initial_bg, m_initial_f from the dbase
   !   1 : leaf biomass from depth, a logistic in log10(g DW/m2) from m_initial
   !       (above m_initial_d0) to m_initial_minb (below m_initial_mind); roots
   !       and fruit as f_bg and f_seed of the leaves (Cockburn Sound)
   ! In 0 and 1, groups locked to another zone (zone_lock) are set to zero.
   IF ( data%mac_initial >= 0 ) THEN
      depth = _STATE_VAR_S_(data%id_depth)      ! water depth above the cell

      DO mi=1,data%num_mphy
         IF( data%mpars(mi)%zone_lock >0 .and. data%mpars(mi)%zone_lock .ne. INT(matz) ) THEN
            MAC_A = zero_ ; MAC_B = zero_ ; MAC_F = zero_
         ELSEIF ( data%mac_initial == 0 ) THEN
            MAC_A = data%mpars(mi)%m_initial
            MAC_B = data%mpars(mi)%m_initial_bg
            MAC_F = data%mpars(mi)%m_initial_f
         ELSE
            ! log10 g DW/m2; the minimum is floored so log10 stays finite
            maxBiomass = LOG10(data%mpars(mi)%m_initial/data%mpars(mi)%X_cdw)
            minBiomass = MAX(data%mpars(mi)%m_initial_minb, data%mpars(mi)%m0, 1e-3)
            minBiomass = MIN(LOG10(minBiomass/data%mpars(mi)%X_cdw), maxBiomass)
            minDepth   = data%mpars(mi)%m_initial_d0
            maxDepth   = data%mpars(mi)%m_initial_mind

            scale      = 12. / (maxDepth-minDepth)
            offset     = -6./scale - minDepth
            a          = (depth+offset)*scale
            biomass    = minBiomass + exp(-a)/(1.+exp(-a))*(maxBiomass-minBiomass)

            MAC_A = data%mpars(mi)%X_cdw * 10.**biomass      ! mmol C/m2
            MAC_B = data%mpars(mi)%f_bg * MAC_A
            MAC_F = data%mpars(mi)%f_seed * MAC_A
         ENDIF

         _STATE_VAR_S_(data%id_mphya(mi)) = MAC_A
         _STATE_VAR_S_(data%id_mphyb(mi)) = MAC_B
         IF( data%mpars(mi)%fruit_model >0 ) _STATE_VAR_S_(data%id_mphyf(mi)) = MAC_F
      ENDDO

      IF ( data%simEpiphytes ) _STATE_VAR_S_(data%id_epi) = data%epi_initial
   ENDIF

   !---------------------------------------------------------------------------
   ! Initial diagnostics, summed as in the time step
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_
   _DIAG_VAR_S_(data%id_canopy_lai) = zero_
   _DIAG_VAR_S_(data%id_canopy_sh_dens) = zero_
   _DIAG_VAR_S_(data%id_canopy_biovolume) = zero_
   _DIAG_VAR_S_(data%id_gpp) = zero_
   _DIAG_VAR_S_(data%id_root_o) = zero_
   _DIAG_VAR_S_(data%id_root_d) = 0.05

   DO mi=1,data%num_mphy
      MAC_A = _STATE_VAR_S_(data%id_mphya(mi))
      MAC_B = _STATE_VAR_S_(data%id_mphyb(mi))
      MAC_F = zero_
      IF( data%mpars(mi)%fruit_model >0 ) MAC_F = _STATE_VAR_S_(data%id_mphyf(mi))
      _DIAG_VAR_S_(data%id_mac)    = _DIAG_VAR_S_(data%id_mac) + MAC_A + MAC_B + MAC_F
      _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + MAC_A
      _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + MAC_B
      _DIAG_VAR_S_(data%id_canopy_lai) = _DIAG_VAR_S_(data%id_canopy_lai) + &
                    (one_ - exp(-data%mpars(mi)%k_omega * MAC_A))
   ENDDO
END SUBROUTINE aed_initialize_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_column_macrophyte(data,column,layer_map)
!-------------------------------------------------------------------------------
! STEP 1 (epi_model 2): works up from the bed through the layers the canopy
! occupies. Sets each layer's canopy properties (blockage, frontal area, ke_mac)
! and epiphyte net production (epi_npp); epi_net is the column total, and
! epi_col = 1 tells step 2c the result is ready. No fluxes are set here: the
! epiphyte state is updated in step 2c, the water exchange in step 3.
! The canopy geometry is from the current leaf C (mac_canopy), as in step 2b.
! layer_map(1) is the top layer, layer_map(SIZE(layer_map)) the bottom one.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_map(:)
!
!LOCALS
   INTEGER  :: layer_idx, layer, mi
   AED_REAL :: dz, matz, hc
   AED_REAL :: hgt_layer_bot, hgt_canopy, hgt_layer_top, layer_frac, canopy_frac
   AED_REAL :: epi, epi_l, epi_prod, epi_resp, epi_net, epi_max, leaf_area
   AED_REAL :: par, extc, temp, Io, canopy_par
   AED_REAL :: sum_net, sum_prod, sum_resp
   AED_REAL :: sh_dens, sh_diam, biovolume
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( data%epi_model /= 2 ) RETURN

   ! Check this column is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   IF ( .NOT. mac_in_zone(data, matz) ) RETURN

   ! Canopy only spans the water column for submerged plants
   mi = 1
   IF ( .NOT. ANY(data%mpars(:)%growth_form == SUBMERGED) ) RETURN

   ! Initialise all layers canopy properties to 0 (these are on the 3D grid)
   DO layer = SIZE(layer_map),1,-1
      layer_idx = layer_map(layer)
      _DIAG_VAR_(data%id_canopy_blockage) = zero_
      _DIAG_VAR_(data%id_canopy_frarea) = zero_
      _DIAG_VAR_(data%id_kemac) = zero_
      _DIAG_VAR_(data%id_epi_npp) = zero_
   END DO

   ! Canopy in this column, and epiphyte capacity
   CALL mac_canopy(data, column, layer_map(SIZE(layer_map)), matz, sh_dens, sh_diam, hgt_canopy, biovolume)
   leaf_area = 3.142*sh_diam*hgt_canopy*sh_dens
   hgt_canopy = MAX(hgt_canopy, 1e-3)
   epi_max = data%epi_max * leaf_area   ! mmol C/m2 benthos
   epi = _STATE_VAR_S_(data%id_epi)     ! mmol C/m2 benthos, spread evenly over the canopy height
   Io = _STATE_VAR_S_(data%id_I_0)

   ! Work up through the layers
   sum_net = zero_ ; sum_prod = zero_ ; sum_resp = zero_
   hgt_layer_top = zero_
   DO layer = SIZE(layer_map),1,-1
      layer_idx = layer_map(layer)
      dz   = _STATE_VAR_(data%id_dz)
      par  = _STATE_VAR_(data%id_par)
      extc = _STATE_VAR_(data%id_extc)
      temp = _STATE_VAR_(data%id_tem)
      hgt_layer_bot = hgt_layer_top
      hgt_layer_top = hgt_layer_top + dz

      IF ( hgt_canopy <= hgt_layer_top ) THEN
         ! canopy ends in this layer
         hc = hgt_canopy - hgt_layer_bot                       ! canopy height in this layer
         canopy_par = MAX( MIN( par * EXP(-extc*(hgt_layer_top-hgt_canopy)), Io), zero_)
         _DIAG_VAR_S_(data%id_d_par) = canopy_par               ! at the top of the canopy
      ELSE
         ! canopy fully spans through layer
         hc = dz
         canopy_par = par
      ENDIF
      layer_frac  = hc/dz             ! vertical fraction of the layer the canopy occupies
      canopy_frac = hc/hgt_canopy     ! fraction of the canopy in this layer

      ! Set layer's macrophyte biovolume/blockage, and frontal area (e.g. for layer specific drag)
      _DIAG_VAR_(data%id_canopy_blockage) = biovolume * canopy_frac
      _DIAG_VAR_(data%id_canopy_frarea) = sh_dens * sh_diam * (dz*layer_frac)
      _DIAG_VAR_(data%id_kemac) = &
           _DIAG_VAR_(data%id_canopy_blockage) * data%mpars(mi)%KeMAC ! or maybe use Aeff?

      ! Epiphytes on the leaves in this layer (mmol C/m2 benthos/s)
      epi_l = epi * canopy_frac
      epi_prod = mac_epi_prod(data, column, layer_idx, canopy_par, extc, Io, hc, temp, epi, epi_max)
      epi_resp = data%R_epir*(data%theta_epi_resp**(temp-20.))
      epi_net  = (epi_prod-epi_resp)*epi_l
      _DIAG_VAR_(data%id_epi_npp) = epi_net/dz * secs_per_day   ! per m3 of this layer
      sum_net  = sum_net  + epi_net
      sum_prod = sum_prod + epi_prod*epi_l
      sum_resp = sum_resp + epi_resp*epi_l

      IF ( hgt_canopy <= hgt_layer_top ) EXIT ! top of canopy was found so stop looping up
   END DO
   ! End vertical/column loop

   _DIAG_VAR_S_(data%id_epi_net) = sum_net * secs_per_day
   _DIAG_VAR_S_(data%id_epi_col) = one_
   IF (diag_level>1) _DIAG_VAR_S_(data%id_epig) = sum_prod * secs_per_day
   IF (diag_level>1) _DIAG_VAR_S_(data%id_epir) = sum_resp * secs_per_day

   ! UPDATE CANOPY POSTURE & DRAG (future)
   ! Based on cell velocity and canopy geometry, compute in-canopy velocity
!  CALL sav_canopy_velocity( U , canopy_height , dz , u1 , u2 , drag )
   ! Use in canopy velocity to update canopy posture, and drag
!  CALL sav_posture_model( u1 , canopy_height , drag)
END SUBROUTINE aed_calculate_column_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_macrophyte(data,column,layer_idx)
!-------------------------------------------------------------------------------
! STEP 3 (epi_model 2): the water exchange of the epiphytes in this layer, from
! epi_npp of step 1. Only when step 2c used the column result (epi_col = 2), so
! the water exchange always matches the change in the epiphyte state.
! Pelagic fluxes are per m3 of the layer.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( data%epi_model /= 2 .OR. .NOT. data%simMacFeedback ) RETURN
   IF ( .NOT. (_DIAG_VAR_S_(data%id_epi_col) > 1.5) ) RETURN   ! also if NaN

   CALL mac_epi_exchange(data, column, layer_idx, _DIAG_VAR_(data%id_epi_npp)/secs_per_day)
END SUBROUTINE aed_calculate_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_macrophyte(data,column,layer_idx)
!-------------------------------------------------------------------------------
! STEP 2: a. plants (production, respiration, losses, and exchange with the
! water and sediment), b. canopy geometry, c. epiphytes.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mi
   AED_REAL :: mphy, MAC_A, MAC_B, MAC_F       ! State
   AED_REAL :: mphy_flux, mphy_flux_a, mphy_flux_b
   AED_REAL :: fT, fI, fSal, fSal_g, fSal_r, fDO, fEpi
   AED_REAL :: extc, dz, par, Io, temp, salinity, day, matz

   AED_REAL :: canopy_hgt, canopy_extc, par_canopy, par_swi

   AED_REAL :: sh_diameter, sh_height, n_shoot, biovolume
   AED_REAL :: sh_dens, sh_diam, hgt_canopy

   AED_REAL :: primprod(data%num_mphy), respiration(data%num_mphy)

   AED_REAL :: landaimid, WLmid
   AED_REAL :: gpp, npp
   AED_REAL :: resp, respiration_A, respiration_B, R_resp_A, R_resp_B, theta_resp, k_resp
   AED_REAL :: mort, mort_A, mort_B, R_mort_A, R_mort_B
   AED_REAL :: f_tran, tau_tran, f_tran_fruit, tau_tran_fruit, f1, f2, f_seed, f_below
   AED_REAL :: Omega_MAC, sine_blade, E_comp,  R_growth, I_K, kA  , A_eff
   AED_REAL :: t_start_g, t_dur_g,t_max_g, t_max_r, t_start_r, t_dur_r
   AED_REAL :: light_int, lterm1, factor, kI, term1, term2, factor2, npp0, x, tmp1, tmp2, tmp12, tmp22
   AED_REAL :: MAC_F_release, r_release, f_release, trigger_fruit_growth
   INTEGER ::  ll
   LOGICAL ::  wet
   AED_REAL :: epi, epi_prod, epi_resp, epi_flux, epi_max, ratio, din, dip, fN, fP, K_epi, leaf_area
   AED_REAL :: c_w, c_s, nh4, no3, pNH4, sw_N, as_N, f_N, sw_P, as_P, f_P, f_nut, dnpp, KN
   AED_REAL :: wup_N, sup_N, wup_P, sup_P
   AED_REAL :: leaf_C, leaf_loss, epi_loss
   LOGICAL  :: use_col
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)

   IF ( .NOT. mac_in_zone(data, matz) ) return

   A_eff = 0. ; f_tran = 0. ; mphy_flux_a = 0. ; mphy_flux_b = 0. ; npp = 0.
   leaf_C = zero_ ; leaf_loss = zero_

   !--- BENTHIC CELL SETUP TASKS
   ! Retrieve current (local) environmental conditions in benthic cell
   salinity = _STATE_VAR_(data%id_sal)         ! salinity
   temp  = _STATE_VAR_(data%id_tem)            ! temperature
   extc = _STATE_VAR_(data%id_extc)            ! extinction
   par = _STATE_VAR_(data%id_par)              ! photosynthetically active radiation
   dz = _STATE_VAR_(data%id_dz)                ! cell thickness
   Io = _STATE_VAR_S_(data%id_I_0)             ! surface short wave radiation

   ! Adjust PAR to top of canopy
   !canopy_extc = MAX( 0.01,_DIAG_VAR_(data%id_kemac))       !0.1
   !canopy_hgt = _DIAG_VAR_S_(data%id_canopy_height)         !0.1 ;
   !par_canopy = _DIAG_VAR_S_(data%id_d_par) ! MAX( MIN( par * exp(-extc*( MAX(dz-canopy_hgt,zero_))),Io), zero_)
   canopy_extc = 0.1
   canopy_hgt = 0.1
   par_canopy = MAX( MIN( par * exp(-extc*( MAX(dz-canopy_hgt,zero_))),Io), zero_)
   par_swi = MAX( par_canopy * exp(-(extc+canopy_extc)*canopy_hgt), zero_)

   ! Initialise epiphyte density (for macrophyte fI limitation)
   epi = zero_
   IF ( data%simEpiphytes ) epi = _STATE_VAR_S_(data%id_epi)
   CALL mac_canopy(data, column, layer_idx, matz, sh_dens, sh_diam, hgt_canopy, biovolume)
   leaf_area = 3.142*sh_diam*hgt_canopy*sh_dens
   epi = epi / MAX(leaf_area, 1e-3)  ! ( mmolC epiphytes/m2 leaf = mmolC epiphytes /m2 benthos / m2leaf/m2 benthos )

   ! Reset cumulative biomass/canopy diagnostics (used to sum over all groups)
   _DIAG_VAR_S_(data%id_canopy_biovolume) = zero_
   _DIAG_VAR_S_(data%id_canopy_sh_dens) = zero_
   _DIAG_VAR_S_(data%id_canopy_sh_diam) = zero_
   _DIAG_VAR_S_(data%id_canopy_height) = zero_
   _DIAG_VAR_S_(data%id_canopy_lai) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_
   _DIAG_VAR_S_(data%id_mac_tr) = zero_
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_gpp) = zero_

   ! Initialise fruit growth control variables
   trigger_fruit_growth = zero_
   IF( data%simFruiting ) THEN
      _DIAG_VAR_S_(data%id_mac_fr) = zero_
      _DIAG_VAR_S_(data%id_mac_ft) = zero_
      trigger_fruit_growth = _DIAG_VAR_S_(data%id_mac_gt) ! =0 when start releasing
   ENDIF

   !--- LOOP THROUGH EACH GROUP/SPECIES
   ! Compute productivity, respiration, translocation, fruiting
   DO mi=1,data%num_mphy
      ! Check if this species/group is locked to grow in a specific zone
      IF( data%mpars(mi)%zone_lock >0 .and. data%mpars(mi)%zone_lock .ne. INT(matz) ) CYCLE

      ! Set local pars for species
      R_growth  = data%mpars(mi)%R_growth
      I_K       = data%mpars(mi)%I_K
      Omega_MAC = data%mpars(mi)%k_omega
      theta_resp= data%mpars(mi)%theta_resp
      f_below   = data%mpars(mi)%f_bg
      tau_tran  = data%mpars(mi)%tau_tran
      sine_blade= data%mpars(mi)%sine_blade
      E_comp    = data%mpars(mi)%E_comp
      kA        = data%mpars(mi)%kA
      R_mort_A  = data%mpars(mi)%R_mort
      R_mort_B  = data%mpars(mi)%R_mort_bg
      R_resp_A  = data%mpars(mi)%R_resp
      R_resp_B  = data%mpars(mi)%R_resp_bg
      K_epi     = data%mpars(mi)%K_epi

      !  Retrieve current (local cell) biomass for AG, BG pools
      MAC_A = _STATE_VAR_S_(data%id_mphya(mi))            ! leaf biomass
      MAC_B = _STATE_VAR_S_(data%id_mphyb(mi))            ! below-ground biomass
      MAC_F = zero_                                       ! fruit biomass (if fruiting)
      mort_A = zero_ ; mort_B = zero_ ; respiration_B = zero_ ; f_release = zero_

      !--- TEMPERATURE AND SALINITY, as for PHY and MAG (1 if temp_method, sal_method = 0)
      !# sal_method < 0 limits growth (fSal_g); > 0 enhances leaf respiration (fSal_r)
      fT = fTemp_function(data%mpars(mi)%temp_method,  &
                          data%mpars(mi)%T_max,        &
                          data%mpars(mi)%T_std,        &
                          data%mpars(mi)%theta_growth, &
                          data%mpars(mi)%aTn,          &
                          data%mpars(mi)%bTn,          &
                          data%mpars(mi)%kTn,temp)
      fSal = phyto_salinity(data%salpars,mi,salinity)
      fSal_g = one_ ; fSal_r = one_
      IF ( data%mpars(mi)%sal_method < 0 ) fSal_g = fSal
      IF ( data%mpars(mi)%sal_method > 0 ) fSal_r = fSal
      IF ( data%id_fT(mi) > 0 ) THEN
         _DIAG_VAR_S_(data%id_fT(mi))   = fT
         _DIAG_VAR_S_(data%id_fSal(mi)) = fSal
      ENDIF

      !--- GROWTH
      IF ( data%mpars(mi)%light_model == 1 ) THEN
         ! use Baird approach with total/bulk light approach

         ! gross primary production
         fEpi         = 1.0                               ! epiphyte limitation
         IF ( data%simEpiphytes ) fEpi = K_epi/(epi+K_epi)  ! epiphyte shading
         x            = (fEpi*par_canopy)/I_K             ! light limitation
         A_eff        = 1 - exp(-Omega_MAC*MAC_A)         ! effective area
         fI           = x/(1 + x) * (kA/(kA+A_eff))       ! light limitation+shelf shading
         primprod(mi) = R_growth*(min(fI,one_)) * fT * fSal_g * MAC_A ! GPP rate (mmol C/m2/s)

         ! respiration; note the definition of respiration in total light
         ! model is different to the one in the spectral light model. Here
         ! the respiration includes pure respiratory fraction, mortality and excretion
         respiration_A = R_resp_A * theta_resp**(temp-20.0) * fSal_r * MAC_A
         respiration_B = R_resp_B * theta_resp**(temp-20.0) * MAC_B

         respiration(mi) = respiration_A + respiration_B  ! Resp rate (mmol C/m2/s)

         ! mortality (leaf senescence, root die-off)
         mort_A = R_mort_A * theta_resp**(temp-20.0) * MAC_A
         mort_B = R_mort_B * theta_resp**(temp-20.0) * MAC_B

         ! translocation
         f_tran = (f_below-(MAC_B)/(MAC_A+MAC_B)) * (MAC_A+MAC_B) * tau_tran

         ! net flux rates (mmol C/m2/s)
         mphy_flux_a = (primprod(mi)  - respiration_A - mort_A - f_tran)
         mphy_flux_b = (- respiration_B - mort_B + f_tran)
         npp = primprod(mi)  - respiration_A
      ELSEIF ( data%mpars(mi)%light_model == 2 ) THEN
         ! use Baird approach with assumed (fixed) light spectra

         ! Integrate over 300-800nm to calculate photons being captured (Eq 8 in Baird et al 2016)
         fEpi         = 1.0                               ! epiphyte limitation
         IF ( data%simEpiphytes ) fEpi = 1.0              ! epiphyte shading
         light_int = zero_
         DO ll = 1,SIZE(WLint)-1
             lterm1    = 1. - exp(-ALlint(ll)*Omega_MAC*MAC_A*sine_blade)
             landaimid = (landaiint(ll)+landaiint(ll+1))/2.
             WLmid     = (WLint(ll)+WLint(ll+1))/2.
             light_int = light_int + landaimid*WLmid*(WLint(ll+1)-WLint(ll))*lterm1
         ENDDO
         light_int = light_int*(par_canopy*fEpi)/landaiSUM ! proportion of incoming par to clear-sky irradiance
         factor    = 1./(h*c*Av*1e9)                      ! conversion constant of photons from W/m2 to photon/m2/s
         kI        = factor*light_int                     ! rate of photon capture, mol photon/m2/s;
         primprod(mi) = kI  !?

         ! respiration
         ! compensation light (sum of respiration and mortality), mol photon/m2/day,
         term1 = E_comp*ALl*Omega_MAC*sine_blade

         term2 = 5500./550./1000. * R_mort_A * 86400.       ! mortality, converted to mol photon/m2/day;
         k_resp = (term1 - term2) * MAC_A /86400.         ! respiration rate in mol photon/m2/s, over 24 h
         IF (data%resp_mode == 0) k_resp = 2.*k_resp      ! daylight only, Eq(9) of Baird et al. 2016

         ! net production
         factor2 = 550./5500.*1000.                       ! factor to convert photon to carbon, mmol C/m2/s
         resp = k_resp * factor2 * fSal_r                 ! respiration rate in mmol C/m2/s
         respiration(mi) = resp
         gpp = kI*factor2 * fT * fSal_g                   ! gross production, mmol C/m2/s
         IF (data%resp_mode == 0) THEN
            npp0 = max(zero_,(gpp-resp))                  ! net production rate (none in the dark)
         ELSE
            npp0 = gpp - resp                             ! net production rate (< 0 in the dark)
         ENDIF
         npp = min(R_growth*fT*fSal_g*MAC_A,npp0)         ! capped at the maximum growth rate

         ! mortality
         mort_A = R_mort_A * theta_resp**(temp-20.0) * MAC_A
         mort_B = R_mort_B * theta_resp**(temp-20.0) * MAC_B

         ! below-ground respiration, separate from the leaf term above; with
         ! R_resp_bg > 0, E_comp should be the leaf-only compensation irradiance
         respiration_B = R_resp_B * theta_resp**(temp-20.0) * MAC_B
         respiration(mi) = resp + respiration_B

         ! translocation between AG/BG
         f_tran = (f_below - (MAC_B)/(MAC_A+MAC_B)) * (MAC_A+MAC_B) * tau_tran

         ! update AG and BG biomass rates, and effective projected area fraction, A_eff
         mphy_flux_a = (npp - mort_A - f_tran)
         mphy_flux_b = (- respiration_B - mort_B + f_tran)

         A_eff = 1 - exp(-Omega_MAC * MAC_A)
         _DIAG_VAR_S_(data%id_aeff) = A_eff

         !PRINT *,'LIGHT2',npp, term1, term2, kI, k_resp !, ALl, E_comp, Omega_MAC, sine_blade
      ELSEIF ( data%mpars(mi)%light_model == 3 ) THEN
         ! use Baird approach with AED spectrally-resolved light

         light_int = _DIAG_VAR_S_(data%id_parc)           ! light capture by plant leaves (returned from OASIM)
         factor    = 1./(h*c*Av*1e9)                      ! conversion constant of photons from W/m2 to photon/m2/s
         kI        = factor*light_int                     ! rate of photon capture, mol photon/m2/s;
         primprod(mi) = kI  !?

         ! respiration
         ! compensation light (sum of respiration and mortality), mol photon/m2/day, divided by 86400 to mol photon/m2/s
         term1 = E_comp*ALl*Omega_MAC*sine_blade/86400.

         term2 = 5500./550./1000. * R_mort_A              ! mortality, converted to mol photon/m2/s;
         k_resp = (term1 - term2) * MAC_A                 ! respiration rate in mol photon/m2/s, over 24 h
         IF (data%resp_mode == 0) k_resp = 2.*k_resp      ! daylight only, Eq(9) of Baird et al. 2016

         ! net production
         factor2 = 550./5500.*1000.                       ! factor to convert photon to carbon, mmol C/m2/s
         resp = k_resp * factor2 * fSal_r                 ! respiration rate in mmol C/m2/s
         respiration(mi) = resp
         gpp = kI*factor2 * fT * fSal_g                   ! gross production, mmol C/m2/s
         IF (data%resp_mode == 0) THEN
            npp0 = max(zero_,(gpp-resp))                  ! net production rate (none in the dark)
         ELSE
            npp0 = gpp - resp                             ! net production rate (< 0 in the dark)
         ENDIF
         npp = min(R_growth*fT*fSal_g*MAC_A,npp0)         ! capped at the maximum growth rate

         ! mortality
         mort_A = R_mort_A * theta_resp**(temp-20.0) * MAC_A
         mort_B = R_mort_B * theta_resp**(temp-20.0) * MAC_B

         ! below-ground respiration, separate from the leaf term above; with
         ! R_resp_bg > 0, E_comp should be the leaf-only compensation irradiance
         respiration_B = R_resp_B * theta_resp**(temp-20.0) * MAC_B
         respiration(mi) = resp + respiration_B

         ! translocation between AG/BG
         f_tran = (f_below - (MAC_B)/(MAC_A+MAC_B)) * (MAC_A+MAC_B) * tau_tran

         ! update AG and BG biomass rates, and effective projected area fraction, A_eff
         mphy_flux_a = (npp - mort_A - f_tran)
         mphy_flux_b = (- respiration_B - mort_B + f_tran)

         A_eff = 1 - exp(-Omega_MAC * MAC_A)
         _DIAG_VAR_S_(data%id_aeff) = A_eff
     ELSEIF ( data%mpars(mi)%light_model == 10 ) THEN
         ! use legacy AED approach with bulk light

         fI   = photosynthesis_irradiance(data%mpars(mi)%light_model, &
                                 I_K, data%mpars(mi)%I_S, par_canopy, extc, Io, dz)

         ! Primary productivity
         primprod(mi) = R_growth * (fI * fT * fSal_g) *  MAC_A

         ! Respiration and general metabolic loss
         respiration(mi) = bio_respiration(R_resp_A, theta_resp, temp)

         ! Salinity stress effect on respiration
         fDO  = 1.0 ! phyto_oxygen(data%mpars,mi,oxygen)

         respiration(mi) = respiration(mi) * (fSal_r * fDO) *  MAC_A

         ! mortality (leaf senescence; below-ground biomass is static here)
         mort_A = R_mort_A * theta_resp**(temp-20.0) * MAC_A

         npp = primprod(mi) - respiration(mi)
         mphy_flux_a = npp - mort_A
         mphy_flux_b = zero_
     ENDIF

     !--- NUTRIENT LIMITATION AND SOURCES OF UPTAKE
     !# For npp > 0 the N and P demand (npp * X_ncon, X_pcon) is met first by the
     !# water, up to its share water_nutrient_frac and less as the water runs low;
     !# the roots cover the rest, from a linked sediment pool or, if none is
     !# linked, from an assumed ample supply. So growth is only nutrient limited
     !# if water_nutrient_frac = 1, or the linked sediment pool is low. Growth is
     !# scaled by the scarcer of N and P. K_N, K_P = 0 switch this off.
     !# wup_*, sup_* : fractions of the actual uptake from the water, sediment
     wup_N = data%water_nutrient_frac ; sup_N = one_ - wup_N
     wup_P = data%water_nutrient_frac ; sup_P = one_ - wup_P
     pNH4 = 0.5
     IF( data%simMacFeedback .AND. npp > zero_ ) THEN
        KN = data%mpars(mi)%K_N
        nh4 = zero_ ; no3 = zero_
        IF(data%id_nh4>0) nh4 = MAX(_STATE_VAR_(data%id_nh4),zero_)
        IF(data%id_nox>0) no3 = MAX(_STATE_VAR_(data%id_nox),zero_)
        c_s = zero_
        IF(data%id_sed_nup>0) c_s = _STATE_VAR_S_(data%id_sed_nup) / (data%sc_nup*data%sed_layer_dz)
        CALL mac_nutrient_supply(data%water_nutrient_frac, (data%id_nh4>0 .OR. data%id_nox>0), &
                                 nh4+no3, data%id_sed_nup>0, c_s, KN, sw_N, as_N, f_N)

        c_w = zero_ ; c_s = zero_
        IF(data%id_po4>0) c_w = _STATE_VAR_(data%id_po4)
        IF(data%id_sed_pup>0) c_s = _STATE_VAR_S_(data%id_sed_pup) / (data%sc_pup*data%sed_layer_dz)
        CALL mac_nutrient_supply(data%water_nutrient_frac, data%id_po4>0, c_w, &
                                 data%id_sed_pup>0, c_s, data%mpars(mi)%K_P, sw_P, as_P, f_P)

        ! growth limited by the scarcer nutrient
        f_nut = MIN(f_N, f_P)
        dnpp = npp * (one_ - f_nut)
        npp = npp - dnpp
        mphy_flux_a = mphy_flux_a - dnpp

        IF( f_N > zero_ ) THEN
           wup_N = sw_N / f_N ; sup_N = (one_-sw_N) * as_N / f_N
        ENDIF
        IF( f_P > zero_ ) THEN
           wup_P = sw_P / f_P ; sup_P = (one_-sw_P) * as_P / f_P
        ENDIF

        ! ammonium preference (as for phytoplankton); concentration-weighted if K_N = 0
        IF( data%id_nh4>0 .AND. data%id_nox>0 ) THEN
           IF( KN > zero_ ) THEN
              pNH4 = zero_
              IF( nh4 > zero_ ) pNH4 = nh4*no3 / ((nh4+KN)*(no3+KN)) + nh4*KN / ((nh4+no3)*(no3+KN))
           ELSEIF( nh4+no3 > zero_ ) THEN
              pNH4 = nh4 / (nh4+no3)
           ENDIF
        ELSEIF( data%id_nh4>0 ) THEN
           pNH4 = one_
        ELSEIF( data%id_nox>0 ) THEN
           pNH4 = zero_
        ENDIF
     ENDIF

     !--- SET GROWTH/RESPIRATION RATES
     IF( .NOT. data%simStaticBiomass ) THEN
          ! Set bottom fluxes for the benthic pools (mmol C/m2/s)
          _FLUX_VAR_B_(data%id_mphya(mi)) = _FLUX_VAR_B_(data%id_mphya(mi)) + mphy_flux_a
          _FLUX_VAR_B_(data%id_mphyb(mi)) = _FLUX_VAR_B_(data%id_mphyb(mi)) + mphy_flux_b
          ! leaf C lost (senescence, or net decline under stress), for epiphyte sloughing
          leaf_C = leaf_C + MAC_A
          leaf_loss = leaf_loss + MAX(mort_A, -mphy_flux_a)
     ENDIF

     !--- FRUITING : seagrass fruit/seed development, & release trigger
     f_tran_fruit = zero_
     IF (data%mpars(mi)%fruit_model == 1) THEN
        f_seed         = data%mpars(mi)%f_seed
        tau_tran_fruit = data%mpars(mi)%tau_tran_fruit
        r_release      = data%mpars(mi)%r_release
        t_start_g      = data%mpars(mi)%t_start_g
        t_dur_g        = data%mpars(mi)%t_dur_g
        t_start_r      = data%mpars(mi)%t_start_r
        t_dur_r        = data%mpars(mi)%t_dur_r

        ! 1. get the day of the year, and prior fruit biomass
        MAC_F = max(zero_, _STATE_VAR_S_(data%id_mphyf(mi)))     ! fruits biomass
        day = FLOOR(_STATE_VAR_S_(data%id_yearday))  ! day number
        ! reset trigger for fruit to grow if the day<t_start_g

        IF (day<t_start_g) trigger_fruit_growth=one_

        ! 2. translocation from AG to fruit/seeds
        t_max_g = t_start_g+t_dur_g
        tmp1    = 12./t_dur_g*day+6.*(t_start_g+t_max_g)/(t_start_g-t_max_g)
        tmp2    = exp(-tmp1)
        f1      = 1./(1.+tmp2)

        f_tran_fruit = max(zero_,(f_seed-MAC_F/(MAC_F+MAC_A)) * (MAC_F+MAC_A)*tau_tran_fruit*f1)

        ! ...assume that if fruit ratio reach the 90% of maximum (mature) value,
        ! ...then stop fruit growth & start releasing at constant rate
        IF ( (MAC_F)/(MAC_F+MAC_A) > f_seed*0.9 ) THEN
           trigger_fruit_growth=zero_
        ENDIF

        ! 3. releasing
        t_max_r = t_start_r+t_dur_r
        tmp12   = 12. /t_dur_r*day + 6.*(t_start_r+t_max_r)/(t_start_r-t_max_r)
        tmp22   = exp(-tmp12)
        f2      = 1./(1.+tmp22)
        ! (with simStaticBiomass the stage is still tracked, but the pools are left unchanged)
        IF (trigger_fruit_growth>0.) THEN
           IF( .NOT. data%simStaticBiomass ) THEN
              _FLUX_VAR_B_(data%id_mphyf(mi)) = _FLUX_VAR_B_(data%id_mphyf(mi)) + f_tran_fruit
              _FLUX_VAR_B_(data%id_mphya(mi)) = _FLUX_VAR_B_(data%id_mphya(mi)) - f_tran_fruit
           ENDIF
        ELSE
            f_tran_fruit = zero_
            MAC_F_release = MAC_F*r_release
            f_release = MAC_F_release*f2
            ! released fruit leaves the plant (to fruit_target_variable, below)
            IF( .NOT. data%simStaticBiomass ) &
               _FLUX_VAR_B_(data%id_mphyf(mi)) = _FLUX_VAR_B_(data%id_mphyf(mi)) - f_release
        ENDIF
        _DIAG_VAR_S_(data%id_mac_gt) = trigger_fruit_growth
        _DIAG_VAR_S_(data%id_mac_fr) = _DIAG_VAR_S_(data%id_mac_fr) + MAC_F
        _DIAG_VAR_S_(data%id_mac_ft) = _DIAG_VAR_S_(data%id_mac_ft) - f_tran_fruit
     ENDIF


     !--- FEEDBACK TO WATER COLUMN (CARBON, OXYGEN, NUTRIENTS)
     !# Benthic fluxes are per m2 of bed (the API call converts them to volumetric 
     !# rates based on the cell's thickness. Net production npp is exchanged at 
     !# the species' X_ncon / X_pcon :
     !#  - npp > 0: photosynthesis, DIC -> plant, O2 released, N and P taken from
     !#    the water (wup_*, above; the rest from the sediment, below); reducing
     !#    nitrate to organic N releases a further 2 O2 per N;
     !#  - npp < 0: respiration, plant -> DIC, O2 used, N and P returned as DON
     !#    and DOP (water_excr_frac; the rest to the sediment, below).
     !# No feedback when the bottom cell is thinner than active_layer_limit (m),
     !# as the host spreads the flux over dz and very thin cells may blow up.
     IF( data%simMacFeedback .AND. dz>data%active_layer_limit ) THEN
        IF(data%id_oxy>0) _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + npp
        IF(data%id_dic>0) _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - npp
        IF(npp>0) THEN
          IF(data%id_nox>0) _FLUX_VAR_(data%id_nox) = &
                   _FLUX_VAR_(data%id_nox) - npp * data%mpars(mi)%X_ncon * wup_N * (one_-pNH4)
          IF(data%id_nh4>0) _FLUX_VAR_(data%id_nh4) = &
                   _FLUX_VAR_(data%id_nh4) - npp * data%mpars(mi)%X_ncon * wup_N * pNH4
          IF(data%id_po4>0) _FLUX_VAR_(data%id_po4) = &
                   _FLUX_VAR_(data%id_po4) - npp * data%mpars(mi)%X_pcon * wup_P
          ! O2 from nitrate reduction, water and sediment nitrate (released at the leaves)
          IF(data%id_oxy>0) THEN
            IF(data%id_nox>0) _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) &
                   + 2. * npp * data%mpars(mi)%X_ncon * wup_N * (one_-pNH4)
            IF(data%sed_nup_no3) _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) &
                   + 2. * npp * data%mpars(mi)%X_ncon * sup_N
          ENDIF
        ELSE
          IF(data%id_don>0) _FLUX_VAR_(data%id_don) = &
                   _FLUX_VAR_(data%id_don) - npp * data%mpars(mi)%X_ncon * data%water_excr_frac
          IF(data%id_dop>0) _FLUX_VAR_(data%id_dop) = &
                   _FLUX_VAR_(data%id_dop) - npp * data%mpars(mi)%X_pcon * data%water_excr_frac
        ENDIF
     ENDIF

     !--- PLANT LOSSES (MORTALITY, FRUIT RELEASE)
     !# Leaf senescence (mort_A) is all released to the water; root die-off
     !# (mort_B) and below-ground respiration (light_model 1) only at
     !# water_excr_frac, the rest going to the sediment (below). Of the mortality,
     !# f_mort_dom goes to DOM and the rest to POM. Released fruit goes to
     !# fruit_target_variable (C), with its N and P to the PON/POP links.
     !# The POM and fruit links may be water or sediment variables (eg to send
     !# leaf litter to the sediment). Unlinked pathways leave the system.
     IF( data%simMacFeedback .AND. .NOT. data%simStaticBiomass ) THEN
        ! mortality (as carbon released from the plant to the water)
        mort = mort_A + mort_B * data%water_excr_frac

        ! POM and fruit (no thin layer check for sediment links)
        wet = dz>data%active_layer_limit
        IF( wet .OR. data%poc_sheet ) CALL mac_add_flux(column, layer_idx, data%id_poc, data%poc_sheet, &
                 mort * (one_-data%f_mort_dom))
        IF( wet .OR. data%pon_sheet ) CALL mac_add_flux(column, layer_idx, data%id_pon, data%pon_sheet, &
                 (mort*(one_-data%f_mort_dom) + f_release) * data%mpars(mi)%X_ncon)
        IF( wet .OR. data%pop_sheet ) CALL mac_add_flux(column, layer_idx, data%id_pop, data%pop_sheet, &
                 (mort*(one_-data%f_mort_dom) + f_release) * data%mpars(mi)%X_pcon)
        IF( wet .OR. data%fruit_sheet ) CALL mac_add_flux(column, layer_idx, data%id_fruit, data%fruit_sheet, &
                 f_release)
     ENDIF
     IF( data%simMacFeedback .AND. .NOT. data%simStaticBiomass .AND. &
                                       dz>data%active_layer_limit ) THEN
        ! below-ground respiration: roots draw all their O2 from the water
        ! (lacunal transport); the DIC and DOM are split by water_excr_frac
        IF(data%id_oxy>0) _FLUX_VAR_(data%id_oxy) = &
                 _FLUX_VAR_(data%id_oxy) - respiration_B
        IF(data%id_dic>0) _FLUX_VAR_(data%id_dic) = &
                 _FLUX_VAR_(data%id_dic) + respiration_B * data%water_excr_frac
        IF(data%id_don>0) _FLUX_VAR_(data%id_don) = &
                 _FLUX_VAR_(data%id_don) + respiration_B * data%mpars(mi)%X_ncon * data%water_excr_frac
        IF(data%id_dop>0) _FLUX_VAR_(data%id_dop) = &
                 _FLUX_VAR_(data%id_dop) + respiration_B * data%mpars(mi)%X_pcon * data%water_excr_frac

        ! mortality, DOM
        IF(data%id_doc>0) _FLUX_VAR_(data%id_doc) = &
                 _FLUX_VAR_(data%id_doc) + mort * data%f_mort_dom
        IF(data%id_don>0) _FLUX_VAR_(data%id_don) = &
                 _FLUX_VAR_(data%id_don) + mort * data%f_mort_dom * data%mpars(mi)%X_ncon
        IF(data%id_dop>0) _FLUX_VAR_(data%id_dop) = &
                 _FLUX_VAR_(data%id_dop) + mort * data%f_mort_dom * data%mpars(mi)%X_pcon
     ENDIF

     !--- EXCHANGE WITH THE SEDIMENT
     !# The complement of the water exchange: the roots' share of N and P uptake
     !# (sup_*, above), and (1-water_excr_frac) of N and P release, root respiration and
     !# root die-off. Targets held per m3 (eg porewater NH4) get the areal flux
     !# divided by sed_layer_dz, areal (per m2) targets get it as is (sc_*,
     !# set from the target's units). Unlinked pathways leave the system.
     IF( data%simMacFeedback ) THEN
        IF(npp>0) THEN
          IF(data%id_sed_nup>0) _FLUX_VAR_B_(data%id_sed_nup) = _FLUX_VAR_B_(data%id_sed_nup) &
                   - npp * data%mpars(mi)%X_ncon * sup_N * data%sc_nup
          IF(data%id_sed_pup>0) _FLUX_VAR_B_(data%id_sed_pup) = _FLUX_VAR_B_(data%id_sed_pup) &
                   - npp * data%mpars(mi)%X_pcon * sup_P * data%sc_pup
        ELSE
          IF(data%id_sed_don>0) _FLUX_VAR_B_(data%id_sed_don) = _FLUX_VAR_B_(data%id_sed_don) &
                   - npp * data%mpars(mi)%X_ncon * (one_-data%water_excr_frac) * data%sc_don
          IF(data%id_sed_dop>0) _FLUX_VAR_B_(data%id_sed_dop) = _FLUX_VAR_B_(data%id_sed_dop) &
                   - npp * data%mpars(mi)%X_pcon * (one_-data%water_excr_frac) * data%sc_dop
        ENDIF

        IF( .NOT. data%simStaticBiomass ) THEN
          ! below-ground respiration, and root die-off (as carbon)
          resp = respiration_B * (one_-data%water_excr_frac)
          mort = mort_B * (one_-data%water_excr_frac)

          IF(data%id_sed_dic>0) _FLUX_VAR_B_(data%id_sed_dic) = _FLUX_VAR_B_(data%id_sed_dic) &
                   + resp * data%sc_dic
          IF(data%id_sed_doc>0) _FLUX_VAR_B_(data%id_sed_doc) = _FLUX_VAR_B_(data%id_sed_doc) &
                   + mort * data%f_mort_dom * data%sc_doc
          IF(data%id_sed_don>0) _FLUX_VAR_B_(data%id_sed_don) = _FLUX_VAR_B_(data%id_sed_don) &
                   + (resp + mort*data%f_mort_dom) * data%mpars(mi)%X_ncon * data%sc_don
          IF(data%id_sed_dop>0) _FLUX_VAR_B_(data%id_sed_dop) = _FLUX_VAR_B_(data%id_sed_dop) &
                   + (resp + mort*data%f_mort_dom) * data%mpars(mi)%X_pcon * data%sc_dop
          IF(data%id_sed_poc>0) _FLUX_VAR_B_(data%id_sed_poc) = _FLUX_VAR_B_(data%id_sed_poc) &
                   + mort * (one_-data%f_mort_dom) * data%sc_poc
          IF(data%id_sed_pon>0) _FLUX_VAR_B_(data%id_sed_pon) = _FLUX_VAR_B_(data%id_sed_pon) &
                   + mort * (one_-data%f_mort_dom) * data%mpars(mi)%X_ncon * data%sc_pon
          IF(data%id_sed_pop>0) _FLUX_VAR_B_(data%id_sed_pop) = _FLUX_VAR_B_(data%id_sed_pop) &
                   + mort * (one_-data%f_mort_dom) * data%mpars(mi)%X_pcon * data%sc_pop
        ENDIF
     ENDIF

     !--- CANOPY CHARACTERISTICS (step 2b; diagnostics only - see mac_canopy)
     CALL mac_shoots(data, mi, MAC_A, n_shoot, sh_height, sh_diameter, biovolume)

     !--- ADD GROUP INFO TO TOTAL MAC COMMUNITY / CANOPY
     ! Increment updates to the plant community level diagnistics for this group
     _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + MAC_A
     _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + MAC_B
     _DIAG_VAR_S_(data%id_mac_tr) = _DIAG_VAR_S_(data%id_mac_tr) - f_tran
     _DIAG_VAR_S_(data%id_mac)    = _DIAG_VAR_S_(data%id_mac) + (MAC_A+MAC_B+MAC_F)
     _DIAG_VAR_S_(data%id_gpp)    = _DIAG_VAR_S_(data%id_gpp) + npp*secs_per_day ! (mmol C/m2/day)

     ! Increment canopy characteristics
     _DIAG_VAR_S_(data%id_canopy_sh_dens) = _DIAG_VAR_S_(data%id_canopy_sh_dens) + n_shoot
     _DIAG_VAR_S_(data%id_canopy_sh_diam) = _DIAG_VAR_S_(data%id_canopy_sh_diam) + sh_diameter
     _DIAG_VAR_S_(data%id_canopy_biovolume) = _DIAG_VAR_S_(data%id_canopy_biovolume) + biovolume
     _DIAG_VAR_S_(data%id_canopy_height) = _DIAG_VAR_S_(data%id_canopy_height) + sh_height
     _DIAG_VAR_S_(data%id_canopy_lai) = _DIAG_VAR_S_(data%id_canopy_lai) + A_eff  ! TBC+       &
                     ! (one_ - exp(-data%mpars(mi)%k_omega * mphy*(one_-data%mpars(mi)%f_bg)))
   ENDDO ! Finish looping through groups

   !--- FINALISE CANOPY & COMMUNITY CALCULATIONS

   ! Finalise canopy averaging
   _DIAG_VAR_S_(data%id_canopy_sh_diam) = _DIAG_VAR_S_(data%id_canopy_sh_diam) / REAL(data%num_mphy)
   _DIAG_VAR_S_(data%id_canopy_height) = _DIAG_VAR_S_(data%id_canopy_height) / REAL(data%num_mphy)

   ! Approximate root depth and oxygen input (eg., for sediment bgc model)
   _DIAG_VAR_S_(data%id_root_o) = _DIAG_VAR_S_(data%id_gpp) * 0.2 ! data%bg_gpp_frac     ! mmolO2/m2/d
   _DIAG_VAR_S_(data%id_root_d) = 0.05 !MAX( MIN(_DIAG_VAR_S_(data%id_mac_ag) * data%coef_bm_hgt,0.25),0.01) !m

   ! Export additional diagnostic variables
   _DIAG_VAR_S_(data%id_d_par)= par_canopy

   !--- EPIPHYTES (step 2c)
   !# Epiphytes live on the leaves and exchange with the water only (see
   !# mac_epi_exchange). Net production comes from step 1 (epi_model 2, water
   !# exchange in step 3), or is computed here for the bottom layer (epi_model 1,
   !# or 2 if the host made no column call). Losses go to the plants' mortality
   !# links (POM, with f_mort_dom as DOM): a background rate R_epib, plus
   !# sloughing with the leaves - epiphytes are lost at the fractional rate the
   !# leaves are lost, by senescence or die-back under stress (no extra parameter).
   IF( data%epi_model > 0 ) THEN
      epi_prod = zero_ ; epi_resp = zero_
      epi = _STATE_VAR_S_(data%id_epi)  ! biomass of epiphytes in the benthos
      wet = dz>data%active_layer_limit

      use_col = .FALSE.
      IF ( data%epi_model == 2 ) use_col = _DIAG_VAR_S_(data%id_epi_col) > 0.5
      IF ( use_col ) THEN
         ! net production of step 1; step 3 applies the water exchange
         epi_flux = _DIAG_VAR_S_(data%id_epi_net) / secs_per_day
         _DIAG_VAR_S_(data%id_epi_col) = 2.
      ELSE
         ! bulk: the canopy as one layer, exchanging with the bottom water layer
         leaf_area = (3.142*(_DIAG_VAR_S_(data%id_canopy_sh_diam))* &
                             _DIAG_VAR_S_(data%id_canopy_height))*_DIAG_VAR_S_(data%id_canopy_sh_dens)
         epi_max = data%epi_max * leaf_area  ! ( mmolC epiphytes /m2 benthos = mmolC epiphytes/m2 leaf * m2leaf/m2 benthos )
         epi_prod = mac_epi_prod(data, column, layer_idx, par_canopy, extc, Io, dz, temp, epi, epi_max)
         epi_resp = data%R_epir*(data%theta_epi_resp**(temp-20.))
         epi_flux = (epi_prod-epi_resp)*epi
         IF (diag_level>1) _DIAG_VAR_S_(data%id_epig) = epi_prod * epi * secs_per_day
         IF (diag_level>1) _DIAG_VAR_S_(data%id_epir) = epi_resp * epi * secs_per_day
         IF ( data%simMacFeedback .AND. wet ) CALL mac_epi_exchange(data, column, layer_idx, epi_flux)
      ENDIF

      ! losses: background, and sloughing with the leaves
      epi_loss = data%R_epib*(data%theta_epi_resp**(temp-20.))
      IF ( leaf_C > zero_ ) epi_loss = epi_loss + leaf_loss/leaf_C
      epi_loss = epi_loss*epi

      _FLUX_VAR_B_(data%id_epi) = _FLUX_VAR_B_(data%id_epi) + epi_flux - epi_loss

      IF (diag_level>0) _DIAG_VAR_S_(data%id_epib) = epi
      IF (diag_level>1) _DIAG_VAR_S_(data%id_epil) = epi_loss * secs_per_day

      IF( data%simMacFeedback ) THEN
        ! losses, POM (the links may be water or sediment variables)
        IF( wet .OR. data%poc_sheet ) CALL mac_add_flux(column, layer_idx, data%id_poc, data%poc_sheet, &
                 epi_loss * (one_-data%f_mort_dom))
        IF( wet .OR. data%pon_sheet ) CALL mac_add_flux(column, layer_idx, data%id_pon, data%pon_sheet, &
                 epi_loss * (one_-data%f_mort_dom) * data%epi_Xnc)
        IF( wet .OR. data%pop_sheet ) CALL mac_add_flux(column, layer_idx, data%id_pop, data%pop_sheet, &
                 epi_loss * (one_-data%f_mort_dom) * data%epi_Xpc)
        ! losses, DOM
        IF( wet ) THEN
          IF (data%id_doc>0) _FLUX_VAR_(data%id_doc) = _FLUX_VAR_(data%id_doc) &
                   + epi_loss * data%f_mort_dom
          IF (data%id_don>0) _FLUX_VAR_(data%id_don) = _FLUX_VAR_(data%id_don) &
                   + epi_loss * data%f_mort_dom * data%epi_Xnc
          IF (data%id_dop>0) _FLUX_VAR_(data%id_dop) = _FLUX_VAR_(data%id_dop) &
                   + epi_loss * data%f_mort_dom * data%epi_Xpc
        ENDIF
      ENDIF
   ENDIF
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
   AED_REAL :: dz, mphy = 0., matz
   INTEGER  :: mi
   AED_REAL, PARAMETER :: max_extc = 2.0
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. mac_in_zone(data, matz) ) return

   ! Recall extinction value (for this cell) recorded during column loop calculations
   ! Accumulate, never assign: extinction holds earlier modules' contributions
   extinction = extinction + MIN( _DIAG_VAR_(data%id_kemac) , max_extc)
END SUBROUTINE aed_light_extinction_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_drag_macrophyte(data, column, layer_idx, drag)
!-------------------------------------------------------------------------------
! Get the effect of macrophyte biomass on benthic drag
!
!
!-------------------------------------------------------------------------------
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
  ! AED_REAL,INTENT(inout) :: blockage
!
!LOCALS
   INTEGER  :: mi, n_shoot
   AED_REAL :: sh_height, sh_diameter, m_biomass, blockage
   AED_REAL :: dz, matz, vel, depth, mphy
   AED_REAL :: Sx = one_, Sy = one_
   AED_REAL :: K_CD, Uc
!
!-------------------------------------------------------------------------------
!BEGIN
   drag = zero_

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. mac_in_zone(data, matz) ) return

   ! get the necessary environment information - cell thickness, water velocity
   vel  = _STATE_VAR_(data%id_vel)
   dz   = _STATE_VAR_(data%id_dz)
   depth= _STATE_VAR_S_(data%id_depth)

   IF( data%drag_model==1 ) THEN

      DO mi=1,data%num_mphy
         !# above ground density of macrophyte group i
         mphy = _STATE_VAR_S_(data%id_mphya(mi))

         !# additional drag due to biomass
         drag = drag + (data%mpars(mi)%K_CD * (mphy/dz) )
      ENDDO

   ELSEIF( data%drag_model==2 ) THEN
      ! set canopy average conditions for drag computation
      n_shoot = INT(_DIAG_VAR_S_(data%id_canopy_stem_dens))
      sh_diameter = _DIAG_VAR_S_(data%id_canopy_stem_diam)
      sh_height = _DIAG_VAR_S_(data%id_canopy_height)

      ! now do computation of velocity and drag within canopy
      Uc = vegetation_drag(n_shoot, data%mac_pattern, &
                          sh_diameter, Sx, Sy, vel, sh_height, K_CD)
      _DIAG_VAR_S_(data%id_canopy_velocity) = Uc

      !# additional drag due to canopy
      drag = K_CD
   ELSEIF( data%drag_model==3 ) THEN
      ! Work out leaf effective length etc and buoyancy
   ENDIF

  _DIAG_VAR_S_(data%id_canopy_drag) = drag
END SUBROUTINE aed_bio_drag_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION vegetation_drag(n,vp,d,Sx,Sy,Ub,hv,Cd) RESULT(Uc)
!-------------------------------------------------------------------------------
! Return the drag on the water flow given vegetation canopy characteristics
!-------------------------------------------------------------------------------
!#ARGUMENTS
   INTEGER, INTENT(in) :: n     !# Number of stems per unit area (1/m^2)
   INTEGER, INTENT(in) :: vp    !# Vegetation pattern/arrangement
   AED_REAL,INTENT(in) :: d     !# Mean stem diameter (m)
   AED_REAL,INTENT(in) :: Sx    !# Distance between elemetns in x-direction (streamwise direction) (m)
   AED_REAL,INTENT(in) :: Sy    !# Distance between elemetns in y-direction (lateral direction) (m)
   AED_REAL,INTENT(in) :: Ub    !# Bulk velocity that is the average cross-sectional velocity outside the canopy (m/s)
   AED_REAL,INTENT(in) :: hv    !# Mean height of vegetation elements
   AED_REAL,INTENT(inout) :: Cd !# Drag coefficient (computed here)
!
!#RETURN VAL
   AED_REAL :: Uc               !# Returned value of within-canopy velocity
!
!#LOCALS
   AED_REAL :: nu, Lambdap, Up
   AED_REAL :: nSqrt, Beta, Rec
!
!-------------------------------------------------------------------------------
!BEGIN
   !Environmental parameters; Units: metric;
   nu = 10**(-6)                ! Kinematic viscosity of water
   Lambdap = n * 3.1418  * d**2 / 4  ! Planar area (dimensionless)
   Up = Ub / (1-Lambdap);       ! Pore velocity (m/s)

   IF (vp == 1) THEN
     !(a) If linear or staggered arrangement, use a constricted velocity
     Beta  = Sx / Sy
     Uc    = (1-Lambdap) * Up/(1-sqrt(4*Lambdap/(3.1418 *Beta)))
   ELSE
     !(b) If random array of stems, use the following relationship for Uc :
     nSqrt = sqrt(REAL(n))
     Uc    = (1-Lambdap) * Up/(1-d*nSqrt)
   ENDIF

   Rec = Uc * d / nu            ! Stem Reynolds number
   Cd  = 1 + Rec**(-2/3)        ! Drag coefficient
END FUNCTION vegetation_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_0d(nsource,x,y,ntarget,targetx,targety)
!-------------------------------------------------------------------------------
! 0D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: nsource, ntarget
   AED_REAL,DIMENSION(nsource),INTENT(in)  :: x, y
   AED_REAL,DIMENSION(ntarget),INTENT(in)  :: targetx
   AED_REAL,DIMENSION(ntarget),INTENT(out) :: targety
!
!LOCALS
   INTEGER  :: i, j
   AED_REAL :: frac
!
!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO j = 1,ntarget
      DO while (i+1<nsource)
         IF (x(i+1)>=targetx(j)) exit
         i = i+1
      ENDDO
      frac = (targetx(j)-x(i))/(x(i+1)-x(i))
      targety(j) = y(i) + frac*(y(i+1)-y(i))
   ENDDO
END SUBROUTINE interp_0d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_macrophyte


#if 0


%% configuration for general parameters

% initial condition
MAC_A_INI = 25000;  % above-ground (AG, leaf without seeds) biomass, mmol C/m2
MAC_B_INI = 20000;  % below-ground (BG) biomass, mmol C/m2

% options
Light_model = 2;   % 1: total light model; 2: spectral-resolved model
Fruiting    = 1;   % 0: no fruiting; 1: include fruiting

% BGC parameters
R_growth = 0.4;    % maximum AG growth rate, /day
theta_resp=1.04;   % theta coefficient for respiration
T_standard=20;     % standard temperature for respiration
f_below=0.5;       % BG/(AG+BG) fraction
Omega_MAC = 1e-4;  %0.0011; % carbon-specific area of seagrass (mmol C/m2)^-1, converted from (240 mgC/m2)^-1
tau_tran = 0.033;  % AG/BG translocation rate, /day
sine_blade = 0.5;  % sine blade shape

% fruiting parameters
MAC_F_INI = 0;    % initial fruit biomass, mmol C/m2
tau_tran_fruit=0.01;  % translocation rate from leaves to fruit, /day
r_release=0.1;    % fruit releasing rate, /day (assume all fruit released in 10 days)

f_seed =0.4;       % seed biomass/(AG+seed biomass) fraction
t_start_g = 120;   % calenda day to grow, days of the year
t_dur_g = 30;      % duration to reach maximum growth, days

t_start_r = 210;   % calenda day to release fruits, days of the year
t_dur_r = 15;       % duration to reach maximum release, days




%% total light model parameters

I_K=80; % Half saturation constant for light limitation of growth
%I_S=120;  % Saturating light intensity for optimum photosynthesis
kA = 1;   % shelf shading effect

R_resp_A = 0.06;   % maximum AG respiration rate, /day
R_resp_B = 0.004;  % maximum BG respiration rate, /day

% fT = 1; % no temperature limiation
% fS = 1; % no salinity limiation
% fN = 1; % no sediment nutrient limiation

%% spectral-resolved parameters

WL_range = [300 800]; % wavelength range
WL_mean  = 550;       % mean wavelength, nm

h = 6.626e-34;  % Planck constant, Js;
c = 2.998e8;    % light speed, m/s;
Av= 6.02e23;    % Avagadro number, /mol, convert light capture to units of mol photon/m2/s;

E_comp = 15; %4.5*3.5; % compensation scalar PAR irradiance, mol photon/m2/d;
R_mort_A = 0.03;   % maximum AG mortality rate, /day
R_mort_B = 0.004;  % maximum BG mortality rate, /day


%% clear-sky wave lengths and irradiance at a particular wave-band W/m2/nm, from EMS model
% wave lengths
wavei = [140.00, 150.00, 160.00, 170.00, ...
    180.00, 190.00, 200.00, 205.00, 210.00, 215.00, 220.00, ...
    225.00, 230.00, 235.00, 240.00, 245.00, 250.00, 255.00, ...
    260.00, 265.00, 270.00, 275.00, 280.00, 285.00, 290.00, ...
    295.00, 300.00, 305.00, 310.00, 315.00, 320.00, 325.00, ...
    330.00, 335.00, 340.00, 345.00, 350.00, 355.00, 360.00, ...
    365.00, 370.00, 375.00, 380.00, 385.00, 390.00, 395.00, ...
    400.00, 405.00, 410.00, 415.00, 420.00, 425.00, 430.00, ...
    435.00, 440.00, 445.00, 450.00, 455.00, 460.00, 465.00, ...
    470.00, 475.00, 480.00, 485.00, 490.00, 495.00, 500.00, ...
    505.00, 510.00, 515.00, 520.00, 525.00, 530.00, 535.00, ...
    540.00, 545.00, 550.00, 555.00, 560.00, 565.00, 570.00, ...
    575.00, 580.00, 585.00, 590.00, 595.00, 600.00, 610.00, ...
    620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, ...
    690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, ...
    800.00, 850.00, 900.00, 950.00, 1000.00, 1100.00, 1200.00, ...
    1300.00, 1400.00, 1500.00, 1600.00, 1700.00, 1800.00, 1900.00, ...
    2000.00, 2100.00, 2200.00, 2300.00, 2400.00, 2500.00, 2600.00, ...
    2700.00, 2800.00, 2900.00, 3000.00, 3100.00, 3200.00, 3300.00, ...
    3400.00, 3500.00, 3600.00, 3700.00, 3800.00, 3900.00, 4000.00, ...
    4100.00, 4200.00, 4300.00, 4400.00, 4500.00, 4600.00, 4700.00, ...
    4800.00, 4900.00, 5000.00, 6000.00, 7000.00, 8000.00 ];

% clear-sky irradiance at a particular wave-band W/m2/nm, from EMS model
landai = [ 0.0000, 0.0001, 0.0000, 0.0004, ...
    0.0009,0.0017,0.0030,0.0050,0.0100,0.0180,0.0300,0.0420,0.0520, ...
    0.0540,0.0580,0.0640,0.0640,0.1000,0.1300,0.2000,0.2500,0.2200, ...
    0.2400,0.3400,0.5200,0.6300,0.6100,0.6700,0.7600,0.8200,0.8500, ...
    1.0200,1.1500,1.1100,1.1100,1.1700,1.1800,1.1600,1.1600,1.2900, ...
    1.3300,1.3200,1.2300,1.1500,1.1200,1.2000,1.5400,1.8800,1.9400, ...
    1.9200,1.9200,1.8900,1.7800,1.8200,2.0300,2.1500,2.2000,2.1900, ...
    2.1600,2.1500,2.1700,2.2000,2.1600,2.0300,1.9900,2.0400,1.9800, ...
    1.9700,1.9600,1.8900,1.8700,1.9200,1.9500,1.9700,1.9800,1.9800, ...
    1.9500,1.9200,1.9000,1.8900,1.8700,1.8700,1.8700,1.8500,1.8400, ...
    1.8300,1.8100,1.7700,1.7400,1.7000,1.6600,1.6200,1.5900,1.5500, ...
    1.5100,1.4800,1.4400,1.4100,1.3700,1.3400,1.3000,1.2700,1.1270, ...
    1.0030,0.8950,0.8030,0.7250,0.6060,0.5010,0.4060,0.3280,0.2670, ...
    0.2200,0.1820,0.1520,0.1274,0.1079,0.0917,0.0785,0.0676,0.0585, ...
    0.0509,0.0445,0.0390,0.0343,0.0303,0.0268,0.0230,0.0214,0.0191, ...
    0.0171,0.0153,0.0139,0.0125,0.0114,0.0103,0.0095,0.0087,0.0080, ...
    0.0073,0.0067,0.0061,0.0056,0.0051,0.0048,0.0044,0.0042,0.0021, ...
    0.0012,0.0006];

% calculate the integrated irradiance over the 300-800 nm wave lengths
WLint=wavei(27:103);              % 300-800 nm
landaiint=landai(27:103);         % corresponding irradiance at 300-800 nm

% total irradiance integrate over the 300-800 nm
landaiSUM=0;
for ww=1:length(WLint)-1
    landaimid = (landaiint(ww)+landaiint(ww+1))/2;
    WLmid = (WLint(ww)+WLint(ww+1))/2;
    landaiSUM=landaiSUM+landaimid*(WLint(ww+1)-WLint(ww));
end

%% leaf light absorbance at wave lengths, from Baird et al., 2016
Wli =[300  390  500  530  640  680  705  800];
ALli=[0.72 0.72 0.68 0.38 0.38 0.60 0.04 0.0];

% interpolate into fine-resolution wave lengths for further processing
ALlint=interp1(Wli,ALli,WLint);
ALl=mean(ALlint);  % overall leaf absorbance





%% This is a Matlab script for the seagrass dynamic model, only consider
%  the light limitation, e.g. no temperature, salinity, and sediment
%  nutrient limitation is included.
%
%  configurations: in Seagrass_model_config.m;
%  forcing data: surface PAR data at Cockburn Sound centre, exported
%                CSIEM model, then bottom PAR is calculated based on
%                pre-set water depth and extinction coefficient

clear; close all;

% load in configuration file
run('./Seagrass_model_config.m');

% load in light data
load('.\environmental_condition\light_CScenter.mat');
PAR_time = output_site.PAR_time;   % model times
timestep = 4/24;                   % timestpes
extc = 0.1;                        % extinction coefficient, /m
d = 5;                             % depth
par=output_site.PAR_surf.*exp(-extc*d);   % calculated bottom PAR


% allocate biomass for AG, BG, and Fruits
MAC_A = zeros(size(PAR_time));  % leaf biomass
MAC_B = zeros(size(PAR_time));  % below-ground biomass
MAC_F = zeros(size(PAR_time));  % fruits biomass

% allocate parameters for spectral-resoved model
npp = zeros(size(PAR_time));    % net production
resp = zeros(size(PAR_time));   % respiration of AG due to compensation light
mort_A = zeros(size(PAR_time)); % mortality of leaf
mort_B = zeros(size(PAR_time)); % mortality of roots
f_tran = zeros(size(PAR_time)); % translocation rates between AG/BG
A_eff = zeros(size(PAR_time));  % effective area

% allocate parameters for total light model
gpp = zeros(size(PAR_time));           % gross production
respiration_A = zeros(size(PAR_time)); % respiraton of leaf
respiration_B = zeros(size(PAR_time)); % respiration of roots

% allocate parameters for fruiting processes
f_tran_fruit = zeros(size(PAR_time));  % translocation rate of leaf to fruits
f1 = zeros(size(PAR_time));            % f1 function for growth
f2 = zeros(size(PAR_time));            % f2 function for release

% allocate parameters for spectral-resoved model

% initialize the biomass for AG, BG, and Fruits
MAC_A(1)=MAC_A_INI;
MAC_B(1)=MAC_B_INI;
MAC_F(1)=MAC_F_INI;

% a control on fruit growth, becomes 0 when start releasing
trigger_fruit_growth = 1;

% loop through the timesteps
for dd=2:length(PAR_time)

    % growth
    if Light_model == 2 % spectral-resolved model

        % integrated over 300-800nm to calculate photons being captured
        % Eq(8) of Baird et al. 2016
        light_int=0;
        for ll=1:length(WLint)-1
            lterm1 = 1 - exp(-ALlint(ll)*Omega_MAC.*MAC_A(dd-1)*sine_blade);
            landaimid = (landaiint(ll)+landaiint(ll+1))/2;
            WLmid = (WLint(ll)+WLint(ll+1))/2;

            light_int=light_int+landaimid*WLmid*(WLint(ll+1)-WLint(ll))*lterm1;
        end
        light_int=light_int*par(dd)/landaiSUM; %proportion of incoming irradiation par(dd) to clear-sky irradiance

        factor    = 1/(h*c*Av*1e9);    % conversion constant of photons from W/m2 to photon/m2/s
        kI = factor*light_int;         % rate of photon capture, mol photon/m2/s;

        % respiration

        term1 = E_comp*ALl*Omega_MAC.*sine_blade; % compensation light
        term2 = 5500/550/1000*R_mort_A;           % respiration, converted to mol photon/m2/s;
        k_resp = 2*(term1 - term2)*MAC_A(dd-1);   % respiration rate in photon, Eq(9) of Baird et al. 2016


        % net production
        factor2 = 550/5500*1000*86400; % factor to converting photon to carbon, mmol C/m2/day

        resp(dd) = k_resp*factor2/86400;        % respiration rate in mmol C/m2/day
        npp0 = max(0,(kI*factor2-resp(dd)));    % net production rate
        npp(dd) = min(R_growth*MAC_A(dd-1),npp0); % cross-check of NPP npp0; %

        % mortality
        mort_A(dd) = MAC_A(dd-1)*R_mort_A * theta_resp^(T_standard-20.0);
        mort_B(dd) = MAC_B(dd-1)*R_mort_B * theta_resp^(T_standard-20.0);

        % translocation between AG/BG
        f_tran(dd)=(f_below - (MAC_B(dd-1))/(MAC_A(dd-1)+MAC_B(dd-1)))*(MAC_A(dd-1)+MAC_B(dd-1))*tau_tran;

        % update AG and BG biomass, and effective projected area fraction
        % A_eff
        MAC_A(dd)= MAC_A(dd-1) + (npp(dd) - mort_A(dd) - f_tran(dd))*timestep;
        MAC_B(dd)= MAC_B(dd-1) + (- mort_B(dd) + f_tran(dd))*timestep;
        A_eff(dd) = 1 - exp(-Omega_MAC.*MAC_A(dd));

    elseif Light_model == 1 % total light model
        A_eff(dd) = 1 - exp(-Omega_MAC.*MAC_A(dd-1)); % effective area
        x = par(dd)/I_K;                             % light limitation
        fI = x ./ (1 + x)*(kA/(kA+A_eff(dd))); % light limitation+shelf shading
        gpp(dd) = MAC_A(dd-1)*R_growth*(min(fI,1)); % GPP rate per day

        % respiration; note the definition of respiration in total light
        % model is different to the one in the spectral light model. Here
        % the respiration includes pure respiratory fraction, mortality and excretion
        respiration_A(dd) = MAC_A(dd-1)*R_resp_A * theta_resp^(T_standard-20.0);
        respiration_B(dd) = MAC_B(dd-1)*R_resp_B * theta_resp^(T_standard-20.0);

        f_tran(dd)=(f_below - (MAC_B(dd-1))/(MAC_A(dd-1)+MAC_B(dd-1)))*(MAC_A(dd-1)+MAC_B(dd-1))*tau_tran;

        MAC_A(dd)= MAC_A(dd-1) + (gpp(dd) - respiration_A(dd) - f_tran(dd))*timestep;
        MAC_B(dd)= MAC_B(dd-1) + (- respiration_B(dd) + f_tran(dd))*timestep;

    else

        error('Light model option can be recognized');
    end

    % seagrass fruiting growth and release
    if Fruiting == 1

        % 1. get the day of the year
        t_vec = datevec(PAR_time(1));
        t = PAR_time(dd)-datenum(t_vec(1),1,1);

        % reset trigger for fruit to grow if t<t_start_g
        if t<t_start_g
            trigger_fruit_growth=1;
        end

        % 1. translocation from AG to seeds
        t_max_g=t_start_g+t_dur_g;
        tmp1=12/t_dur_g*t+6*(t_start_g+t_max_g)/(t_start_g-t_max_g);
        tmp2=exp(-tmp1);
        f1(dd)=1/(1+tmp2);

        f_tran_fruit(dd)=(f_seed-MAC_F(dd-1)/(MAC_F(dd-1)+MAC_A(dd-1)))*...
            (MAC_F(dd-1)+MAC_A(dd-1))*tau_tran_fruit*f1(dd);

        % assume if fruit ratio reach the 90% of maximum value, then
        % stop fruit growing and start releasing at constant rate
        if MAC_F(dd-1)/(MAC_F(dd-1)+MAC_A(dd-1))>f_seed*0.9
            trigger_fruit_growth=0;
            MAC_F_release=MAC_F(dd-1)*r_release;
        end

        % 2. releasing

        t_max_r=t_start_r+t_dur_r;
        tmp12=12/t_dur_r*t+6*(t_start_r+t_max_r)/(t_start_r-t_max_r);
        tmp22=exp(-tmp12);
        f2(dd)=1/(1+tmp22);
        if trigger_fruit_growth
            MAC_F(dd)=MAC_F(dd-1)+f_tran_fruit(dd);
            MAC_A(dd)=MAC_A(dd)-f_tran_fruit(dd);
        else
            f_tran_fruit(dd)=0;
            f_release=MAC_F_release*f2(dd);
            MAC_F(dd)=max(0,MAC_F(dd-1)-f_release);
            MAC_A(dd)=MAC_A(dd);
        end
    end



end

#endif
