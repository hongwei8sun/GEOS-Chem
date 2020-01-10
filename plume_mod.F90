!------------------------------------------------------------------------------
!         GEOS-Chem Global Chemical Stratospherical Transport Model                  !
!-----------------------------------------------------------------------------!

MODULE Plume_Mod

  USE precision_mod
  USE ERROR_MOD
  USE ErrCode_Mod
  USE PhysConstants,   ONLY : PI, Re, g0
  USE CMN_SIZE_Mod,    ONLY : IIPAR, JJPAR, LLPAR

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: plume_init
  PUBLIC :: plume_run
  PUBLIC :: plume_write_std
  PUBLIC :: plume_cleanup

  PUBLIC :: Plume_I, Plume_J, Plume_L ! (n_boxes_max)
  PUBLIC :: n_boxes_max, n_rings_max
  PUBLIC :: Plume_Species ! (n_boxes_max, n_rings_max, State_Chm%nSpecies)


  integer, allocatable  :: Plume_I(:), Plume_J(:), Plume_L(:) ! (n_boxes_max)
  real(fp), allocatable :: Plume_Species(:,:,:)

  integer, parameter    :: n_boxes_max = 3  ! 50*40*40
  integer               :: tt

  ! The center of the cylinder
  real(fp), allocatable :: box_lon(:)
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)

  ! Consider the box as a cylinder

  ! box_radius(n_boxes_max,N_rings)
  real(fp), allocatable :: box_radiusA(:,:) ! vertical radius
  real(fp), allocatable :: box_radiusB(:,:) ! horizontal radius of the the cross-section


  ! Theta is the clockwise angle between z-axis (P) and vertical radiusA
  real(fp), allocatable :: box_theta(:)    ! 0 ~ 180 degree
  real(fp), allocatable :: box_length(:)

  ! D_radius should only be used at the beginning!
  real(fp), parameter   :: Init_radius = 100.0e+0_fp ! [m], the width of each ring
  real(fp), parameter   :: D_radius    = 100.0e+0_fp ! [m], the width of each ring
  integer, parameter    :: n_rings_max = 500 ! Degine the number of rings in one box

  ! medical concentration of each ring
  real(fp), allocatable :: box_concnt(:,:,:) ! [kg/m3],(n_boxes_max,N_rings,N_species)
  real(fp), allocatable :: box_concnt_K(:)     ! [kg/m3], box_concnt_K(N_rings_max)
  real(fp), allocatable :: RK(:,:), AA_env(:)  ! for Runge-Kutta method
  real(fp), allocatable :: eddy_h(:)  ! 
  real(fp), allocatable :: eddy_A(:)  ! 
  real(fp), allocatable :: eddy_B(:)  ! 
  real(fp), allocatable :: kA(:), kB(:)  ! 

  real(fp), allocatable :: env_amount(:)
  real(fp), allocatable :: backgrd_concnt(:,:)

CONTAINS


!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  SUBROUTINE plume_init( am_I_root, Input_Opt, State_Chm, State_Met, RC )

    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Met_Mod, ONLY : MetState
    USE State_Chm_Mod, ONLY : ChmState
    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE
    USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON

    USE UnitConv_Mod,    ONLY : Convert_Spc_Units

    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    integer                       :: i_box, i_ring

    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt

    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    integer                       :: N_species
    integer                       :: i_species

    integer                       :: i_lon, i_lat, i_lev
    real(fp)                      :: curr_lon, curr_lat, curr_pressure

    real(fp)                      :: X_edge2, Y_edge2
    real(fp), pointer             :: X_edge(:), Y_edge(:), P_edge(:)

    real(fp)                      :: Dx, Dy

    CHARACTER(LEN=255)     :: FILENAME2
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg

    RC        =  GC_SUCCESS
    ErrMsg    =  ''

    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Plume Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'

    allocate(Plume_I(n_boxes_max))
    allocate(Plume_J(n_boxes_max))
    allocate(Plume_L(n_boxes_max))

    allocate(Plume_Species(n_boxes_max, n_rings_max, State_Chm%nSpecies))


    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))

    allocate(box_radiusA(n_boxes_max,n_rings_max))
    allocate(box_radiusB(n_boxes_max,n_rings_max))

    allocate(box_theta(n_boxes_max))
    allocate(box_length(n_boxes_max))

    N_species = State_Chm%nSpecies

    allocate(box_concnt(n_boxes_max,n_rings_max,N_species))
    allocate(box_concnt_K(n_rings_max))
    allocate(RK(4,n_rings_max))
    allocate(AA_env(4))

    allocate(eddy_h(n_rings_max))
    allocate(eddy_A(n_rings_max))
    allocate(eddy_B(n_rings_max))
    allocate(kA(n_rings_max))
    allocate(kB(n_rings_max))

    allocate(env_amount(n_boxes_max))
    allocate(backgrd_concnt(n_boxes_max,State_Chm%nSpecies))


    box_lon    = (/5.0e+0_fp,  5.1e+0_fp,  5.2e+0_fp/)
    box_lat    = (/4.0e+0_fp, 4.1e+0_fp, 4.2e+0_fp/)
    box_lev    = (/20.0e+0_fp, 20.0e+0_fp, 20.0e+0_fp/)      ! hPa

!    box_radiusA(:,1)  = (/5.0e+0_fp,  10.0e+0_fp,  10.0e+0_fp/)     
     ! the value of the innest ring for every box
!    box_radiusB(:,1)  = (/5.0e+0_fp,  10.0e+0_fp,  10.0e+0_fp/)     ! m


    ! Set the initial value of max/min radius for each ring
    do i_ring=1,n_rings_max
      box_radiusA(:,i_ring) = i_ring * D_radius
      box_radiusB(:,i_ring) = i_ring * D_radius
    enddo

    box_theta         = (/0.0e+0_fp,    0.0e+0_fp,    0.0e+0_fp/)             ! degree
    box_length        = (/1000.0e+0_fp, 1000.0e+0_fp, 1000.0e+0_fp/)    ! m



    ! put the Eulerian grid concentration into plume at the beginning
    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!

    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)


    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'plume_mod.F90')
       RETURN
    ENDIF


    do i_box = 1,n_boxes_max

       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)        ! hPa

       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(curr_pressure,P_edge)

       do i_species = 1,N_species
       do i_ring = 1,N_species
          box_concnt(i_box,i_ring,i_species) = State_Chm%Species(i_lon,i_lat,i_lev,i_species) 
          ! [kg/kg] ???
       enddo ! do i_ring = 1,N_species
       enddo ! do i_species = 1,N_species

    enddo ! i_box = 1,n_boxes_max


    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'plume_mod.F90' )
       RETURN
    ENDIF


    do i_ring=1,n_rings_max
!    do i_species = 1,N_species
!      box_concnt(:,i_ring,i_species) = (/0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp/)   ! [kg/m3]
!    enddo  
 
       box_concnt_K(i_ring) = 0.0e+0_fp

    enddo


        ! Sigma = 1000m, Dr = 100m
       box_concnt(1,1:100,1) = (/ 49.93753904622905, 49.44065223056165, 48.46166172381721, 47.0294031682171, 45.1853538936598, 42.98163818012711, 40.478582433394344, 37.741980099450366, 34.84023877480174, 31.84158071858716, 28.811453683589995, 25.810283697274816, 22.891668088580712, 20.10106915473274, 17.475030009987826, 15.040897717863777, 12.817007570753681, 10.813258341494365, 9.0319925809447, 7.469088762520901, 6.11517266734503, 4.95686262553737, 3.977975435911384, 3.1606351537644426, 2.48624367061748, 1.9362885175832174, 1.4929783320557766, 1.1397090441806172, 0.8613735655817554, 0.6445344572000735, 0.4774828697510154, 0.3502083574681211, 0.25430346155063505, 0.1828247940026448, 0.13012926263627014, 0.09170054150365799, 0.06397729689499619, 0.0441913153467525, 0.030220735187977774, 0.020461156511408944, 0.01371552337467006, 0.009102310518350318, 0.0059806441790512185, 0.003890463430776882, 0.002505600144312384, 0.001597641188721152, 0.0010085647538279085, 0.0006303552588524261, 0.000390053365115542, 0.00023895698661023193, 0.00014493473873164686, 8.703266982466485e-05, 5.1742710555468766e-05, 3.0456016577577623e-05, 1.77481904348566e-05, 1.0239815105987296e-05, 5.849075200675953e-06, 3.3078008188488504e-06, 1.8520322753971736e-06, 1.0266320224615827e-06, 5.634275999136617e-07, 3.061388760736478e-07, 1.6468570551530405e-07, 8.771025650048737e-08, 4.624895376357952e-08, 2.414407099189213e-08, 1.2478895312221815e-08, 6.3855577725641675e-09, 3.2350320733131727e-09, 1.6226146706350309e-09, 8.057665991536951e-10, 3.9615047298176933e-10, 1.928271364234862e-10, 9.292513057582789e-11, 4.433587289933671e-11, 2.094278344422346e-11, 9.79423670500225e-12, 4.5348595885495465e-12, 2.078806895816242e-12, 9.434556900632728e-13, 4.2392196790663004e-13, 1.8858510687236676e-13, 8.305885063498716e-14, 3.6217752542073234e-14, 1.5635584907363262e-14, 6.682882459159197e-15, 2.8279424573025658e-15, 1.184770873670088e-15, 4.914227126564088e-16, 2.0180556398551702e-16, 8.204801948335148e-17, 3.302631537603281e-17, 1.316161553421132e-17, 5.192964243507581e-18, 2.028516381535571e-18, 7.845105014892859e-19, 3.003834899930149e-19, 1.1387028598265411e-19, 4.273678215867316e-20, 1.5879992099449743e-20 /)

    env_amount = (/0.0e+0_fp, 0.0e+0_fp, 0.0e+0_fp/)


    ! Create output file
    FILENAME2   = 'Plume_theta_max_min_radius.txt'
    tt = 0

    OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!    Do i_ring = 1, n_rings_max
!       WRITE(262,'(I0.4,3(x,E16.5E4))') i_ring, box_theta(1), box_radiusA(1,i_ring), box_radiusB(1,i_ring)
!    End Do
     
     write(262,*)box_concnt(1,:,1)


  END SUBROUTINE plume_init


!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  SUBROUTINE plume_run(am_I_Root, State_Chm, State_Met, Input_Opt, RC)

    USE Input_Opt_Mod,   ONLY : OptInput

    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Met_Mod,   ONLY : MetState

    USE TIME_MOD,        ONLY : GET_TS_DYN

    USE GC_GRID_MOD,     ONLY : XEDGE, YEDGE, XMID, YMID                 
    USE CMN_SIZE_Mod,    ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON 
    ! DLAT( IIPAR, JJPAR, LLPAR ), DLON( IIPAR, JJPAR, LLPAR )
    ! XEDGE  ( IM+1, JM,   L ), YEDGE  ( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

    USE UnitConv_Mod,    ONLY : Convert_Spc_Units

    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt

    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    REAL :: Dt          ! = 600.0e+0_fp          
    REAL :: Dt2               

    integer  :: i_box, i_ring, i_species

    integer  :: i_lon, i_lat, i_lev
    
    integer  :: t1s, Ki

    real(fp) :: curr_lon, curr_lat, curr_pressure
    real(fp) :: next_lon, next_lat, next_pressure

    real(fp) :: X_edge2, Y_edge2

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    real(fp), pointer :: Ptemp(:,:,:)

    real(fp), pointer :: P_BXHEIGHT(:,:,:)

    real(fp)          :: curr_u, curr_v, curr_omeg

    real(fp)          :: curr_Ptemp, Ptemp_shear  ! potential temperature
    real(fp)          :: U_shear, V_shear, UV_shear  

    real(fp), pointer :: P_edge(:)
    real(fp), pointer :: P_mid(:)

    real(fp)          :: Dx, Dy
    real(fp)          :: wind_s_shear
    real(fp)          :: theta_previous       
    real(fp)          :: dbox_theta, dbox_radiusA, dbox_radiusB
    ! angle between travel direction and lon (-90 ~ 90 deg)
    real(fp)          :: box_alpha        

    real(fp), pointer :: X_edge(:)    
    real(fp), pointer :: Y_edge(:)       
    real(fp), pointer :: X_mid(:)    
    real(fp), pointer :: Y_mid(:)       

    real(fp), pointer :: P_I, R_e     

    real(fp)          :: AA(n_rings_max), BB(n_rings_max), DD(n_rings_max) 
    ! used to calculate concentration distribution

    real(fp)          :: eddy_v
    real(fp)          :: Cv, Ch, Omega_N

    real(fp)          :: grid_volumn, exchange_amount

    !real(fp) :: RK(4,n_rings_max)
    !real(fp) :: box_concnt_K(n_rings_max)

    CHARACTER(LEN=255)     :: FILENAME2
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg

    RC        =  GC_SUCCESS
    ErrMsg    =  ''

    FILENAME2   = 'Plume_theta_max_min_radius.txt'


    Dt = GET_TS_DYN()

    u => State_Met%U   ! figure out state_met%U is based on lat/lon or modelgrid(i,j)
    v => State_Met%V   ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    Ptemp => State_Met%THETA  ! Updraft velocity [Pa/s]

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]
    P_mid => State_Met%PMID(1,1,:)  ! Pressure (w/r/t moist air) at level centers (hPa)

    P_BXHEIGHT => State_Met%BXHEIGHT  ![IIPAR,JJPAR,KKPAR]

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!
    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)  
    X_mid => XMID(:,1,1)   ! IIPAR+1
    Y_mid => YMID(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)

!    N_species = State_Chm%nSpecies

    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90')
       RETURN
    ENDIF

    !=====================================================================
    ! Run Plume advection HERE
    !=====================================================================
    do i_box = 1,n_boxes_max

       ! make sure the location is not out of range
       do while (box_lat(i_box) > Y_edge(JJPAR+1))
          box_lat(i_box) = Y_edge(JJPAR+1) - ( box_lat(i_box)-Y_edge(JJPAR+1) )
       end do

       do while (box_lat(i_box) < Y_edge(1))
          box_lat(i_box) = Y_edge(1) + ( box_lat(i_box)-Y_edge(1) )
       end do

       do while (box_lon(i_box) > X_edge(IIPAR+1))
          box_lon(i_box) = box_lon(i_box) - 360.0
       end do

       do while (box_lon(i_box) < X_edge(1))
          box_lon(i_box) = box_lon(i_box) + 360.0
       end do


       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)        ! hPa

       ! calcualte the box_alpha
       ! Next adjacent point
       if(i_box==n_boxes_max)then
         next_lon      = box_lon(i_box-1)
         next_lat      = box_lat(i_box-1)
         next_pressure = box_lev(i_box-1)        ! hPa

       else
         next_lon      = box_lon(i_box+1)
         next_lat      = box_lat(i_box+1)
         next_pressure = box_lev(i_box+1)        ! hPa

       endif

       ! Judge the sign of box_alpha
       if((curr_lon-next_lon).NE.0.0)then
         box_alpha = ATAN( (next_lat - curr_lat) / (next_lon - curr_lon) )
       else if((next_lat-curr_lat).LT.0.0)then
         box_alpha = 0.5*PI
       else if((next_lat-curr_lat).GT.0.0)then
         box_alpha = -0.5*PI
       endif


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2) 
       i_lev = Find_iPLev(curr_pressure,P_edge)

       backgrd_concnt(i_box,:) = State_Chm%Species(i_lon,i_lat,i_lev,:)


       Plume_I(i_box) = i_lon
       Plume_J(i_box) = i_lat
       Plume_L(i_box) = i_lev


       curr_u    = Interplt_wind(u,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_v    = Interplt_wind(v,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_omeg = Interplt_wind(omeg,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

       curr_Ptemp = Interplt_wind(Ptemp,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)


       i_lon = Find_iLonLat(next_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(next_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(next_pressure,P_edge)

       !!!!!!
       ! For deformation of cross-section caused by wind shear (A.D.Naiman et al., 2010):
       !!!!!!

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! attention ***
       wind_s_shear = 0.0


       theta_previous   = box_theta(i_box)
       box_theta(i_box) = ATAN( TAN(box_theta(i_box)) + wind_s_shear*Dt )


       do i_ring=1,n_rings_max
       ! make sure use box_theta or TAN(box_theta)  ??? 
       box_radiusA(i_box,i_ring) = box_radiusA(i_box,i_ring) &
                        * (TAN(box_theta(i_box))**2+1)**0.5 / (TAN(theta_previous)**2+1)**0.5
       box_radiusB(i_box,i_ring) = box_radiusB(i_box,i_ring) &
                        * (TAN(box_theta(i_box))**2+1)**(-0.5) / (TAN(theta_previous)**2+1)**(-0.5)
       enddo


       !!!!!!
       ! For the concentration change caused by eddy diffusion:
       ! box_concnt(n_boxes_max,N_rings,N_species)
       !!!!!!

         ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
         Cv = 0.2
         Omega_N = 0.1
         Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
         eddy_v = Cv * Omega_N**2 / sqrt( (Ptemp_shear*g0/curr_Ptemp) )

         ! attenttion ***
         eddy_v = 1.0


         ! Calculate horizontal eddy diffusivity:
         Ch = 0.1
         U_shear = Vertical_shear(u, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         V_shear = Vertical_shear(v, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         UV_shear = sqrt( U_shear**2 + V_shear**2 )
         do i_ring=1,n_rings_max
          !eddy_h(i_ring) = Ch*UV_shear*(Init_radius+(i_ring-1)*D_radius)**2
          ! attention ***
          eddy_h(i_ring) = 1.0

          eddy_A(i_ring) = eddy_v*cos(box_theta(i_box)) + eddy_h(i_ring)*sin(abs(box_theta(i_box))) ! a
          eddy_B(i_ring) = eddy_v*sin(abs(box_theta(i_box))) + eddy_h(i_ring)*cos(box_theta(i_box)) ! b
         enddo


         ! For the innest ring (i_ring = 1)
         ! kB should be rewrite in a more accurate equation !!!
         kB(1) = eddy_B(1) / ( box_radiusB(i_box,2) / 2.0 )
         kA(1) = eddy_A(1) / ( box_radiusA(i_box,2) / 2.0 )

         do i_ring = 2, n_rings_max-1

         kB(i_ring) = eddy_B(i_ring) / ((box_radiusB(i_box,i_ring+1)-box_radiusB(i_box,i_ring-1)) / 2.0 )
         kA(i_ring) = eddy_A(i_ring) / ((box_radiusA(i_box,i_ring+1)-box_radiusA(i_box,i_ring-1)) / 2.0 )

         enddo ! i_ring

         ! For outest ring (i_ring = n_rings_max)
         kB(n_rings_max) = eddy_B(n_rings_max) &
                /( box_radiusB(i_box,n_rings_max)-box_radiusB(i_box,n_rings_max-1) )
         kA(n_rings_max) = eddy_A(n_rings_max) &
                /( box_radiusA(i_box,n_rings_max)-box_radiusA(i_box,n_rings_max-1) )


       do i_species = 1,State_Chm%nSpecies

       Dt2 = 0.1
       do t1s=1,int(Dt/Dt2)

       !==========================================================================
       ! Use classical Runge-Kutta method (RK4) to solve the diferential equation
       !==========================================================================
       do Ki = 1,4

         do i_ring = 1, n_rings_max
           if(Ki==1)then
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species)
           else if(Ki==4)then
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species) + RK(3,i_ring)*Dt2 !Dt
           else
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species) + RK(Ki-1,i_ring)*Dt2*0.5 !Dt
           endif
         enddo ! i_ring

         ! For innest ring:
         AA(1) = kA(1) * 2.0* ( box_concnt_K(2) - box_concnt_K(1) ) 
         BB(1) = 0.0
         DD(1) = D_radius

         RK(Ki,1)          = ( AA(1)+BB(1) ) / DD(1)


         ! For rings from 2 to (n_rings_max - 1)
         do i_ring = 2, n_rings_max-1

           AA(i_ring) = kA(i_ring) *  box_radiusA(i_box,i_ring) &
                        * ( box_concnt_K(i_ring+1) - box_concnt_K(i_ring) ) 
           BB(i_ring) = kA(i_ring-1) * box_radiusA(i_box,i_ring-1) &
                        * ( box_concnt_K(i_ring-1) - box_concnt_K(i_ring) )
           DD(i_ring) = 0.5*( box_radiusA(i_box,i_ring) + box_radiusA(i_box,i_ring-1) ) * D_radius

           RK(Ki,i_ring) = ( AA(i_ring)+BB(i_ring) ) / DD(i_ring)

         enddo ! i_ring
 
         ! For outest ring:
         AA(n_rings_max) = kA(n_rings_max) * box_radiusA(i_box,n_rings_max) &
            * ( State_Chm%Species(i_lon,i_lat,i_lev,i_species) - box_concnt_K(n_rings_max) ) 
         BB(n_rings_max) = kA(n_rings_max-1) * box_radiusA(i_box,n_rings_max-1) &
            * ( box_concnt_K(n_rings_max-1) - box_concnt_K(n_rings_max) )
         DD(n_rings_max) = 0.5 * D_radius &
            * ( box_radiusA(i_box,n_rings_max) + box_radiusA(i_box,n_rings_max-1) )


         RK(Ki,n_rings_max) = (AA(n_rings_max)+BB(n_rings_max))/DD(n_rings_max)

         AA_env(Ki) = -1.0*AA(n_rings_max)


       enddo ! Ki

       if(MOD(t1s,10)==0)then
       if(i_box==1)then
!        write(6,*) '= RK 1 =>', RK(1,1), RK(2,1), RK(3,1), RK(4,1)
!        write(6,*) '= RK 2 =>', RK(1,10), RK(2,10), RK(3,10), RK(4,10)
       endif
       endif

       do i_ring = 1,n_rings_max
         box_concnt(i_box,i_ring,i_species) = box_concnt(i_box,i_ring,i_species) &
                + Dt2*( RK(1,i_ring)+2.0*RK(2,i_ring)+2.0*RK(3,i_ring)+RK(4,i_ring) )/6.0 ! Dt
! box_concnt(i_box,i_ring,i_species) = box_concnt(i_box,i_ring,i_species) + Dt2 * RK(1,i_ring) ! Dt
       enddo !i_ring

       ! env_amount is used to evalue whether mass is conserved or not
        env_amount(i_box) = env_amount(i_box) &
               + Dt2*( AA_env(1)+2.0*AA_env(2)+2.0*AA_env(3)+AA_env(4) )/6.0

       ! Update the concentration in the background grid cell
       ! after the interaction with plume
       grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
       exchange_amount = Dt2* ( AA_env(1)+2.0*AA_env(2)+2.0*AA_env(3)+AA_env(4) ) /6.0
       backgrd_concnt(i_box,i_species) = &
         ( backgrd_concnt(i_box,i_species)*grid_volumn + exchange_amount ) / grid_volumn

       State_Chm%Species(i_lon,i_lat,i_lev,i_species) = backgrd_concnt(i_box,i_species)

       Plume_Species(i_box,i_ring,i_species) = box_concnt(i_box,i_ring,i_species)

        if(MOD(t1s,10)==0)then
        if(i_box==1)then

         OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='OLD', &
               FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!         write(262,*)box_concnt(1,:,1)

        endif
        endif

       enddo ! t1s
       enddo ! do i_Species = 1,nSpecies

    enddo ! i_box

    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'plume_mod.F90' )
       RETURN
    ENDIF


    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)
    nullify(omeg)
    nullify(X_edge)
    nullify(Y_edge)
    nullify(P_edge)
    nullify(X_mid)
    nullify(Y_mid)
    nullify(P_mid)

  END SUBROUTINE plume_run



!--------------------------------------------------------------------
! functions to find location (i,j,k) of boxes 

  integer function Find_iLonLat(curr_xy,  Dxy,  XY_edge2)
    implicit none
    real(fp) :: curr_xy, Dxy, XY_edge2
    Find_iLonLat = INT( (curr_xy - (XY_edge2 - Dxy)) / Dxy )+1
    ! Notice the difference between INT(), FLOOR(), AINT()
    ! for lon: Xedge_Sec - Dx = Xedge_first
    return
  end function



!-------------------------------------------------------------------
  integer function Find_iPLev(curr_pressure,P_edge)
    implicit none
    real(fp) :: curr_pressure
    real(fp), pointer :: P_edge(:)
    integer :: i_lon, i_lat
    integer :: locate(1)
    locate = MINLOC(abs( P_edge(:)-curr_pressure ))

    if(P_edge(locate(1))-curr_pressure >= 0 )then
       Find_iPLev = locate(1)
    else
       Find_iPLev = locate(1) - 1
    endif

    return
  end function

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) into a specific location (curr_lon, curr_lat, curr_pressure)
! based on the surrounding 4 points.

  real function Interplt_wind(wind, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          ::  distance(2,2), Weight(2,2)
    real(fp)          ::wind_lonlat(2), wind_lonlat_lev

    ! Identify wether particle is exactly located on the grid point
    if(curr_lon==X_mid(i_lon))then
      if(curr_lat==Y_mid(i_lat))then
        if(curr_pressure==P_mid(i_lev))then

          Interplt_wind = wind(i_lon, i_lat, i_lev)
          return

        endif
      endif
    endif

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southeast of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
        ii = i + init_lon - 1
        jj = j + init_lat - 1

        ! For some special circumstance:
        if(ii==0)then
          distance(i,j) = Distance_Circle(curr_lon+360.0, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
        else if(ii==(IIPAR+1))then
          distance(i,j) = Distance_Circle(curr_lon-360.0, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
        else
          distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
        endif

    enddo
    enddo

    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
    enddo
    enddo


    do k = 1,2
      kk = k + init_lev - 1
      wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk)      &
                         + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                         + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                         + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
    enddo


    ! second interpolate vertically (Linear)

    wind_lonlat_lev = wind_lonlat(1)    &
                    + (wind_lonlat(2)-wind_lonlat(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) &
                    * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev),
    !P_mid(i_lev+1), curr_pressure )

    Interplt_wind = wind_lonlat_lev

    return
  end function



  !-------------------------------------------------------------------
  ! calculation the great-circle distance between two points on the earth surface

  real function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2  ! unit is degree
    !real(fp) :: PI, Re

    Distance_Circle = Re * 2.0 * ASIN( (sin( (y1-y2)*PI/180.0 ))**2.0   & 
                      + cos(x1*PI/180.0) * cos(x2*PI/180.0) * (sin( 0.5*(x1-x2)*PI/180.0 ))**2.0 )
    return
  end function


  !-------------------------------------------------------------------  
  ! calculate the wind_s (inside a plume sross-section) shear along pressure direction

  real function Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge,  i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: box_alpha
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: u(:,:,:), v(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:), P_edge(:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: u_lonlat(2), v_lonlat(2)
    real(fp)          :: wind_s(2),  Delt_height


    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southeast of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    do i = 1,2
      do j = 1,2

        ii = i + init_lon - 1
        jj = j + init_lat - 1
        distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))

      enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2

      kk = k + init_lev - 1

      u_lonlat(k) =  Weight(1,1) * u(i_lon,i_lat,kk)   + Weight(1,2) * u(i_lon,i_lat+1,kk)   &
                       + Weight(2,1) * u(i_lon+1,i_lat,kk) + Weight(2,2) * u(i_lon+1,i_lat+1,kk)

      v_lonlat(k) =  Weight(1,1) * v(i_lon,i_lat,kk)   + Weight(1,2) * v(i_lon,i_lat+1,kk)   &
                       + Weight(2,1) * v(i_lon+1,i_lat,kk) + Weight(2,2) * v(i_lon+1,i_lat+1,kk)


      if(box_alpha>=0.0)then
        wind_s(k) = -1.0*u_lonlat(k)*SIN(box_alpha) + v_lonlat(k)*COS(box_alpha)
      else 
        wind_s(k) = u_lonlat(k)*SIN(box_alpha) + v_lonlat(k)*COS(box_alpha)
      endif

    enddo


    ! second vertical shear of wind_s

    ! This code should be changed !!!
    ! Because it is the poressure center in [hPa] instead of height center in [m]
    ! Delt_height    = 0.5 * ( P_BXHEIGHT(init_lon,init_lat,init_lev) + P_BXHEIGHT(init_lon,init_lat,init_lev+1) )
    Delt_height =  Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), P_edge(init_lev), P_edge(init_lev+1), 1 )   &
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1), P_edge(init_lev), P_edge(init_lev+1), 0 )
    ! find the z height of each pressure level in GEOS-Chem

    Wind_shear_s = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function


!-------------------------------------------------------------------
  ! transform the pressure level [Pa] to the height level [m]

  real function Pa2meter(Box_height, P1, P2, Judge) 
  ! Judge: 0 for bottom, 1 for top
  ! 
    implicit none
    real(fp) :: Box_height, P1, P2
    integer  :: Judge
      
      if(Judge==1)then
        ! calculate the height of the top half of the grid box
        Pa2meter = Box_height * ( DLOG(P2) - DLOG(0.5*(P2+P1)) ) / ( DLOG(P2) - DLOG(P1) )
      else
        ! calculate the height of the bottom half of the grid box
        Pa2meter = Box_height * ( DLOG(0.5*(P2+P1)) - DLOG(P1) ) / ( DLOG(P2) - DLOG(P1) )
      endif

    return
  end function


!-------------------------------------------------------------------
! Calculate eddy diffusivity in stratisphere. (U.Schumann, 2012)



  !-------------------------------------------------------------------  
  ! calculate the vertical shear for calculating eddy difussivity

  real function Vertical_shear(var, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge,  i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: var(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:), P_edge(:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: var_lonlat(2)
    real(fp)          :: Delt_height


    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southeast of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    do i = 1,2
      do j = 1,2

        ii = i + init_lon - 1
        jj = j + init_lat - 1
        distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii),Y_mid(jj))

      enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2
      kk            = k + init_lev - 1
      var_lonlat(k) =  Weight(1,1) * var(i_lon,i_lat,kk)   &
                       + Weight(1,2) * var(i_lon,i_lat+1,kk)   &
                       + Weight(2,1) * var(i_lon+1,i_lat,kk) &
                       + Weight(2,2) * var(i_lon+1,i_lat+1,kk)
    enddo


    ! second vertical shear of wind

    Delt_height =  Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev),P_edge(init_lev), P_edge(init_lev+1), 1 )   &
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),P_edge(init_lev+1), P_edge(init_lev+2), 0 )
    ! find the z height of each pressure level in GEOS-Chem

    Vertical_shear = ( var_lonlat(2) - var_lonlat(1) ) / Delt_height

    return
  end function


!-------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------


  SUBROUTINE plume_write_std( am_I_Root, RC )


    USE m_netCDF_io_define
    USE m_netcdf_io_read
    USE m_netcdf_io_open
    USE Ncdf_Mod,            ONLY : NC_Open
    USE Ncdf_Mod,            ONLY : NC_Read_Time
    USE Ncdf_Mod,            ONLY : NC_Read_Arr
    USE Ncdf_Mod,            ONLY : NC_Create
    USE Ncdf_Mod,            ONLY : NC_Close
    USE Ncdf_Mod,            ONLY : NC_Var_Def
    USE Ncdf_Mod,            ONLY : NC_Var_Write
    USE Ncdf_Mod,            ONLY : NC_Get_RefDateTime
    USE CHARPAK_Mod,         ONLY : TRANLC
    USE JulDay_Mod,          ONLY : JulDay

    USE Input_Opt_Mod, ONLY : OptInput

    ! Parameters for netCDF routines
    include "netcdf.inc"


    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success

    REAL(f8), POINTER         :: Arr1D(:,:)
    REAL(f8), POINTER         :: Arr2D(:,:)
    REAL(f8), POINTER         :: Arr2DOld(:,:)
    REAL(f4), POINTER         :: Arr3D(:,:)
    REAL(f4), POINTER         :: Arr4D(:,:)
    REAL*8,   POINTER         :: timeVec(:)
    CHARACTER(LEN=255)        :: NcFile
    CHARACTER(LEN=255)        :: Pfx, title, Reference, Contact
    CHARACTER(LEN=255)        :: timeunit
    INTEGER                   :: fId, lonId, latId, levId, TimeId
    INTEGER                   :: VarCt
    INTEGER                   :: nLon, nLat, nLev, nTime
    INTEGER                   :: i_box, i_ring
!    INTEGER                   :: L
    LOGICAL                   :: IsOldFile

    CHARACTER(LEN=255)            :: FILENAME2

    FILENAME2   = 'Plume_theta_max_min_radius.txt'
    tt = tt +1

!    IF(mod(tt,6)==0)THEN     ! output once every hour
!    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)

!       OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='OLD', &
!             FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!       Do i_ring = 1, n_rings_max
!        WRITE(262,'(I0.4,3(x,E16.5E4))') i_ring, box_theta(2), box_radiusA(2,i_ring), box_radiusB(2,i_ring), box_concnt(2,i_ring,1)
!       End Do


!    ENDIF

  END SUBROUTINE plume_write_std

!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  subroutine plume_cleanup()

    if (allocated(box_lon))      deallocate(box_lon)
    if (allocated(box_lat))      deallocate(box_lat)
    if (allocated(box_lev))      deallocate(box_lev)
    if (allocated(box_length))   deallocate(box_length)
    if (allocated(box_radiusA))  deallocate(box_radiusA)
    if (allocated(box_radiusB))  deallocate(box_radiusB)
    if (allocated(box_theta))    deallocate(box_theta)

    WRITE(6,'(a)') '--> Plume Module Cleanup <--'

  end subroutine plume_cleanup


END MODULE Plume_Mod
