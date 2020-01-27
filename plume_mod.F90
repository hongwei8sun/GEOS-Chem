!------------------------------------------------------------------------------
!   GEOS-Chem Global Chemical Stratospherical Transport Model                  
!------------------------------------------------------------------------------

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
  PUBLIC :: box_concnt ! (n_boxes_max, n_rings_max, State_Chm%nSpecies)


  integer, allocatable  :: Plume_I(:), Plume_J(:), Plume_L(:) ! (n_boxes_max)

  integer, parameter    :: n_boxes_max = 3  
  integer               :: tt
  integer               :: N_inject

  ! Consider the box as a cylinder
  ! The center of the cylinder
  real(fp), allocatable :: box_lon(:)
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)

  ! box_radius(n_boxes_max,N_rings) of the cross-section
  real(fp), allocatable :: box_radiusA(:,:) ! vertical radius
  real(fp), allocatable :: box_radiusB(:,:) ! horizontal radius

  ! Theta is the clockwise angle between z-axis (P) and vertical radiusA
  real(fp), allocatable :: box_theta(:)     ! 0 ~ 180 degree
  real(fp), allocatable :: box_length(:)

  ! D_radius should only be used at the beginning!
  real(fp), parameter   :: Init_radius = 100.0e+0_fp ! [m]
  real(fp), parameter   :: D_radius    = 100.0e+0_fp ! [m], the width of ring
  integer, parameter    :: n_rings_max = 500 ! number of rings in one box

  ! medical concentration of each ring [molec/cm3]
  real(fp), allocatable :: box_concnt(:,:,:) !(n_boxes_max,N_rings,N_species)
  real(fp), allocatable :: box_concnt_K(:) ! box_concnt_K(N_rings_max)

  real(fp), allocatable :: RK(:,:), Outer2env(:)  ! for Runge-Kutta method

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
    allocate(Outer2env(4))

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

    box_theta         = (/0.0e+0_fp,    0.0e+0_fp,    0.0e+0_fp/) ! degree
    box_length        = (/1000.0e+0_fp, 1000.0e+0_fp, 1000.0e+0_fp/) ! m



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


    !======================================================================
    ! identify the location of the plume box, 
    ! put the corresponding concentration of backgound grid into the initial 
    ! concentration inside plume in unit of [molec/cm3].
    !======================================================================
    do i_box = 1,n_boxes_max

       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)        ! hPa

       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(curr_pressure,P_edge)

       do i_ring = 1,N_species
       do i_species = 1,N_species
          box_concnt(i_box,i_ring,i_species) = &
                             State_Chm%Species(i_lon,i_lat,i_lev,i_species) 
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
    WRITE(6,*)'=== the unit for species currently is:', State_Chm%Spc_Units


    box_concnt_K(:) = 0.0e+0_fp ! [n_ring_max]


    ! State_Chm%nAdvect: the last one is PASV
    N_inject = State_Chm%nAdvect

    ! test:
    ! Sigma = 1000m, Dr = 100m, 
    box_concnt(1,1:100,N_inject) = 1.0e-20_fp * (/ 49.75062395963412, 47.799874091655, &
        44.124845129229776, 39.13522691209341, 33.34884054292372, 27.30372133198547, &
        21.477867910536958, 16.232623367917487, 11.787303827793176, 8.223722828857746, &
        5.512526265224261, 3.550267686981849, 2.196846681170371, 1.3060704926959117, &
        0.7460393034533921, 0.4094350507187042, 0.2158920003816539, 0.10937455590914426, &
        0.05323831183339602, 0.024897771075163673, 0.011187289686031028, &
        0.0048296706861092, 0.0020032648696475534, 0.0007983391948902373, &
        0.0003056783983185701, 0.00011245279835161713, 3.97469680767457e-05, &
        1.349789251681507e-05, 4.404089598230278e-06, 1.3806212284140175e-06, &
        4.158351228414166e-07, 1.203361218150526e-07, 3.345793045646391e-08, &
        8.937794355639906e-09, 2.293981243569646e-09, 5.656888100032663e-10, &
        1.3402738186563035e-10, 3.050968338802662e-11, 6.672830424109946e-12, &
        1.4022023691113985e-12, 2.8309977584744754e-13, 5.491570649142986e-14, &
        1.0234858565821022e-14, 1.8327166978005226e-15, 3.1530949469932235e-16, &
        5.2120308919508e-17, 8.277613310493872e-18, 1.2630818904628463e-18, &
        1.8517659888261724e-19, 2.6083683310135188e-20, 3.530042668631253e-21, &
        4.590068975479606e-22, 5.734382911127981e-23, 6.883073193309613e-24, &
        7.937912491322717e-25, 8.795457748423832e-26, 9.363512772281097e-27, &
        9.577394759835739e-28, 9.412049237974876e-29, 8.88687789784163e-30, &
        8.061993740251162e-31, 7.026902388814562e-32, 5.884554719608361e-33, &
        4.734689586357957e-34, 3.660139495529953e-35, 2.7185165709661418e-36, &
        1.9399679381421983e-37, 1.3301032079719932e-38, 8.762022218692594e-40, &
        5.545638806275141e-41, 3.3723064310643125e-42, 1.970292490713008e-43, &
        1.1060190168004981e-44, 5.965168380036516e-46, 3.091085662064773e-47, &
        1.5389609318814287e-48, 7.361602611243622e-50, 3.383337837584202e-51, &
        1.4939861871492439e-52, 6.338347115385447e-54, 2.5836499800303364e-55, &
        1.0118579679040722e-56, 3.8074451468681845e-58, 1.37649923807487e-59, &
        4.7813057471765284e-61, 1.5956780746337997e-62, 5.116491534769963e-64, &
        1.5762585775725126e-65, 4.665636357469077e-67, 1.326852145501891e-68, &
        3.6254539752288515e-70, 9.517666478947988e-72, 2.400638620885891e-73, &
        5.817699914767302e-75, 1.3545781439450351e-76, 3.0302958619542277e-78, &
        6.513197121451957e-80, 1.345028951291385e-81, 2.6686846856493707e-83, &
        5.0873441170807215e-85 /)

    env_amount = (/0.0e+0_fp, 0.0e+0_fp, 0.0e+0_fp/)

    ! Create output file
    FILENAME2   = 'Plume_theta_max_min_radius.txt'
    tt = 0

    OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       WRITE(262,*)box_concnt(1,:,N_inject)

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
    ! XEDGE( IM+1, JM,   L ), YEDGE( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

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

    real(fp)          :: Outer(n_rings_max), Inner(n_rings_max), V_ring(n_rings_max) 
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

    u => State_Met%U ! [m/s]
    v => State_Met%V ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    Ptemp => State_Met%THETA  ! Updraft velocity [Pa/s]

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]
    P_mid => State_Met%PMID(1,1,:) ![hPa] Pressure (moist air)@level centers

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
          box_lat(i_box) = Y_edge(JJPAR+1) &
                         - ( box_lat(i_box)-Y_edge(JJPAR+1) )
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
       ! [molec/cm3]


       Plume_I(i_box) = i_lon
       Plume_J(i_box) = i_lat
       Plume_L(i_box) = i_lev


       curr_u = Interplt_wind(u, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_v = Interplt_wind(v, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_omeg = Interplt_wind(omeg, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

       curr_Ptemp = Interplt_wind(Ptemp, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)


       i_lon = Find_iLonLat(next_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(next_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(next_pressure,P_edge)

       !=====================================================================
       ! For deformation of cross-section caused by wind shear 
       ! (A.D.Naiman et al., 2010):
       !=====================================================================

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! *** attention *** ???
       wind_s_shear = 0.0


       theta_previous   = box_theta(i_box)
       box_theta(i_box) = ATAN( TAN(box_theta(i_box)) + wind_s_shear*Dt )


       do i_ring=1,n_rings_max
       ! make sure use box_theta or TAN(box_theta)  ??? 
       box_radiusA(i_box,i_ring) = box_radiusA(i_box,i_ring) &
                                  * (TAN(box_theta(i_box))**2+1)**0.5 &
                                  / (TAN(theta_previous)**2+1)**0.5
       box_radiusB(i_box,i_ring) = box_radiusB(i_box,i_ring) &
                                  * (TAN(box_theta(i_box))**2+1)**(-0.5) &
                                  / (TAN(theta_previous)**2+1)**(-0.5)
       enddo


       !=====================================================================
       ! For the concentration change caused by eddy diffusion:
       ! box_concnt(n_boxes_max,N_rings,N_species)
       !=====================================================================

       ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
       Cv = 0.2
       Omega_N = 0.1
       Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! *** attention *** ???
       eddy_v = Cv * Omega_N**2 / sqrt( (Ptemp_shear*g0/curr_Ptemp) )
       eddy_v = 1.0


       ! Calculate horizontal eddy diffusivity:
       Ch = 0.1
       U_shear = Vertical_shear(u, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       V_shear = Vertical_shear(v, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       UV_shear = sqrt( U_shear**2 + V_shear**2 )
       do i_ring=1,n_rings_max
        ! *** attention *** ???
          eddy_h(i_ring) = Ch*UV_shear*(Init_radius+(i_ring-1)*D_radius)**2
          eddy_h(i_ring) = 1.0

          eddy_A(i_ring) = eddy_v*cos(box_theta(i_box)) &
                          + eddy_h(i_ring)*sin(abs(box_theta(i_box))) ! a
          eddy_B(i_ring) = eddy_v*sin(abs(box_theta(i_box))) &
                          + eddy_h(i_ring)*cos(box_theta(i_box)) ! b
       enddo


       !=================================================================
       ! Calculate the transport rate
       ! kB should be rewrite in a more accurate equation !!!
       !=================================================================
       ! For innest ring:
!       kB(1) = eddy_B(1) / ( box_radiusB(i_box,1) &
!                   + 0.5*(box_radiusB(i_box,2)-box_radiusB(i_box,1)) )
!       kA(1) = eddy_A(1) / ( box_radiusA(i_box,1) &
!                   + 0.5*(box_radiusA(i_box,2)-box_radiusA(i_box,1)) )

       kB(1) = eddy_B(1) / ( 0.5*box_radiusB(i_box,2) )
       kA(1) = eddy_A(1) / ( 0.5*box_radiusA(i_box,2) )


       do i_ring = 2, n_rings_max-1

       kB(i_ring) = eddy_B(i_ring) &
         / ( 0.5*(box_radiusB(i_box,i_ring+1)-box_radiusB(i_box,i_ring-1)) )
       kA(i_ring) = eddy_A(i_ring) &
         / ( 0.5*(box_radiusA(i_box,i_ring+1)-box_radiusA(i_box,i_ring-1)) )

       enddo ! 

       ! For outest ring (i_ring = n_rings_max)
       kB(n_rings_max) = eddy_B(n_rings_max) &
        /( box_radiusB(i_box,n_rings_max)-box_radiusB(i_box,n_rings_max-1) )
       kA(n_rings_max) = eddy_A(n_rings_max) &
        /( box_radiusA(i_box,n_rings_max)-box_radiusA(i_box,n_rings_max-1) )


       !===================================================================
       ! Calculate the volume of each ring
       !===================================================================
       V_ring(1) = PI * box_radiusA(i_box,1) * box_radiusB(i_box,1) &
                * box_length(i_box) * 1.0e+6_fp ! [cm3]

       DO i_ring = 2, n_rings_max
          V_ring(i_ring) = 1.0e+6_fp* box_length(i_box) * PI* & ! [cm3]
              ( box_radiusA(i_box,i_ring)*box_radiusB(i_box,i_ring) &
              - box_radiusA(i_box,i_ring-1)*box_radiusB(i_box,i_ring-1) )
       ENDDO


       !==================================================================
       ! Begin to calculate the dilution inside plume 
       !==================================================================
       do i_species = 1,State_Chm%nSpecies

       Dt2 = 1.0
       do t1s=1,int(Dt/Dt2)

       ! Use Runge-Kutta method (RK4) to solve diferential equation
       do Ki = 1,4

       do i_ring = 1, n_rings_max
         if(Ki==1)then
           box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species)
         else if(Ki==4)then
           box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species) &
                                + RK(3,i_ring)*Dt2 !Dt
         else
           box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species) &
                                + RK(Ki-1,i_ring)*Dt2*0.5 !Dt
         endif
       enddo ! i_ring

         ! For innest ring:
         Outer(1) = 1.0e+6_fp * box_length(i_box) * PI &
          * ( (box_radiusA(i_box,1)+0.5*kA(1)*Dt2) &
            * (box_radiusB(i_box,1)+0.5*kB(1)*Dt2) &
            - (box_radiusA(i_box,1)-0.5*kA(1)*Dt2) &
            * (box_radiusB(i_box,1)-0.5*kB(1)*Dt2) ) &
          * ( box_concnt_K(2) - box_concnt_K(1) ) ! [molec]

         Inner(1) = 0.0

         RK(Ki,1) = ( Outer(1)+Inner(1) ) / V_ring(1) ! [molec/cm3]

         ! For rings from 2 to (n_rings_max - 1)
         do i_ring = 2, n_rings_max-1
           Outer(i_ring) = 1.0e+6_fp * box_length(i_box) * PI &
              * ( (box_radiusA(i_box,i_ring) + 0.5*kA(i_ring)*Dt2) &
                * (box_radiusB(i_box,i_ring) + 0.5*kB(i_ring)*Dt2) &
                - (box_radiusA(i_box,i_ring) - 0.5*kA(i_ring)*Dt2) &
                * (box_radiusB(i_box,i_ring) - 0.5*kB(i_ring)*Dt2) ) &
              * (box_concnt_K(i_ring+1)-box_concnt_K(i_ring)) ! [molec]


           Inner(i_ring) = 1.0e+6_fp * box_length(i_box) * PI &
              * ( (box_radiusA(i_box,i_ring-1) + 0.5*kA(i_ring-1)*Dt2) &
                * (box_radiusB(i_box,i_ring-1) + 0.5*kB(i_ring-1)*Dt2) &
                - (box_radiusA(i_box,i_ring-1) - 0.5*kA(i_ring-1)*Dt2) &
                * (box_radiusB(i_box,i_ring-1) - 0.5*kB(i_ring-1)*Dt2) ) &
              * (box_concnt_K(i_ring-1)-box_concnt_K(i_ring)) ! [molec]


           RK(Ki,i_ring) = ( Outer(i_ring)+Inner(i_ring) ) / V_ring(i_ring)

        
         enddo ! i_ring
 
         ! For outest ring:
         Outer(n_rings_max) = 1.0e+6_fp * box_length(i_box) * PI &
           * ( (box_radiusA(i_box,n_rings_max) + 0.5*kA(n_rings_max)*Dt2) &
             * (box_radiusB(i_box,n_rings_max) + 0.5*kB(n_rings_max)*Dt2) &
             - (box_radiusA(i_box,n_rings_max) - 0.5*kA(n_rings_max)*Dt2) &
             * (box_radiusB(i_box,n_rings_max) - 0.5*kB(n_rings_max)*Dt2) ) &
          * (backgrd_concnt(i_box,i_species) - box_concnt_K(n_rings_max) )

         Inner(n_rings_max) = 1.0e+6_fp * box_length(i_box) * PI &
         * ( (box_radiusA(i_box,n_rings_max-1)+0.5*kA(n_rings_max-1)*Dt2) &
          * (box_radiusB(i_box,n_rings_max-1)+0.5*kB(n_rings_max-1)*Dt2) &
          - (box_radiusA(i_box,n_rings_max-1)-0.5*kA(n_rings_max-1)*Dt2) & 
          * (box_radiusB(i_box,n_rings_max-1)-0.5*kB(n_rings_max-1)*Dt2) ) &
         * (box_concnt_K(n_rings_max-1)-box_concnt_K(n_rings_max)) ! [molec]

         RK(Ki,n_rings_max) = ( Outer(n_rings_max) + Inner(n_rings_max) ) &
                              / V_ring(n_rings_max)

         Outer2env(Ki) = -1.0 * Outer(n_rings_max) ! [molec]

       enddo ! Ki

       ! test for plume model:
       if(MOD(t1s,600)==0)then
       if(i_box==1)then
!        write(6,*) '= success in plume model =>', t1s
!        write(6,*) '= RK 2 =>', RK(1,10), RK(2,10), RK(3,10), RK(4,10)
       endif
       endif

       do i_ring = 1,n_rings_max
          box_concnt(i_box,i_ring,i_species) = &
                box_concnt(i_box,i_ring,i_species) &
                        + ( RK(1,i_ring) + 2.0*RK(2,i_ring) &
                           + 2.0*RK(3,i_ring) + RK(4,i_ring) ) /6.0 
       enddo !i_ring

       ! env_amount is used to evalue whether mass is conserved or not
       env_amount(i_box) = env_amount(i_box) & ! [molec]
          + ( Outer2env(1) + 2.0*Outer2env(2) &
             + 2.0*Outer2env(3) + Outer2env(4) ) / 6.0

       !================================================================
       ! Update the concentration in the background grid cell
       ! after the interaction with plume
       !================================================================
       grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]

       exchange_amount = ( Outer2env(1) + 2.0*Outer2env(2) &
                          + 2.0*Outer2env(3) + Outer2env(4) ) / 6.0

       backgrd_concnt(i_box,i_species) = ( exchange_amount + &
                backgrd_concnt(i_box,i_species)*grid_volumn ) / grid_volumn

       State_Chm%Species(i_lon,i_lat,i_lev,i_species) = backgrd_concnt(i_box,i_species)


       !================================================================
       ! test: Output the concentration into a txt file
       !================================================================
!       if(MOD(t1s,10)==0)then
!       if(i_box==1)then
!          OPEN( 262, FILE=TRIM( FILENAME2   ), STATUS='OLD', &
!                         FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!             WRITE(262,*)box_concnt(1,:,N_inject)
!       endif
!       endif

       enddo ! t1s
       enddo ! do i_Species = 1,nSpecies

       ! test
!       if(i_box==1)then
!          WRITE(6,*)'= plume =>',N_inject, box_concnt(1,1:20,N_inject)
!       endif

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

    IF(mod(tt,6)==0)THEN     ! output once every hour
    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)

       OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='OLD', &
             FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       Do i_ring = 1, n_rings_max
          WRITE(262,*)box_concnt(1,:,N_inject)
!          WRITE(262,'(I0.4,4(x,E16.5E4))') i_ring, box_theta(2), box_radiusA(2,i_ring), box_radiusB(2,i_ring), box_concnt(2,i_ring,N_inject)
       End Do

    ENDIF
    ENDIF

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
