!-------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!--------------------------------------------------------------------------
MODULE Lagrange_Mod

  USE precision_mod
  USE ERROR_MOD
  USE ErrCode_Mod
  USE PhysConstants, ONLY : PI, Re, g0, AIRMW, AVO
  USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR

  USE TIME_MOD,        ONLY : GET_YEAR
  USE TIME_MOD,        ONLY : GET_MONTH
  USE TIME_MOD,        ONLY : GET_DAY
  USE TIME_MOD,        ONLY : GET_HOUR
  USE TIME_MOD,        ONLY : GET_MINUTE
  USE TIME_MOD,        ONLY : GET_SECOND

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: lagrange_init
  PUBLIC :: lagrange_run
  PUBLIC :: plume_run
  PUBLIC :: lagrange_write_std
  PUBLIC :: lagrange_cleanup

  ! PUBLIC VARIABLES:
  PUBLIC :: box_lon, box_lat, box_lev
  PUBLIC :: n_boxes_max
  PUBLIC :: N_curr, N_prev
  PUBLIC :: box_u, box_v, box_omeg
  PUBLIC :: box_Ptemp

  PUBLIC :: Plume_I, Plume_J, Plume_L ! (n_boxes_max)
  PUBLIC :: n_rings_max
  PUBLIC :: Max_rings
  PUBLIC :: box_concnt ! (n_boxes_max, n_rings_max, N_species)

  PUBLIC :: Extra_amount        ! [molec]

  real(fp), allocatable :: box_lon(:), box_lat(:), box_lev(:)
  integer, parameter    :: n_boxes_max = 9999 !6904224     
  integer               :: N_curr, N_prev
  real(fp), allocatable :: box_u(:), box_v(:), box_omeg(:)
  real(fp), allocatable :: box_Ptemp(:)

  integer, allocatable  :: Plume_I(:), Plume_J(:), Plume_L(:) !(n_boxes_max)
  integer, allocatable  :: Max_rings(:)
  integer, parameter    :: n_rings_max = 128 ! 2^7, number of rings in one box
  ! medical concentration of each ring [molec/cm3]
  real(fp), allocatable :: box_concnt(:,:,:) !(n_boxes_max,N_rings,N_species)


  real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)
  real(fp), pointer :: P_edge(:)

  REAL(fp), allocatable :: V_ring(:,:) !(n_rings_max)

  integer, parameter    :: N_parcels   = 10 ! 131        
  integer               :: tt     
  ! Aircraft would release 131 aerosol parcels every time step

!  integer               :: i_rec
  integer               :: N_species


  ! box_radius(n_boxes_max,N_rings) of the cross-section
  real(fp), allocatable :: box_radiusA(:,:) ! vertical radius
  real(fp), allocatable :: box_radiusB(:,:) ! horizontal radius

  ! Theta is the clockwise angle between z-axis (P) and vertical radiusA
  real(fp), allocatable :: box_theta(:)     ! 0 ~ 180 degree
  real(fp), allocatable :: box_length(:)

  ! D_radius should only be used at the beginning!
  real(fp), parameter   :: Init_radius = 100.0e+0_fp ! [m]
  real(fp), parameter   :: D_radius    = 100.0e+0_fp ! [m], the width of ring


  ! Used for plume_run subroutine:
  real(fp), allocatable :: box_concnt_K(:) ! box_concnt_K(N_rings_max)
  real(fp), allocatable :: RK(:,:), Outer2env(:)  ! for Runge-Kutta method

  real(fp), allocatable :: kA(:), kB(:)

  real(fp), allocatable :: env_amount(:)
  real(fp), allocatable :: backgrd_concnt(:,:)

  real(fp), allocatable :: Extra_amount(:,:)

CONTAINS


!-----------------------------------------------------------------

  SUBROUTINE lagrange_init(am_I_root, Input_Opt, State_Chm, State_Met, RC)

    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Met_Mod, ONLY : MetState
    USE State_Chm_Mod, ONLY : ChmState

    USE GC_GRID_MOD,   ONLY : XMID, YMID
    USE GC_GRID_MOD,     ONLY : GET_AREA_M2

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR
    USE TIME_MOD,        ONLY : GET_MINUTE
    USE TIME_MOD,        ONLY : GET_SECOND


    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    INTEGER                       :: i_box, i_ring
    INTEGER                       :: ii, jj, kk
    CHARACTER(LEN=255)            :: FILENAME, FILENAME_INIT
    CHARACTER(LEN=255)            :: FILENAME2

    integer :: i_lon, i_lat, i_lev            !1:IIPAR

    REAL(fp) :: lon1, lat1, lon2, lat2
        
    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND

    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C

    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Lagrnage Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'

    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))
    allocate(box_u(n_boxes_max))
    allocate(box_v(n_boxes_max))
    allocate(box_omeg(n_boxes_max))
    allocate(box_Ptemp(n_boxes_max))

    allocate(Max_rings(n_boxes_max))
    allocate(V_ring(n_boxes_max,n_rings_max))
!--------------------------------------------------
! Allocate following variables for plume_run
!--------------------------------------------------
    allocate(Plume_I(n_boxes_max))
    allocate(Plume_J(n_boxes_max))
    allocate(Plume_L(n_boxes_max))

    allocate(box_radiusA(n_boxes_max,n_rings_max))
    allocate(box_radiusB(n_boxes_max,n_rings_max))

    allocate(box_theta(n_boxes_max))
    allocate(box_length(n_boxes_max))

!    N_species = State_Chm%nAdvect
    N_species = 1

    allocate(box_concnt(n_boxes_max,n_rings_max,N_species))
    allocate(box_concnt_K(n_rings_max))
    allocate(RK(4,n_rings_max))
    allocate(Outer2env(4))

    allocate(kA(n_rings_max))
    allocate(kB(n_rings_max))

    allocate(env_amount(n_boxes_max))
    allocate(backgrd_concnt(n_boxes_max,N_species))

    allocate(Extra_amount(n_boxes_max,N_species))

    X_mid  => XMID(:,1,1)   ! IIPAR
    Y_mid  => YMID(1,:,1)
    P_mid  => State_Met%PMID(1,1,:)  ! Pressure at level centers (hPa)

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]

    Max_rings = n_rings_max

!--------------------------------------------------
! Set the initail value of box location, length
!--------------------------------------------------

    ! (1) run month by month
!    FILENAME_INIT   = '/n/home12/hongwei/hongwei/merra2_2x25_standard_Dec/Lagrange_Init_box_i_lon_lat_lev.txt'
!    OPEN( 361,      FILE=TRIM( FILENAME_INIT   ), STATUS='old', &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!    Do i_box = 1, n_boxes_max
!       READ(361,'( 2(x,I8) , 3(x,E12.5) )') N_prev, iibox, box_lon(i_box), box_lat(i_box), box_lev(i_box)
!    End Do
    
!--------------------------------------------------

    ! (2)
    do i_box = 1,n_boxes_max,1
        box_lon(i_box) = -141.0e+0_fp   ! 0   
        box_lat(i_box) = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) &
                        * (-1.0)**FLOOR(i_box/6000.0) ! -29.95S : 29.95N : 0.1
        box_lev(i_box) = 52.0e+0_fp       ! about 20 km
    enddo

!--------------------------------------------------------
! Set the initial value of plume character:
! box_theta, box_length,
!  max/min radius for each ring
!--------------------------------------------------------

    box_theta    = 0.0e+0_fp     ! [degree]

    i_box = 1
    lon1 = box_lon(i_box) - 0.5 * ( box_lon(i_box+1) - box_lon(i_box) )
    lon2 = box_lon(i_box) + 0.5 * ( box_lon(i_box+1) - box_lon(i_box) )
    lat1 = box_lat(i_box) - 0.5 * ( box_lat(i_box+1) - box_lat(i_box) )
    lat2 = box_lat(i_box) + 0.5 * ( box_lat(i_box+1) - box_lat(i_box) )
    box_length(i_box) = Distance_Circle(lon1, lat1, lon2, lat2)

    DO i_box = 2, n_boxes_max-1
      lon1 = 0.5 * ( box_lon(i_box) + box_lon(i_box-1) )
      lon2 = 0.5 * ( box_lon(i_box) + box_lon(i_box+1) )
      lat1 = 0.5 * ( box_lat(i_box) + box_lat(i_box-1) )
      lat2 = 0.5 * ( box_lat(i_box) + box_lat(i_box+1) )
      box_length(i_box)   = Distance_Circle(lon1, lat1, lon2, lat2)  ! [m]
    ENDDO

    i_box = n_boxes_max
    lon1 = box_lon(i_box) - 0.5 * ( box_lon(i_box) - box_lon(i_box-1) )
    lon2 = box_lon(i_box) + 0.5 * ( box_lon(i_box) - box_lon(i_box-1) )
    lat1 = box_lat(i_box) - 0.5 * ( box_lat(i_box) - box_lat(i_box-1) )
    lat2 = box_lat(i_box) + 0.5 * ( box_lat(i_box) - box_lat(i_box-1) )
    box_length(i_box) = Distance_Circle(lon1, lat1, lon2, lat2)


    do i_ring=1,n_rings_max
      box_radiusA(:,i_ring) = i_ring * D_radius
      box_radiusB(:,i_ring) = i_ring * D_radius
    enddo

    tt   = 0
    N_prev = 0
    N_curr = N_prev + N_parcels

!--------------------------------------------------------
! Set the initial concentration of injected aerosols
!--------------------------------------------------------
    ! State_Chm%nAdvect: the last one is PASV
    box_concnt = 0.0e+0_fp

    ! Injection rate is 30 [kg/km] ~ 30 [g/m] 
    DO i_box = 1, n_boxes_max

      box_concnt(i_box,1,N_species) = 3.0e+1_fp * box_length(i_box) &
       /( PI*box_radiusA(i_box,1)*box_radiusB(i_box,1)*box_length(i_box) )

      ! Change from [g/m3] to [molec/cm3],
      ! 98.0 g/mol is the molar mass of H2SO4
       box_concnt(i_box,1,N_species) = box_concnt(i_box,1,N_species)/1.0e+6_fp &
                                    / 98.0 * AVO
    ENDDO


    box_concnt_K(:) = 0.0e+0_fp ! [n_ring_max]
    env_amount      = 0.0e+0_fp

!--------------------------------------------------
! Get time info to create output txt file
!--------------------------------------------------

    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()
    MINUTE      = GET_MINUTE()
    SECOND      = GET_SECOND()

    WRITE(YEAR_C,*) YEAR
    WRITE(MONTH_C,*) MONTH
    WRITE(DAY_C,*) DAY
    WRITE(HOUR_C,*) HOUR
    WRITE(MINUTE_C,*) MINUTE
    WRITE(SECOND_C,*) SECOND


!-----------------------------------------------------------------
! Output box location
!-----------------------------------------------------------------

    FILENAME   = 'Lagrange_Init_xyz_' // TRIM(ADJUSTL(YEAR_C)) &
         // '-' //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) &
         // '-' // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
         // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'

    OPEN( 261,      FILE=TRIM( FILENAME ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    Do i_box = 1, n_boxes_max
       WRITE(261,'( 2(x,I8) , 3(x,E12.5) )')N_prev, i_box, box_lon(i_box), &
                                            box_lat(i_box), box_lev(i_box)
    End Do


!===================================================================
! Calculate the volume of each ring
! V_ring would never change for the whole simulation
!===================================================================
    DO i_box = 1, n_boxes_max

      V_ring(i_box,1) = 1.0e+6_fp* box_length(i_box)* PI &
             * box_radiusA(i_box,1) * box_radiusB(i_box,1)   ! [cm3]

      DO i_ring = 2, Max_rings(i_box)
         V_ring(i_box,i_ring) = 1.0e+6_fp* box_length(i_box)* PI & ! [cm3]
            * ( box_radiusA(i_box,i_ring)*box_radiusB(i_box,i_ring) &
              - box_radiusA(i_box,i_ring-1)*box_radiusB(i_box,i_ring-1) )
      ENDDO

    ENDDO

!-----------------------------------------------------------------
! Output box concentration
!-----------------------------------------------------------------
    FILENAME2   = 'Plume_Init_concentration_molec_' // TRIM(ADJUSTL(YEAR_C)) // '-' &
        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' // &
        TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) // ':' // &
        TRIM(ADJUSTL(SECOND_C)) // '.txt'

    OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    DO i_box = 1, 1001, 100
      WRITE(262,*) i_box
      WRITE(262,*) SUM( box_concnt(i_box,:,N_species)*V_ring(i_box,:) )       ! [molec]
      WRITE(262,*) box_concnt(i_box,:,N_species)
    ENDDO

!     WRITE(262,*) box_concnt(1,:,N_species)

!-----------------------------------------------------------------
! Output State_Met%AD(i_lon,i_lat,i_lev) into State_Met_AD.txt
!-----------------------------------------------------------------
    OPEN( 314,      FILE='State_Met_AD.txt', STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    Do i_lon = 1, IIPAR
    Do i_lat = 1, JJPAR
    Do i_lev = 1, LLPAR
       WRITE(314,'(x,E12.5)') State_Met%AD(i_lon,i_lat,i_lev)
    End Do
    End Do
    End Do

!-----------------------------------------------------------------------
! Output State_Met%AREA_M2(i_lon,i_lat,i_lev) into State_Met_AREA_M2.txt
!-----------------------------------------------------------------------
    OPEN( 315,      FILE='State_Met_AREA_M2.txt', STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    Do i_lon = 1, IIPAR
    Do i_lat = 1, JJPAR
    Do i_lev = 1, LLPAR
       WRITE(315,'(x,E12.5)') GET_AREA_M2(i_lon,i_lat,i_lev)
    End Do
    End Do
    End Do

  END SUBROUTINE lagrange_init

!-----------------------------------------------------------------
!=================================================================

  SUBROUTINE lagrange_run(am_I_Root, State_Chm, State_Met, Input_Opt, RC)

    USE Input_Opt_Mod, ONLY : OptInput

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE TIME_MOD,      ONLY : GET_TS_DYN

    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE              
    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON 

    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    REAL :: Dt, RK_Dt(4)                  ! = 600.0e+0_fp          

    real(fp), pointer  :: PASV_EU           
    integer            :: nAdv        
    REAL(fp)           :: MW_g

    integer :: i_box, N_box

    integer :: i_lon            !1:IIPAR
    integer :: i_lat
    integer :: i_lev

    integer :: ii, jj, kk

    integer :: Ki

    real(fp) :: curr_lon, curr_lat, curr_pressure
    real(fp) :: RK_lon, RK_lat
    real(fp) :: X_edge2
    real(fp) :: Y_edge2

    real(fp) :: dbox_lon, dbox_lat, dbox_lev
    real(fp) :: dbox_x_PS, dbox_y_PS
    real(fp) :: box_x_PS, box_y_PS

    real(fp) :: RK_Dlon(4), RK_Dlat(4), RK_Dlev(4)
    real(fp) :: RK_x_PS, RK_y_PS
    real(fp) :: RK_Dx_PS, RK_Dy_PS

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)
    real(fp), pointer :: Ptemp(:,:,:)

    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: curr_u_PS, curr_v_PS

    real(fp), pointer :: X_edge(:), Y_edge(:)

    real(fp) :: Dx, Dy

    real(fp) :: length0
    real(fp) :: lon1, lon2, lat1, lat2

    real(fp), pointer :: P_I, R_e     

    Dt = GET_TS_DYN()

    RK_Dt(1) = Dt
    RK_Dt(2) = 0.5*Dt
    RK_Dt(3) = 0.5*Dt
    RK_Dt(4) = Dt

    u => State_Met%U   ! figure out state_met%U is based on lat/lon or modelgrid(i,j)
    v => State_Met%V   ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    Ptemp => State_Met%THETA ! Potential temperature [K]

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!
    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)
     
    
    if(N_curr<=n_boxes_max)then
      N_box = N_curr
    else
      N_box = n_boxes_max
    endif

    !-----------------------------------------------------------------------
    ! Run Lagrangian trajectory-track HERE
    !-----------------------------------------------------------------------
    do i_box = 1,N_box,1

       ! make sure the location is not out of range
       do while (box_lat(i_box) > Y_edge(JJPAR+1))
          box_lat(i_box) = Y_edge(JJPAR+1) - ( curr_lat-Y_edge(JJPAR+1) )
       end do

       do while (box_lat(i_box) < Y_edge(1))
          box_lat(i_box) = Y_edge(1) + ( curr_lat-Y_edge(1) )
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

      DO Ki = 1,4,1

       ! check this grid box is the neast one or the one located in left-bottom
       ! ???
       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       if(i_lon==IIPAR+1) i_lon=1

       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2) 
       i_lev = Find_iPLev(curr_pressure,P_edge)

       !-----------------------------------------------------------------
       ! interpolate temperature for Plume module:
       !-----------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
          box_Ptemp(i_box) = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       else
          box_Ptemp(i_box) = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       endif

       !------------------------------------------------------------------
       ! For vertical wind speed:
       ! pay attention for the polar region * * *
       !------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       else
          curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       endif

       box_omeg(i_box) = curr_omeg

       ! For 
       dbox_lev = RK_Dt(Ki) * curr_omeg / 100.0     ! Pa => hPa
       curr_pressure = box_lev(i_box) + dbox_lev

       RK_Dlev(Ki) = RK_Dt(1) * curr_omeg / 100.0     ! Pa => hPa

       if(curr_pressure<P_mid(LLPAR)) &
             curr_pressure = P_mid(LLPAR) + ( P_mid(LLPAR) - curr_pressure )
       if(curr_pressure>P_mid(1)) &
             curr_pressure = P_mid(1) - ( curr_pressure - P_mid(1) )


       !------------------------------------------------------------------
       ! For the region where lat<72, use Regualr Longitude-Latitude Mesh:
       !------------------------------------------------------------------
       if(abs(curr_lat)<=72.0)then

           curr_u    = Interplt_wind_RLL(u, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
           curr_v    = Interplt_wind_RLL(v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

           box_u(i_box) = curr_u
           box_v(i_box) = curr_v

           dbox_lon  = (RK_Dt(Ki)*curr_u) &
                        / (2.0*PI*Re*cos(box_lat(i_box)*PI/180.0)) * 360.0
           dbox_lat  = (RK_Dt(Ki)*curr_v) / (PI*Re) * 180.0

           curr_lon  = box_lon(i_box) + dbox_lon
           curr_lat  = box_lat(i_box) + dbox_lat
 
           RK_Dlon(Ki) = (RK_Dt(1)*curr_u) &
                        / (2.0*PI*Re*cos(box_lat(i_box)*PI/180.0)) * 360.0
           RK_Dlat(Ki) = (RK_Dt(1)*curr_v) / (PI*Re) * 180.0

       endif

       !------------------------------------------------------------------
       ! For the polar region (lat>=72), use polar sterographic
       !------------------------------------------------------------------
       if(abs(curr_lat)>72.0)then

         if(abs(curr_lat)>Y_mid(JJPAR))then 
            curr_u_PS = Interplt_uv_PS_polar(1, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
            curr_v_PS = Interplt_uv_PS_polar(0, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         else
            curr_u_PS = Interplt_uv_PS(1, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
            curr_v_PS = Interplt_uv_PS(0, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
         endif

         box_u(i_box) = curr_u_PS
         box_v(i_box) = curr_v_PS

         dbox_x_PS = RK_Dt(Ki)*curr_u_PS
         dbox_y_PS = RK_Dt(Ki)*curr_v_PS

         RK_Dx_PS = RK_Dt(1)*curr_u_PS
         RK_Dy_PS = RK_Dt(1)*curr_v_PS


         !------------------------------------------------------------------
         ! change from (lon,lat) in RLL to (x,y) in PS: 
         !------------------------------------------------------------------
         if(box_lat(i_box)<0)then
           box_x_PS = -1.0* Re* cos(box_lon(i_box)*PI/180.0) &
                                / tan(box_lat(i_box)*PI/180.0)
           box_y_PS = -1.0* Re* sin(box_lon(i_box)*PI/180.0) &
                                / tan(box_lat(i_box)*PI/180.0)
         else
           box_x_PS = Re* cos(box_lon(i_box)*PI/180.0) &
                        / tan(box_lat(i_box)*PI/180.0)
           box_y_PS = Re* sin(box_lon(i_box)*PI/180.0) &
                        / tan(box_lat(i_box)*PI/180.0)
         endif

         RK_x_PS  = box_x_PS + RK_Dx_PS
         RK_y_PS  = box_y_PS + RK_Dy_PS

         box_x_PS  = box_x_PS + dbox_x_PS
         box_y_PS  = box_y_PS + dbox_y_PS


         !------------------------------------------------------------------
         ! change from (x,y) in PS to (lon,lat) in RLL
         !------------------------------------------------------------------
         if(box_x_PS>0.0)then
           curr_lon = atan( box_y_PS / box_x_PS )*180.0/PI 
         endif
         if(box_x_PS<0.0 .and. box_y_PS<=0.0)then
           curr_lon = atan( box_y_PS / box_x_PS )*180.0/PI -180.0
         endif
         if(box_x_PS<0.0 .and. box_y_PS>0.0)then
           curr_lon = atan( box_y_PS / box_x_PS )*180.0/PI +180.0
         endif
           
         if(curr_lat<0.0)then
           curr_lat= -1* atan( Re/sqrt(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         else
           curr_lat= atan( Re / sqrt(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         endif

         !------------------------------------------------------------------
         ! For 4th order Runge Kutta
         !------------------------------------------------------------------
         if(RK_x_PS>0.0)then
           RK_lon = atan( RK_y_PS / RK_x_PS )*180.0/PI
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS<=0.0)then
           RK_lon = atan( RK_y_PS / RK_x_PS )*180.0/PI -180.0
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS>0.0)then
           RK_lon = atan( RK_y_PS / RK_x_PS )*180.0/PI +180.0
         endif

         if(RK_lat<0.0)then
           RK_lat = -1 * atan( Re / sqrt(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         else
           RK_lat = atan( Re / sqrt(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         endif

         RK_Dlon(Ki) = RK_lon - box_lon(i_box)
         RK_Dlat(Ki) = RK_lat - box_lat(i_box)

       endif

      ENDDO ! Ki = 1,4,1

      box_lon(i_box) = box_lon(i_box) + &
                (RK_Dlon(1)+2.0*RK_Dlon(2)+2.0*RK_Dlon(3)+RK_Dlon(4))/6.0
      box_lat(i_box) = box_lat(i_box) + &
                (RK_Dlat(1)+2.0*RK_Dlat(2)+2.0*RK_Dlat(3)+RK_Dlat(4))/6.0
      box_lev(i_box) = box_lev(i_box) + &
                (RK_Dlev(1)+2.0*RK_Dlev(2)+2.0*RK_Dlev(3)+RK_Dlev(4))/6.0


      ! ====================================================================
      ! Add concentraion of PASV into conventional Eulerian GEOS-Chem in 
      ! corresponding with injected parcels in Lagrangian model
      ! For conventional GEOS-Chem for comparison with Lagrangian Model:
      !
      !  AD(I,J,L) = grid box dry air mass [kg]
      !  AIRMW     = dry air molecular wt [g/mol]
      !  MW_G(N)   = species molecular wt [g/mol]
      !     
      ! the conversion is:
      ! 
      !====================================================================
      if(i_box>N_prev)then
         nAdv = State_Chm%nAdvect         ! the last one is PASV_EU
         PASV_EU => State_Chm%Species(i_lon,i_lat,i_lev,nAdv)     ! Unit: [kg/kg]

         MW_g = State_Chm%SpcData(nAdv)%Info%emMW_g
         ! Here assume the injection rate is 30 kg/km for H2SO4: 
         PASV_EU = PASV_EU + box_length(i_box)*1.0e-3_fp*30.0 /State_Met%AD(i_lon,i_lat,i_lev)

         ! write(6,*)'== test 1 ==>', State_Chm%Spc_Units
      endif


    end do  !do i_box = 1,n_boxes_max

    !------------------------------------------------------------------
    ! Adjust the length/radius of the box based on new location
    !------------------------------------------------------------------
    i_box = 1
    lon1 = box_lon(i_box) - 0.5 * ( box_lon(i_box+1) - box_lon(i_box) )
    lon2 = box_lon(i_box) + 0.5 * ( box_lon(i_box+1) - box_lon(i_box) )
    lat1 = box_lat(i_box) - 0.5 * ( box_lat(i_box+1) - box_lat(i_box) )
    lat2 = box_lat(i_box) + 0.5 * ( box_lat(i_box+1) - box_lat(i_box) )

    length0           = box_length(i_box)
    box_length(i_box) = Distance_Circle(lon1, lat1, lon2, lat2)

    box_radiusA(i_box,:) = box_radiusA(i_box,:)*SQRT(length0/box_length(i_box))
    box_radiusB(i_box,:) = box_radiusB(i_box,:)*SQRT(length0/box_length(i_box))

    DO i_box = 2, n_boxes_max-1
      lon1 = 0.5 * ( box_lon(i_box) + box_lon(i_box-1) )
      lon2 = 0.5 * ( box_lon(i_box) + box_lon(i_box+1) )
      lat1 = 0.5 * ( box_lat(i_box) + box_lat(i_box-1) )
      lat2 = 0.5 * ( box_lat(i_box) + box_lat(i_box+1) )

      length0            = box_length(i_box)
      box_length(i_box)  = Distance_Circle(lon1, lat1, lon2, lat2)  ! [m]

      box_radiusA(i_box,:) = box_radiusA(i_box,:)*SQRT(length0/box_length(i_box))
      box_radiusB(i_box,:) = box_radiusB(i_box,:)*SQRT(length0/box_length(i_box))
    ENDDO

    i_box = n_boxes_max
    lon1 = box_lon(i_box) - 0.5 * ( box_lon(i_box) - box_lon(i_box-1) )
    lon2 = box_lon(i_box) + 0.5 * ( box_lon(i_box) - box_lon(i_box-1) )
    lat1 = box_lat(i_box) - 0.5 * ( box_lat(i_box) - box_lat(i_box-1) )
    lat2 = box_lat(i_box) + 0.5 * ( box_lat(i_box) - box_lat(i_box-1) )

    length0           = box_length(i_box)
    box_length(i_box) = Distance_Circle(lon1, lat1, lon2, lat2)

    box_radiusA(i_box,:) = box_radiusA(i_box,:)*SQRT(length0/box_length(i_box))
    box_radiusB(i_box,:) = box_radiusB(i_box,:)*SQRT(length0/box_length(i_box))

    !------------------------------------------------------------------
    ! Everything is done, clean up pointers
    !------------------------------------------------------------------
    nullify(u)
    nullify(v)
    nullify(omeg)
    nullify(X_edge)
    nullify(Y_edge)

  END SUBROUTINE lagrange_run


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real function Interplt_wind_RLL(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          ::  distance(2,2), Weight(2,2)
    real(fp)          ::wind_lonlat(2), wind_lonlat_lev

    ! Identify wether particle is exactly located on the grid point
    if(curr_pressure==P_mid(i_lev))then
    if(curr_lon==X_mid(i_lon))then
    if(curr_lat==Y_mid(i_lat))then

          Interplt_wind_RLL = wind(i_lon, i_lat, i_lev)
          return

    endif
    endif
    endif


    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
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

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo

    ! Calculate the inverse distance weight

    do i = 1,2
    do j = 1,2

       if(distance(i,j)==0)then
          Weight(:,:) = 0
          Weight(i,j) = 1
          GOTO 100
       endif

          Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )

    enddo
    enddo

 100 CONTINUE


    do k = 1,2     
        kk = k + init_lev - 1 
        if(init_lon==0)then
            wind_lonlat(k) =  Weight(1,1) * wind(IIPAR,init_lat,kk) &
                          + Weight(1,2) * wind(IIPAR,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        else if(init_lon==IIPAR)then
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk) &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                          + Weight(2,1) * wind(1,init_lat,kk)   &
                          + Weight(2,2) * wind(1,init_lat+1,kk)
        else
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk) &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        endif
    enddo


    ! second interpolate vertically (Linear)

    wind_lonlat_lev = wind_lonlat(1) + (wind_lonlat(2)-wind_lonlat(1)) &
                                 / (P_mid(init_lev+1)-P_mid(init_lev)) &
                                     * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev), P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL = wind_lonlat_lev

    return
  end function


! functions to interpolate vertical wind speed (w)
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The w value at polar point is the average of all surrounding grid points.

  real function Interplt_wind_RLL_polar(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(3), Weight(3)
    real(fp)          :: wind_lonlat(2), wind_lonlat_lev
    real(fp)          :: wind_polar

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
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
    
    if(init_lat==0)then
      init_lat = 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! calculate the distance between particle and grid point
    j = 1
    do i = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(1), Y_mid(jj))
      else
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo


   if(ii==0)then
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), 90.0e+0_fp)
   else if(ii==(IIPAR+1))then
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(1), 90.0e+0_fp)
   else
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), 90.0e+0_fp)
   endif


   IF(distance(3)==0.0)THEN
       do k=1,2
         kk = k + init_lev - 1
         wind_lonlat(k) = SUM(wind(:,init_lat,kk))/IIPAR
       enddo
   ELSE
       ! Calculate the inverse distance weight
       do i=1,3
          Weight(i) = 1.0/distance(i) / sum( 1.0/distance(:) )
       enddo

       do k=1,2
         kk = k + init_lev - 1

         wind_polar = SUM(wind(:,init_lat,kk))/IIPAR      

         if(init_lon==0)then    
             wind_lonlat(k) =  Weight(1)*wind(IIPAR,init_lat,kk)   &
                             + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         else if(init_lon==IIPAR)then
             wind_lonlat(k) =  Weight(1)*wind(init_lon,init_lat,kk)   &
                             + Weight(2)*wind(1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         else
             wind_lonlat(k) =  Weight(1)*wind(init_lon,init_lat,kk)   &
                             + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         endif
     enddo
   ENDIF

    ! second interpolate vertically (Linear)
    wind_lonlat_lev =   wind_lonlat(1) &
                      + (wind_lonlat(2)-wind_lonlat(1)) &
                      / (P_mid(init_lev+1)-P_mid(init_lev)) &
                      * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev),
    !P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL_polar = wind_lonlat_lev

    return
  end function

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v) 
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The u/v value at polar point is the average of all surrounding grid points.

  real function Interplt_uv_PS_polar(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)

    integer           :: i_uv
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk

    real(fp)          :: x_PS(3), y_PS(3)  ! the third value x_PS(3) is the polar point
    real(fp)          :: uv_PS(3,2)
    real(fp)          :: uv_polars(IIPAR)
    real(fp)          :: distance_PS(3), Weight_PS(3)
    real(fp)          :: uv_xy(2), uv_xy_lev

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
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

    ! For pressure level, P_mid(1) is about surface pressure, has biggerst
    ! value.
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! change from (lon,lat) in RLL to (x,y) in PS: 
    if(curr_lat<0)then
      curr_x = -1.0* Re* cos(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
      curr_y = -1.0* Re* sin(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
    else
      curr_x = Re* cos(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
      curr_y = Re* sin(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
    endif

    ! i,j means the four grid point value that surround around the particle
!    do i=1,2
!    do j=1,2

!      ii = i + init_lon - 1
!      jj = j + init_lat - 1
!
      ! For lon=180 deg:
!      if(ii==IIPAR+1)then
!        ii = 1
!      endif
!      if(ii==0)then
!        ii = IIPAR
!      endif

      ! Get the right ii and jj for interpolation at polar point:

      ! For South Polar Point:
 !     if(jj==0)then
 !       jj = jj+1
 !     endif

      ! For North Polar Point:
!      if(jj==JJPAR+1)then
!        jj = jj-1
!      endif

    if(init_lat==0)then
           jj = init_lat + 1    ! south polar
       else
           jj = init_lat
    endif
    
    do i=1,2  ! Interpolate location and wind of grid points into Polar Stereographic Plane
   
       ii = i + init_lon - 1

       ! For lon=180 deg:
       if(ii==IIPAR+1)then
          ii = 1
       endif
       if(ii==0)then
          ii = IIPAR
       endif

      ! Interpolate location and wind into Polar Stereographic Plane
      if(Y_mid(jj)>0)then
        x_PS(i)= Re* cos(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
        y_PS(i)= Re* sin(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
      else
        x_PS(i)= -1.0* Re* cos(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
        y_PS(i)= -1.0* Re* sin(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
      endif


       do k=1,2
          kk = k + init_lev - 1
          IF(i_uv==1)THEN ! i_ux==1 for u
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = -1.0* ( u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_PS(i,k) = u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                        + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF

          IF(i_uv==0)THEN ! for v
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                       - v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
          else
            uv_PS(i,k) = -1* u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                           + v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF
       enddo
    enddo

    ! Third grid point: south/north polar point
    x_PS(3) = 0.0
    y_PS(3) = 0.0
    
    do k=1,2
      kk = k + init_lev - 1
      IF(i_uv==1)THEN ! i_ux==1 for u
      ! interpolate all the grid points surrounding the polar point:
      do ii = 1,IIPAR
        if(Y_mid(jj)>0)then
          uv_polars(ii) = -1.0* ( u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                         / sin(Y_mid(jj)*PI/180.0) &
                         + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                         / (sin(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_polars(ii) = u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
          endif
        enddo

            uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF

          IF(i_uv==0)THEN ! for v
          do ii = 1,IIPAR
             if(Y_mid(jj)>0)then
               uv_polars(ii) = u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                             - v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
             else
               uv_polars(ii) = -1* u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                                 + v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
             endif
          enddo
             uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF
    enddo



    ! calculate the distance between particle and grid point
    do i = 1,3
       distance_PS(i)= sqrt( (x_PS(i)-curr_x)**2.0 + (y_PS(i)-curr_y)**2.0 )
    enddo

    ! Calculate the inverse distance weight
    do i = 1,3
        Weight_PS(i) = 1.0/distance_PS(i) / sum( 1.0/distance_PS(:) )
    enddo


    do k = 1,2
      uv_xy(k) = Weight_PS(1) * uv_PS(1,k) &
               + Weight_PS(2) * uv_PS(2,k) &
               + Weight_PS(3) * uv_PS(3,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ (uv_xy(2)-uv_xy(1)) &
     / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS_polar = uv_xy_lev

    return
  end function


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real function Interplt_uv_PS(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)

    integer           :: i_uv
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk

    real(fp)          :: x_PS(2,2), y_PS(2,2)
    real(fp)          :: uv_PS(2,2,2)
    real(fp)          :: distance_PS(2,2), Weight_PS(2,2)
    real(fp)          :: uv_xy(2), uv_xy_lev


    ! first interpolate horizontally (Inverse Distance Weighting)

 ! identify the grid point located in the southwest of the particle or under
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

!For pressure level, P_mid(1) is about surface pressure, has biggerst value.
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1

    
    ! change from (lon,lat) in RLL to (x,y) in PS: 
    if(curr_lat<0)then
      curr_x = -1.0* Re* cos(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
      curr_y = -1.0* Re* sin(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
    else
      curr_x = Re* cos(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
      curr_y = Re* sin(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
    endif

    ! i,j means the four grid point value that surround around the particle
    do i=1,2
    do j=1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! Get the right ii and jj for interpolation at polar point:
      ! For South Polar Point:
      if(jj==0)then
        jj = jj+1
      endif
      ! For North Polar Point:
      if(jj==JJPAR+1)then
        jj = jj-1
      endif

    
      ! For lon=180 deg:
      if(ii==IIPAR+1)then
        ii = 1
      endif
      if(ii==0)then
        ii = IIPAR
      endif


      ! Interpolate location and wind into Polar Stereographic Plane
      if(Y_mid(jj)>0)then

        x_PS(i,j) = Re* cos(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)  
        y_PS(i,j) = Re* sin(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
          
        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then ! i_ux==1 for u
             uv_PS(i,j,k)= -1.0* ( u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                         / sin(Y_mid(jj)*PI/180.0) &
                                + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2) )
           endif

           if(i_uv==0)then ! for v
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                          - v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      else

      x_PS(i,j)= -1.0* Re* cos(X_mid(ii)*PI/180.0) /tan(Y_mid(jj)*PI/180.0)
      y_PS(i,j)= -1.0* Re* sin(X_mid(ii)*PI/180.0) /tan(Y_mid(jj)*PI/180.0)

        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
           endif

           if(i_uv==0)then
             uv_PS(i,j,k) = -1* u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) &
                                        / sin(Y_mid(jj)*PI/180.0) &
                              + v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) &
                                        / (sin(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      endif

    enddo
    enddo

    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
       distance_PS(i,j) = sqrt( (x_PS(i,j)-curr_x)**2.0 &
                                                + (y_PS(i,j)-curr_y)**2.0 )
    enddo
    enddo

    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2
        Weight_PS(i,j) = 1.0/distance_PS(i,j) / sum( 1.0/distance_PS(:,:) )
    enddo
    enddo


    do k = 1,2
      uv_xy(k) =  Weight_PS(1,1) * uv_PS(1,1,k) &
                 + Weight_PS(1,2) * uv_PS(1,2,k) &
                 + Weight_PS(2,1) * uv_PS(2,1,k) &
                 + Weight_PS(2,2) * uv_PS(2,2,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ &
              (uv_xy(2)-uv_xy(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) &
                                        * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS = uv_xy_lev

    return
  end function


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
    integer  :: N_box

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

    real(fp)          :: Dx, Dy
    real(fp)          :: wind_s_shear
    real(fp)          :: theta_previous
    real(fp)          :: dbox_theta, dbox_radiusA, dbox_radiusB
    ! angle between travel direction and lon (-90 ~ 90 deg)
    real(fp)          :: box_alpha

    real(fp), pointer :: X_edge(:)
    real(fp), pointer :: Y_edge(:)

    real(fp), pointer :: P_I, R_e

    real(fp)  :: Outer(n_rings_max), Inner(n_rings_max)

    real(fp)  :: D_concnt
    ! used to calculate concentration distribution

    real(fp)  :: eddy_v, eddy_h
    real(fp)  :: eddy_A, eddy_B
    real(fp)  :: Cv, Ch, Omega_N, N_BV

    real(fp)  :: grid_volumn, exchange_amount

    real(fp) ::L_b, L_O, Ee
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

    Ptemp => State_Met%THETA ! ! Potential temperature [K]

    P_BXHEIGHT => State_Met%BXHEIGHT  ![IIPAR,JJPAR,KKPAR]

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!
    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)

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


    if(N_curr<=n_boxes_max)then
      N_box = N_curr
    else
      N_box = n_boxes_max
    endif

    !======================================================================
    ! Set up the new injected plumes/boxes
    ! identify the location of the plume box, 
    ! put the corresponding concentration of backgound grid into the initial 
    ! concentration inside plume in unit of [molec/cm3].
    !======================================================================
    IF(N_prev<n_boxes_max)THEN
    DO i_box = N_prev+1, N_box, 1       ! Only for new injected box

       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)   ! hPa

       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(curr_pressure,P_edge)

       do i_species = 1,N_species
        DO i_ring = 2,n_rings_max ! innest ring contains injectred aerosols
          box_concnt(i_box,i_ring,i_species) = &
                             State_Chm%Species(i_lon,i_lat,i_lev,i_species)
        ENDDO
          Extra_amount(i_box,i_species) = SUM(V_ring(i_box,2:n_rings_max)  &
                      * box_concnt(i_box,2:n_rings_max,i_species) ) ![molec]
       enddo 

    ENDDO
    ENDIF

    !=====================================================================
    ! Run Plume distortion & dilution HERE
    !=====================================================================
    do i_box = 1,N_box

       IF(Max_rings(i_box)==0) GOTO 200

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

       !------------------------------------------------------------------
       ! calcualte the box_alpha
       ! Next adjacent point
       !------------------------------------------------------------------
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

!       backgrd_concnt(i_box,:) = State_Chm%Species(i_lon,i_lat,i_lev,:)
       backgrd_concnt(i_box,1) = State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
       ! [molec/cm3]

       WRITE(6,*)'shw test:'
       WRITE(6,*)State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
       WRITE(6,*)State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect)

       Plume_I(i_box) = i_lon
       Plume_J(i_box) = i_lat
       Plume_L(i_box) = i_lev

       curr_u     = box_u(i_box)
       curr_v     = box_v(i_box)
       curr_omeg  = box_omeg(i_box)
       curr_Ptemp = box_Ptemp(i_box)


       !====================================================================
       ! For deformation of cross-section caused by wind shear 
       ! (A.D.Naiman et al., 2010):
       !====================================================================

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! *** attention *** ???
!       wind_s_shear = 0.0


       theta_previous   = box_theta(i_box)
       box_theta(i_box) = ATAN( TAN(box_theta(i_box)) + wind_s_shear*Dt )


       do i_ring=1,Max_rings(i_box)
       ! make sure use box_theta or TAN(box_theta)  ??? 
       box_radiusA(i_box,i_ring) = box_radiusA(i_box,i_ring) &
                                  * (TAN(box_theta(i_box))**2+1)**0.5 &
                                  / (TAN(theta_previous)**2+1)**0.5
       box_radiusB(i_box,i_ring) = box_radiusB(i_box,i_ring) &
                                  * (TAN(box_theta(i_box))**2+1)**(-0.5) &
                                  / (TAN(theta_previous)**2+1)**(-0.5)
       enddo

       !====================================================================
       ! For the concentration change caused by eddy diffusion:
       ! box_concnt(n_boxes_max,N_rings,N_species)
       !====================================================================

       ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
       Cv = 0.2
       Omega_N = 0.1
       Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)

       N_BV = sqrt(Ptemp_shear*g0/curr_Ptemp)
       IF(N_BV<=0.001) N_BV = 0.001

       eddy_v = Cv * Omega_N**2 / N_BV
!       eddy_v = 1.0

       ! Calculate horizontal eddy diffusivity:
!       Ch = 0.1
!       U_shear = Vertical_shear(u, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
!       V_shear = Vertical_shear(v, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
!       UV_shear = sqrt( U_shear**2 + V_shear**2 )
!       eddy_h(i_ring) = Ch*UV_shear*(Init_radius+(i_ring-1)*D_radius)**2
       eddy_h = 10.0  ! ???

! test for horizontal diffusivity
!       L_b = 2.0*PI*sqrt(curr_u**2+curr_v**2)/N_BV

!       Ee = 0.5 * Omega_N**2 * wind_s_shear**2

!       L_O = 2.0*PI*sqrt( Ee/(N_BV**3) )

!       WRITE(6,*)"= L_b, L_O =>", L_b, L_O
!       WRITE(6,*)"= Dh =>", eddy_h


       eddy_A = eddy_v*cos(box_theta(i_box)) &
                          + eddy_h*sin(abs(box_theta(i_box))) ! a
       eddy_B = eddy_v*sin(abs(box_theta(i_box))) &
                          + eddy_h*cos(box_theta(i_box)) ! b

       !=================================================================
       ! Calculate the transport rate
       ! kB should be rewrite in a more accurate equation !!!
       !=================================================================
       ! For innest ring:
       kB(1) = eddy_B / ( 0.5*box_radiusB(i_box,2) )
       kA(1) = eddy_A / ( 0.5*box_radiusA(i_box,2) )


       do i_ring = 2, Max_rings(i_box)-1

         kB(i_ring) = eddy_B &
           /( 0.5*(box_radiusB(i_box,i_ring+1)-box_radiusB(i_box,i_ring-1)) )
         kA(i_ring) = eddy_A &
           /( 0.5*(box_radiusA(i_box,i_ring+1)-box_radiusA(i_box,i_ring-1)) )

       enddo ! 

       ! For outest ring (i_ring = n_rings_max)
       kB(Max_rings(i_box)) = eddy_B &
        /( box_radiusB(i_box,Max_rings(i_box))-box_radiusB(i_box,Max_rings(i_box)-1) )
       kA(Max_rings(i_box)) = eddy_A &
        /( box_radiusA(i_box,Max_rings(i_box))-box_radiusA(i_box,Max_rings(i_box)-1) )


       !==================================================================
       ! Begin to calculate the dilution inside plume 
       !==================================================================

       !---------------------------------------------------------------------
       ! Decide the time step based on CFL condition
       ! 2*k*Dt/Dr<1 (or Dr-2*k*Dt>0) for diffusion
       !---------------------------------------------------------------------
       Dt2 = Dt

       IF(MINVAL(box_radiusA(i_box,:)-2.0*kA(:)*Dt2)<0.0 &
                        .or. MINVAL(box_radiusB(i_box,:)-2.0*kB(:)*Dt2)<0.0 )THEN
          Dt2 = Dt*0.1
       ENDIF


       if(MINVAL(box_radiusA(i_box,:)-2.0*kA(:)*Dt2)<0.0 &
                        .or. MINVAL(box_radiusB(i_box,:)-2.0*kB(:)*Dt2)<0.0 )then
         Dt2 = Dt*0.01
         WRITE(6,*)"Plume Reporting: it's time to merge rings!"

 300     CONTINUE
         !-------------------------------------------------------------------
         ! Decrease half of rings by combining 2 rings into 1 ring
         !-------------------------------------------------------------------
         DO i_ring=1,Max_rings(i_box)/2
           box_radiusA(i_box,i_ring) = box_radiusA(i_box,i_ring*2)
           box_radiusB(i_box,i_ring) = box_radiusB(i_box,i_ring*2)

           box_concnt(i_box,i_ring,N_species) = &
                   (  box_concnt(i_box,i_ring*2-1,N_species)*V_ring(i_box,i_ring*2-1) &
                     +box_concnt(i_box,i_ring*2,N_species)*V_ring(i_box,i_ring*2) ) &
                 / ( V_ring(i_box,i_ring*2-1)+V_ring(i_box,i_ring*2) )

           V_ring(i_box,i_ring) =V_ring(i_box,i_ring*2-1)+V_ring(i_box,i_ring*2)
         ENDDO

         Max_rings(i_box) = Max_rings(i_box)/2

         !=================================================================
         ! Calculate the transport rate
         ! kB should be rewrite in a more accurate equation !!!
         !=================================================================
         ! For innest ring:
         kB(1) = eddy_B / ( 0.5*box_radiusB(i_box,2) )
         kA(1) = eddy_A / ( 0.5*box_radiusA(i_box,2) )

         do i_ring = 2, Max_rings(i_box)-1
           kB(i_ring) = eddy_B &
             /( 0.5*(box_radiusB(i_box,i_ring+1)-box_radiusB(i_box,i_ring-1)) )
           kA(i_ring) = eddy_A &
             /( 0.5*(box_radiusA(i_box,i_ring+1)-box_radiusA(i_box,i_ring-1)) )
         enddo ! 

         ! For outest ring (i_ring = n_rings_max)
         kB(Max_rings(i_box)) = eddy_B /( box_radiusB(i_box,Max_rings(i_box)) &
                                         -box_radiusB(i_box,Max_rings(i_box)-1) )
         kA(Max_rings(i_box)) = eddy_A /( box_radiusA(i_box,Max_rings(i_box)) &
                                         -box_radiusA(i_box,Max_rings(i_box)-1) )

       endif

       !===================================================================
       ! Only 1 ring left for the whole plume,
       ! dissolve the plume into the background grid cell
       !===================================================================
       IF(Max_rings(i_box)==1)THEN
         backgrd_concnt(i_box,1) = &
            ( SUM(box_concnt( i_box, 1:Max_rings(i_box), 1 )*V_ring(i_box,1:Max_rings(i_box) )) &
             +backgrd_concnt(i_box,1)*grid_volumn - Extra_amount(i_box,1) ) &
           / grid_volumn

         State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = backgrd_concnt(i_box,1)

         Max_rings(i_box) = 0
         WRITE(6,*)'Plume Report: Only 1 ring left in plume!'
         GOTO 200 ! skip this box, go to next box
       ENDIF


       IF(MINVAL(box_radiusA(i_box,:)-2.0*kA(:)*Dt2)<0.0 &
                        .or. MINVAL(box_radiusB(i_box,:)-2.0*kB(:)*Dt2)<0.0 )THEN
         WRITE(6,*)"Plume Reporting: Still need to merge rings!" 
         GOTO 300
       ENDIF


!       DO i_species = 1,N_species
       DO i_species = 1, 1
       do t1s=1,int(Dt/Dt2)

       !---------------------------------------------------------------------
       ! Use Runge-Kutta method (RK4) to solve diferential equation
       !---------------------------------------------------------------------
       do Ki = 1,4

         DO i_ring = 1, Max_rings(i_box)
           if(Ki==1)then
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species)
           else if(Ki==4)then
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species) &
                                  + RK(3,i_ring) !Dt
           else
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring,i_species) &
                                  + 0.5*RK(Ki-1,i_ring) !Dt ???
           endif
         ENDDO ! i_ring

! ??? the 4th order runge kutta method used here is not correct, 
! should use 0.5*Dt instead of 0.5*RK(Ki-1,i_ring)

         ! For innest ring:
           D_concnt = box_concnt_K(2) - box_concnt_K(1)
         Outer(1) = Amount_Dilute(D_concnt, &
                         box_radiusA(i_box,1), box_radiusB(i_box,1), &
                                 kA(1), kB(1), Dt2, box_length(i_box))

         Inner(1) = 0.0

         RK(Ki,1) = ( Outer(1)+Inner(1) ) / V_ring(i_box,1) ! [molec/cm3]

         ! For rings from 2 to (n_rings_max - 1)
         do i_ring = 2, Max_rings(i_box)-1
           D_concnt      = box_concnt_K(i_ring+1)-box_concnt_K(i_ring)
           Outer(i_ring) = Amount_Dilute(D_concnt, &
                         box_radiusA(i_box,i_ring), box_radiusB(i_box,i_ring), &
                                 kA(i_ring), kB(i_ring), Dt2, box_length(i_box))

           D_concnt      = box_concnt_K(i_ring-1)-box_concnt_K(i_ring)
           Inner(i_ring) = Amount_Dilute(D_concnt, &
                         box_radiusA(i_box,i_ring-1), box_radiusB(i_box,i_ring-1), &
                                 kA(i_ring-1), kB(i_ring-1), Dt2, box_length(i_box))

           RK(Ki,i_ring) = ( Outer(i_ring)+Inner(i_ring) ) / V_ring(i_box,i_ring)
         enddo ! i_ring

         ! For outest ring:
         D_concnt = backgrd_concnt(i_box,i_species) - box_concnt_K(Max_rings(i_box))
         Outer(Max_rings(i_box)) = Amount_Dilute(D_concnt , &
              box_radiusA(i_box,Max_rings(i_box)), box_radiusB(i_box,Max_rings(i_box)), &
                       kA(Max_rings(i_box)), kB(Max_rings(i_box)), Dt2, box_length(i_box))

         D_concnt = box_concnt_K(Max_rings(i_box)-1)-box_concnt_K(Max_rings(i_box))
         Inner(Max_rings(i_box)) = Amount_Dilute(D_concnt , &
              box_radiusA(i_box,Max_rings(i_box)), box_radiusB(i_box,Max_rings(i_box)), &
                       kA(Max_rings(i_box)), kB(Max_rings(i_box)), Dt2, box_length(i_box))

         RK(Ki,Max_rings(i_box)) = ( Outer(Max_rings(i_box)) + Inner(Max_rings(i_box)) ) &
                              / V_ring(i_box,Max_rings(i_box))

         Outer2env(Ki) = -1.0 * Outer(Max_rings(i_box)) ! [molec]

       enddo ! Ki


       do i_ring = 1,Max_rings(i_box)
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

!       State_Chm%Species(i_lon,i_lat,i_lev,i_species) = backgrd_concnt(i_box,i_species)
       State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = backgrd_concnt(i_box,i_species)

       !===================================================================
       ! Once the concentration in different rings don't have big difference
       ! dissolve the plume into the background grid cell
       !===================================================================

       IF( box_concnt(i_box,1,1)<box_concnt(i_box,2,1) )THEN ! ??? this criteria should be changed

         backgrd_concnt(i_box,i_species) = &
            ( SUM(box_concnt( i_box, 1:Max_rings(i_box), 1 )*V_ring( i_box,1:Max_rings(i_box) )) &
             +backgrd_concnt(i_box,i_species)*grid_volumn - Extra_amount(i_box,1) ) / grid_volumn

         State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = backgrd_concnt(i_box,i_species) 

         Max_rings(i_box) = 0
         GOTO 200 ! skip this box, go to next box
       ENDIF


       enddo ! t1s
       ENDDO ! do i_Species = 1,nSpecies

 200 CONTINUE

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


    N_prev = N_curr
    N_curr = N_prev + N_parcels
    ! N_curr is used to add particles, here add 5 parcles every time step


    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)
    nullify(omeg)
    nullify(X_edge)
    nullify(Y_edge)

  END SUBROUTINE plume_run

!--------------------------------------------------------------------
! functions to which grid cell (i,j,k) contains the box
!--------------------------------------------------------------------

  integer function Find_iLonLat(curr_xy,  Dxy,  XY_edge2)
    implicit none
    real(fp) :: curr_xy, Dxy, XY_edge2
    Find_iLonLat = INT( (curr_xy - (XY_edge2 - Dxy)) / Dxy )+1
    ! Notice the difference between INT(), FLOOR(), AINT()
    ! for lon: Xedge_Sec - Dx = Xedge_first
    return
  end function


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


  !-------------------------------------------------------------------
  ! calculation the great-circle distance between two points on the earth
  ! surface
  !-------------------------------------------------------------------

  real function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2  ! unit is degree
    !real(fp) :: PI, Re

    Distance_Circle = Re * 2.0 * &
                ASIN( (sin( (y1-y2)*PI/180.0 ))**2.0 &
                       + cos(x1*PI/180.0) * cos(x2*PI/180.0) &
                         * (sin( 0.5*(x1-x2)*PI/180.0 ))**2.0 )
    return
  end function


  !-------------------------------------------------------------------  
  ! calculate the wind_s (inside a plume sross-section) shear along pressure
  ! direction
  !-------------------------------------------------------------------  

  real function Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: box_alpha
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: u(:,:,:), v(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
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

      u_lonlat(k) =  Weight(1,1) * u(i_lon,i_lat,kk) &
                   + Weight(1,2) * u(i_lon,i_lat+1,kk) &
                   + Weight(2,1) * u(i_lon+1,i_lat,kk) &
                   + Weight(2,2) * u(i_lon+1,i_lat+1,kk)

      v_lonlat(k) =  Weight(1,1) * v(i_lon,i_lat,kk) &
                    + Weight(1,2) * v(i_lon,i_lat+1,kk) &
                    + Weight(2,1) * v(i_lon+1,i_lat,kk) &
                    + Weight(2,2) * v(i_lon+1,i_lat+1,kk)

      if(box_alpha>=0.0)then
        wind_s(k) = -1.0 * u_lonlat(k) * SIN(box_alpha) &
                        + v_lonlat(k) * COS(box_alpha)
      else
        wind_s(k) = u_lonlat(k)*SIN(box_alpha) + v_lonlat(k)*COS(box_alpha)
      endif

    enddo


    ! second vertical shear of wind_s

    ! This code should be changed !!!
   ! Because it is the poressure center in [hPa] instead of height center in
    ! [m]
    ! Delt_height    = 0.5 * ( P_BXHEIGHT(init_lon,init_lat,init_lev) +
    ! P_BXHEIGHT(init_lon,init_lat,init_lev+1) )
    Delt_height =  Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), P_edge(init_lev), P_edge(init_lev+1), 1 )   &
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1), P_edge(init_lev), P_edge(init_lev+1), 0 )
    ! find the z height of each pressure level in GEOS-Chem

    Wind_shear_s = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function

!-------------------------------------------------------------------
! transform the pressure level [Pa] to the height level [m]
!-------------------------------------------------------------------

  real function Pa2meter(Box_height, P1, P2, Judge)
  ! Judge: 0 for bottom, 1 for top
  ! 
    implicit none
    real(fp) :: Box_height, P1, P2
    integer  :: Judge

      if(Judge==1)then
        ! calculate the height of the top half of the grid box
        Pa2meter = Box_height * ( DLOG(P2) - DLOG(0.5*(P2+P1)) ) &
                        / ( DLOG(P2) - DLOG(P1) )
      else
        ! calculate the height of the bottom half of the grid box
        Pa2meter = Box_height * ( DLOG(0.5*(P2+P1)) - DLOG(P1) ) &
                        / ( DLOG(P2) - DLOG(P1) )
      endif

    return
  end function


!-------------------------------------------------------------------
! Calculate eddy diffusivity in stratisphere. (U.Schumann, 2012)


  !-------------------------------------------------------------------  
  ! calculate the vertical shear for calculating eddy difussivity

  real function Vertical_shear(var, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: var(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
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

    Delt_height = &
      Pa2meter(P_BXHEIGHT(init_lon,init_lat,init_lev), P_edge(init_lev), P_edge(init_lev+1), 1 ) &
     +Pa2meter(P_BXHEIGHT(init_lon,init_lat,init_lev+1), P_edge(init_lev+1), P_edge(init_lev+2), 0)
    ! find the z height of each pressure level in GEOS-Chem

    Vertical_shear = ( var_lonlat(2) - var_lonlat(1) ) / Delt_height

    return
  end function


!===================================================================
!
!===================================================================

  REAL FUNCTION Amount_Dilute(D_concnt, Ra, Rb, kA, kB, Dt, length)

    IMPLICIT NONE
     
    REAL(fp)    :: D_concnt, Ra, Rb, kA, kB, length
    REAL        :: Dt

      Amount_Dilute = D_concnt * 1.0e+6_fp * length * PI &
                    * (  (Ra+0.5*kA*Dt) * (Rb+0.5*kB*Dt) &
                       - (Ra-0.5*kA*Dt) * (Rb-0.5*kB*Dt)  ) 

    return

  END FUNCTION


  SUBROUTINE lagrange_write_std( am_I_Root, RC )
    
    
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

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR
    USE TIME_MOD,        ONLY : GET_MINUTE
    USE TIME_MOD,        ONLY : GET_SECOND

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
    LOGICAL                   :: IsOldFile

    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND

    CHARACTER(LEN=255)            :: FILENAME
    CHARACTER(LEN=255)            :: FILENAME2

    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C

    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()
    MINUTE      = GET_MINUTE()
    SECOND      = GET_SECOND()

    WRITE(YEAR_C,*) YEAR
    WRITE(MONTH_C,*) MONTH
    WRITE(DAY_C,*) DAY
    WRITE(HOUR_C,*) HOUR
    WRITE(MINUTE_C,*) MINUTE
    WRITE(SECOND_C,*) SECOND

!-------------------------------------------------------------------
! Output box location
!-------------------------------------------------------------------

    FILENAME   = 'Lagrange_xyz_' // TRIM(ADJUSTL(YEAR_C)) // '-' &
        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'
!    FILENAME   = 'Lagrange_1day_box_i_lon_lat_lev.txt'

    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)
!       OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='OLD', &
!             FORM='UNFORMATTED', ACCESS='Direct', Recl=18 )
    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      Do i_box = 1, n_boxes_max
        WRITE(261,'(2(x,I8),3(x,E12.5))') N_prev, i_box, &
                 box_lon(i_box), box_lat(i_box), box_lev(i_box)
      End Do
    ENDIF


!-------------------------------------------------------------------
! Output plume location
!-------------------------------------------------------------------
!    IF(mod(tt,6)==0)THEN     ! output once every hour
!    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)

       FILENAME2   = 'Plume_concentration_molec_' // TRIM(ADJUSTL(YEAR_C)) // &
          '-' //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) //      &
          '-' //TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) //    &
          ':' //TRIM(ADJUSTL(SECOND_C)) // '.txt'

       OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
             FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       DO i_box = 1,1001,100
          WRITE(262,*) Max_rings(i_box)
          WRITE(262,*) box_concnt(i_box,:,N_species)*V_ring(i_box,:)
          WRITE(262,*) box_concnt(i_box,:,N_species)
       ENDDO

!       DO i_box = 1, n_boxes_max
!          WRITE(262,*) i_box, SUM(box_concnt(i_box,:,N_species)*V_ring(i_box,:))
!       ENDDO

!        WRITE(262,*) box_concnt(1,:,N_species)

!    ENDIF
!    ENDIF

    tt = tt + 1

  END SUBROUTINE lagrange_write_std

!-------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  subroutine lagrange_cleanup()

    if (allocated(box_lon))      deallocate(box_lon)
    if (allocated(box_lat))      deallocate(box_lat)
    if (allocated(box_lev))      deallocate(box_lev)
    if (allocated(box_length))   deallocate(box_length)
    if (allocated(box_radiusA))  deallocate(box_radiusA)
    if (allocated(box_radiusB))  deallocate(box_radiusB)
    if (allocated(box_theta))    deallocate(box_theta)

    WRITE(6,'(a)') '--> Lagrange and Plume Module Cleanup <--'


  end subroutine lagrange_cleanup


END MODULE Lagrange_Mod
