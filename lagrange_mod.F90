!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!-----------------------------------------------------------------------------!

MODULE Lagrange_Mod

  USE precision_mod
  USE PhysConstants, ONLY : PI, Re, AIRMW
  USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: lagrange_init
  PUBLIC :: lagrange_run
  PUBLIC :: lagrange_write_std
  PUBLIC :: lagrange_cleanup


  integer, parameter    :: n_boxes_max = 7000000 !  6904224           ! 24*200*1 : lat*lon*lev
  integer, parameter    :: N_parcels   = 131        
  integer               :: tt, N_Dt, N_Dt_previous         ! Aircraft would release 131 aerosol parcels every time step
  integer               :: i_rec

  real(fp), allocatable :: box_lon(:)    
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)
  real(fp), allocatable :: box_depth(:)
  real(fp), allocatable :: box_width(:)
  real(fp), allocatable :: box_length(:)


CONTAINS


!-----------------------------------------------------------------

  SUBROUTINE lagrange_init( am_I_root, State_Met )

    USE State_Met_Mod,   ONLY : MetState
    USE GC_GRID_MOD,     ONLY : GET_AREA_M2

    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met

    INTEGER                       :: i_box, iibox
    INTEGER                       :: ii, jj, kk
    CHARACTER(LEN=255)            :: FILENAME, FILENAME_INIT

    integer :: i_lon            !1:IIPAR
    integer :: i_lat
    integer :: i_lev

    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Lagrnage Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'

    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))
    allocate(box_width(n_boxes_max))
    allocate(box_depth(n_boxes_max))
    allocate(box_length(n_boxes_max))


!--------------------------------------------------

    ! (1)
    FILENAME_INIT   = '/n/home12/hongwei/hongwei/merra2_2x25_standard_Mar/Lagrange_Init_box_i_lon_lat_lev.txt'
    OPEN( 361,      FILE=TRIM( FILENAME_INIT   ), STATUS='old', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    Do i_box = 1, n_boxes_max
       READ(361,'( 2(x,I8) , 3(x,E12.5) )') N_Dt_previous, iibox, box_lon(i_box), box_lat(i_box), box_lev(i_box)
    End Do
    
!--------------------------------------------------

    ! (2)
!    do i_box = 1,n_boxes_max,1
!        box_lon(i_box) = -141.0      
!        box_lat(i_box) = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) * (-1.0)**FLOOR(i_box/6000.0)      ! -29.95S : 29.95N : 0.1
!        box_lev(i_box) = 52.0e+0_fp       ! about 20 km
!    enddo
!    N_Dt_previous = 0

!-------------------------------------------------

    box_width  = 0.0e+0_fp
    box_depth  = 0.0e+0_fp
    box_length = 0.0e+0_fp


    FILENAME   = 'Lagrange_1day_box_i_lon_lat_lev.txt'
    tt   = 0
    N_Dt = N_Dt_previous + N_parcels 

!    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE', &
!          FORM='UNFORMATTED',    ACCESS='Direct', Recl=18 )
    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
! Integer is 2 bytes, 2+4+4+4=18 bytes.


!    i_rec = 0
    Do i_box = 1, n_boxes_max
!       i_rec = i_rec + 1
       WRITE(261,'( 2(x,I8) , 3(x,E12.5) )') N_Dt_previous, i_box, box_lon(i_box), box_lat(i_box), box_lev(i_box)
!       WRITE(261,'( 1(x,I8) , 3(x,E12.5) )') i_box, box_lon(i_box), box_lat(i_box), box_lev(i_box)
    End Do


! Output State_Met%AD(i_lon,i_lat,i_lev) into State_Met_AD.txt
    OPEN( 314,      FILE='State_Met_AD.txt', STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    Do i_lon = 1, IIPAR
    Do i_lat = 1, JJPAR
    Do i_lev = 1, LLPAR
       WRITE(314,'(x,E12.5)') State_Met%AD(i_lon,i_lat,i_lev)
    End Do
    End Do
    End Do


! Output State_Met%AREA_M2(i_lon,i_lat,i_lev) into State_Met_AREA_M2.txt
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

  SUBROUTINE lagrange_run(am_I_Root, State_Chm, State_Met, Input_Opt)

    USE Input_Opt_Mod, ONLY : OptInput
    ! USE PhysConstants, ONLY : PI, Re 
    ! Re    : Radius of Earth [m]

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE TIME_MOD,      ONLY : GET_TS_DYN

    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE, XMID, YMID                 
    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON 
    ! USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON 
    ! DLAT( IIPAR, JJPAR, LLPAR ), DLON( IIPAR, JJPAR, LLPAR )
    ! XEDGE  ( IM+1, JM,   L ), YEDGE  ( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt

    REAL :: Dt                  ! = 600.0e+0_fp          

    real(fp), pointer  :: PASV           
    integer            :: nAdv        
    REAL(fp)           :: MW_g

    integer :: i_box, N_box

    integer :: i_lon            !1:IIPAR
    integer :: i_lat
    integer :: i_lev

    integer :: ii, jj, kk

    real(fp) :: curr_lon
    real(fp) :: curr_lat
    real(fp) :: curr_pressure
    real(fp) :: X_edge2
    real(fp) :: Y_edge2

    real(fp) :: dbox_lon, dbox_lat, dbox_lev
    real(fp) :: dbox_x_PS, dbox_y_PS
    real(fp) :: box_x_PS, box_y_PS

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: curr_u_PS, curr_v_PS

    real(fp), pointer :: P_edge(:)
    real(fp), pointer :: P_mid(:)

    real(fp) :: Dx
    real(fp) :: Dy
    real(fp), pointer :: X_edge(:)    
    real(fp), pointer :: Y_edge(:)       
    real(fp), pointer :: X_mid(:)    
    real(fp), pointer :: Y_mid(:)       

    real(fp), pointer :: P_I, R_e     

    Dt = GET_TS_DYN()

    ! Establish pointers
    ! P_I = PI
    ! R_e = Re

    u => State_Met%U   ! figure out state_met%U is based on lat/lon or modelgrid(i,j)
    v => State_Met%V   ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]
    P_mid => State_Met%PMID(1,1,:)  ! Pressure (w/r/t moist air) at level centers (hPa)

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!
    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)  
    X_mid => XMID(:,1,1)   ! IIPAR
    Y_mid => YMID(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)
     
    
    if(N_Dt<=n_boxes_max)then
      N_box = N_Dt
    else
      N_box = n_boxes_max
    endif

    ! Run Lagrangian advection HERE
    do i_box = 1,N_box,1

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


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2) 
       i_lev = Find_iPLev(curr_pressure,P_mid)


       if(i_lon==IIPAR+1) i_lon=1


       ! For vertical direction:
       ! pay attention for the polar region * * *
       if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_omeg = Interplt_wind_RLL_polar(omeg,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       else
          curr_omeg = Interplt_wind_RLL(omeg,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       endif

       ! ==================================================
       ! test !
!       curr_omeg = 0.0
       ! ==================================================

       dbox_lev = Dt * curr_omeg / 100.0     ! Pa => hPa
       box_lev(i_box) = box_lev(i_box) + dbox_lev

       if(box_lev(i_box)<P_mid(LLPAR)) box_lev(i_box) = P_mid(LLPAR) + ( P_mid(LLPAR) - box_lev(i_box) )
       if(box_lev(i_box)>P_mid(1)) box_lev(i_box) = P_mid(1) - ( box_lev(i_box) - P_mid(1) )

       ! test: replace the same wind around polar with adjacent different values =======================
!       do ii=1,IIPAR,1
!       do kk=1,LLPAR,1

!        u(ii,1,kk) = u(ii,2,kk)*0.6
!        v(ii,1,kk) = v(ii,2,kk)*0.6

!        u(ii,JJPAR,kk) = u(ii,JJPAR-1,kk)*0.6
!        v(ii,JJPAR,kk) = v(ii,JJPAR-1,kk)*0.6

!       enddo
!       enddo

       ! test for solid-body rotation ===================================================================
!       do ii=1,IIPAR,1
!       do jj=1,JJPAR,1
!       do kk=i_lev-2,i_lev+2,1
!         u(ii, jj, kk) = sin(Y_mid(jj)*PI/180.0) *( -1.0*sin(X_mid(ii)*PI/180.0)*0.0 +cos(X_mid(ii)*PI/180.0)*3.0 )
!         v(ii, jj, kk) = sin(Y_mid(jj)*PI/180.0)**2 *( -1.0*cos(X_mid(ii)*PI/180.0)*0.0 -sin(X_mid(ii)*PI/180.0)*3.0 )
!       enddo
!       enddo
!       enddo
       ! test:====================================================================


       if(abs(curr_lat)<=72.0)then
       ! Regualr Longitude-Latitude Mesh:

         curr_u    = Interplt_wind_RLL(u,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         curr_v    = Interplt_wind_RLL(v,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

         dbox_lon = (Dt*curr_u) / (2.0*PI*Re*cos(curr_lat*PI/180.0)) * 360.0
         dbox_lat = (Dt*curr_v) / (PI*Re) * 180.0

         box_lon(i_box) = box_lon(i_box) + dbox_lon
         box_lat(i_box) = box_lat(i_box) + dbox_lat
 
         ! write(6,*)'= RLL =>', i_box, dbox_lon, curr_lon, box_lon(i_box)
         ! write(6,*)'= RLL =>', dbox_lat, curr_lat, box_lat(i_box)
      endif

       if(abs(curr_lat)>72.0)then
       ! Polar Stereographic plane (Dong and Wang, 2012):

         if(abs(curr_lat)>Y_mid(JJPAR))then 
            curr_u_PS = Interplt_uv_PS_polar(1, u, v, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
            curr_v_PS = Interplt_uv_PS_polar(0, u, v, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         else
            curr_u_PS = Interplt_uv_PS(1, u, v, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
            curr_v_PS = Interplt_uv_PS(0, u, v, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
         endif

         dbox_x_PS = Dt*curr_u_PS
         dbox_y_PS = Dt*curr_v_PS

         ! change from (lon,lat) in RLL to (x,y) in PS: 
         if(curr_lat<0)then
           box_x_PS = -1.0* Re* cos(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
           box_y_PS = -1.0* Re* sin(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
         else
           box_x_PS = Re* cos(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
           box_y_PS = Re* sin(curr_lon*PI/180.0) / tan(curr_lat*PI/180.0)
         endif

         box_x_PS  = box_x_PS + dbox_x_PS
         box_y_PS  = box_y_PS + dbox_y_PS

         if(box_x_PS>0.0)then
           box_lon(i_box) = atan( box_y_PS / box_x_PS )*180.0/PI 
         endif
         if(box_x_PS<0.0 .and. box_y_PS<=0.0)then
           box_lon(i_box) = atan( box_y_PS / box_x_PS )*180.0/PI -180.0
         endif
         if(box_x_PS<0.0 .and. box_y_PS>0.0)then
           box_lon(i_box) = atan( box_y_PS / box_x_PS )*180.0/PI +180.0
         endif
           
         if(curr_lat<0.0)then
           box_lat(i_box) = -1 * atan( Re / sqrt(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         else
           box_lat(i_box) = atan( Re / sqrt(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         endif

         ! write(6,*)'= PS =>', i_box, curr_lon, box_lon(i_box), dbox_x_PS, box_x_PS
         ! write(6,*)'= PS =>', curr_lat, box_lat(i_box), dbox_x_PS, box_x_PS

       endif

       ! Add concentraion of PASV into conventional Eulerian GEOS-Chem in corresponding with injected parcels in Lagrangian model
       if(i_box>N_Dt_previous)then
          ! ============================================================================
          ! For conventional GEOS-Chem for comparison with Lagrangian Model:
          !
          !  AD(I,J,L) = grid box dry air mass [kg]
          !  AIRMW     = dry air molecular wt [g/mol]
          !  MW_G(N)   = species molecular wt [g/mol]
          !     
          ! the conversion is:
          ! 
          !====================================================================
          nAdv = State_Chm%nAdvect         ! the last one is PASV
          PASV => State_Chm%Species(i_lon,i_lat,i_lev,nAdv)     ! Unit: [kg/kg]
 
          MW_g = State_Chm%SpcData(nAdv)%Info%emMW_g
          ! Assume every parcle has 110kg H2SO4
          ! write(6,*)'== test 1 ==>', State_Chm%Spc_Units
          PASV = PASV + 110.0/State_Met%AD(i_lon,i_lat,i_lev)

          !WRITE(6,*) '= Lagrange for Eulerian =>', i_box
       endif


       !call
       !grow_box(box_length(i_box),box_width(i_box),box_depth(i_box),i_lon,i_lat,i_lev)
    end do  !do i_box = 1,n_boxes_max

        N_Dt_previous = N_Dt
        N_Dt = N_Dt_previous + N_parcels        ! N_Dt is used to add particles, here add 5 parcles every time step


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

  END SUBROUTINE lagrange_run

!--------------------------------------------------------------------
! functions to find the grid point (i,j,k) which is nearest to location of boxes 

  integer function Find_iLonLat(curr_xy,  Dxy,  XY_edge2)
    implicit none
    real(fp) :: curr_xy, Dxy, XY_edge2
    Find_iLonLat = INT( (curr_xy - (XY_edge2 - Dxy)) / Dxy )+1
    ! Notice the difference between INT(), FLOOR(), AINT()
    ! for lon: Xedge_Sec - Dx = Xedge_first
    return
  end function


  !integer function find_latitude(curr_lat, Dy, Y_edge2)
  !  implicit none
  !  real :: curr_lat, Dy, Y_edge2
  !  find_latitude = INT( (curr_lat - (Y_edge2 - Dy)) / Dy )+1
  !  ! for lat: (Yedge_Sec - Dy) may be different from Yedge_first
  !  ! region for 4*5: (-89,-86:4:86,89))
  !  return
  !end function


  integer function Find_iPLev(curr_pressure,P_mid)
    implicit none
    real(fp) :: curr_pressure
    real(fp), pointer :: P_mid(:)
    integer :: i_lon, i_lat
    integer :: locate(1)
    locate = MINLOC(abs( P_mid(:)-curr_pressure ))

    if(P_mid(locate(1))-curr_pressure >= 0 )then
       Find_iPLev = locate(1)
    else
       Find_iPLev = locate(1) - 1
    endif

    if(Find_iPLev==0) Find_iPLev=1
    if(Find_iPLev==LLPAR) Find_iPLev=LLPAR-1

    return
  end function

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real function Interplt_wind_RLL(wind, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
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
          distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
        else if(ii==(IIPAR+1))then
          distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
        else
          distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
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
            wind_lonlat(k) =  Weight(1,1) * wind(IIPAR,init_lat,kk)     &
                          + Weight(1,2) * wind(IIPAR,init_lat+1,kk)   &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk)   &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        else if(init_lon==IIPAR)then
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk)     &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk)   &
                          + Weight(2,1) * wind(1,init_lat,kk)   &
                          + Weight(2,2) * wind(1,init_lat+1,kk)
        else
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk)     &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk)   &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk)   &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        endif
    enddo


    ! second interpolate vertically (Linear)

    wind_lonlat_lev = wind_lonlat(1)    &
                     + (wind_lonlat(2)-wind_lonlat(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev), P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL = wind_lonlat_lev

    return
  end function


! functions to interpolate vertical wind speed (w)
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The w value at polar point is the average of all surrounding grid points.

  real function Interplt_wind_RLL_polar(wind, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)
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
          distance(i) = Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
        else if(ii==(IIPAR+1))then
          distance(i) = Distance_Circle(curr_lon, curr_lat, X_mid(1), Y_mid(jj))
        else
          distance(i) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
        endif

    enddo


   if(ii==0)then
       distance(3) = Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), 90.0e+0_fp)
   else if(ii==(IIPAR+1))then
       distance(3) = Distance_Circle(curr_lon, curr_lat, X_mid(1), 90.0e+0_fp)
   else
       distance(3) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), 90.0e+0_fp)
   endif


   if(distance(3)==0.0)then
       do k=1,2
         kk = k + init_lev - 1
         wind_lonlat(k) = SUM(wind(:,init_lat,kk))/IIPAR
       enddo
   else
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
    endif

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

  real function Interplt_uv_PS_polar(i_uv, u_RLL, v_RLL, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)

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



    if(init_lat==0)then ! south polar
       jj = init_lat + 1
    endif

    if(init_lat==JJPAR)then ! north polar
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

       x_PS(i) = Re* cos(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
       y_PS(i) = Re* sin(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)

       do k=1,2
          kk = k + init_lev - 1
          if(i_uv==1)then ! i_ux==1 for u
             uv_PS(i,k) = -1.0* ( u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2) )
          endif

          if(i_uv==0)then ! for v
             uv_PS(i,k) = u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) - v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2)
          endif
       enddo
    enddo

    ! Third grid point: south/north polar point
    x_PS(3) = 0.0
    y_PS(3) = 0.0
    
    do k=1,2
       kk = k + init_lev - 1
          if(i_uv==1)then ! i_ux==1 for u
             ! interpolate all the grid points surrounding the polar point:
             do ii = 1,IIPAR
                uv_polars(ii) = -1.0* ( u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2) )
             enddo

             uv_PS(3,k) = SUM(uv_polars)/IIPAR
          endif

          if(i_uv==0)then ! for v
             do ii = 1,IIPAR
             uv_polars(ii) = u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) - v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2)
             enddo
             uv_PS(3,k) = SUM(uv_polars)/IIPAR
          endif
    enddo



    ! calculate the distance between particle and grid point
    do i = 1,3
          distance_PS(i) = sqrt( (x_PS(i)-curr_x)**2.0 + (y_PS(i)-curr_y)**2.0 )
    enddo

    ! Calculate the inverse distance weight
    do i = 1,3
        Weight_PS(i) = 1.0/distance_PS(i) / sum( 1.0/distance_PS(:) )
    enddo


    do k = 1,2
      uv_xy(k) =  Weight_PS(1) * uv_PS(1,k)   + Weight_PS(2) * uv_PS(2,k)   &
                 + Weight_PS(3) * uv_PS(3,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ (uv_xy(2)-uv_xy(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS_polar = uv_xy_lev

    return
  end function


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real function Interplt_uv_PS(i_uv, u_RLL, v_RLL, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)

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

    ! For pressure level, P_mid(1) is about surface pressure, has biggerst value.
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
             uv_PS(i,j,k) = -1.0* ( u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2) )
           endif

           if(i_uv==0)then ! for v
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) - v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      else

        x_PS(i,j) = -1.0* Re* cos(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)
        y_PS(i,j) = -1.0* Re* sin(X_mid(ii)*PI/180.0) / tan(Y_mid(jj)*PI/180.0)

        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) + v_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2)
           endif

           if(i_uv==0)then
             uv_PS(i,j,k) = -1* u_RLL(ii,jj,kk)*cos(X_mid(ii)*PI/180.0) / sin(Y_mid(jj)*PI/180.0) + v_RLL(ii,jj,kk)*sin(X_mid(ii)*PI/180.0) / (sin(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      endif

    enddo
    enddo

    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
          distance_PS(i,j) = sqrt( (x_PS(i,j)-curr_x)**2.0 + (y_PS(i,j)-curr_y)**2.0 )
    enddo
    enddo

    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2
        Weight_PS(i,j) = 1.0/distance_PS(i,j) / sum( 1.0/distance_PS(:,:) )
    enddo
    enddo


    do k = 1,2
      uv_xy(k) =  Weight_PS(1,1) * uv_PS(1,1,k)   + Weight_PS(1,2) * uv_PS(1,2,k)   &
                 + Weight_PS(2,1) * uv_PS(2,1,k) + Weight_PS(2,2) * uv_PS(2,2,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ (uv_xy(2)-uv_xy(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS = uv_xy_lev

    return
  end function


  ! calculation the great-circle distance between two points on the earth surface
  real function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2  ! unit is degree
    !real(fp) :: PI, Re
    
    ! Original equation, better for higher precision (double)
    Distance_Circle = Re * 2.0 * asin( (sin( (y1-y2)*PI/180.0 ))**2.0   & 
                      + cos(y1*PI/180.0) * cos(y2*PI/180.0) * (sin( 0.5*(x1-x2)*PI/180.0 ))**2.0 )

    ! This equation is better for lower precision (float4)
    !Distance_Circle = Re * ACOS( sin(y1*PI/180.0) * sin(y2*PI/180.0)   & 
    !                  + cos(y1*PI/180.0) * cos(y2*PI/180.0) * cos( (x1-x2)*PI/180.0 ) )

    return
  end function


  ! the linear interpolation function for vertical direction
  ! real function Line_Interplt(wind1, wind2, L1, L2, Lx)
  !  implicit none
  !  real(fp) :: wind1, wind2,  L1, L2
  !  real     :: Lx
  !
  !    Line_Interplt = wind1 + (wind2-wind1) / (L2-L1) * (Lx-L1)
  !       
  !    return
  !  end function

  !-------------------------------------------------------------------
  !=================================================================
  ! LAGRANGE_WRITE_STD begins here!
  !=================================================================

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
    INTEGER                   :: i_box
!    INTEGER                   :: L
    LOGICAL                   :: IsOldFile

    CHARACTER(LEN=255)            :: FILENAME

    FILENAME   = 'Lagrange_1day_box_i_lon_lat_lev.txt'
    tt = tt +1

    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)

!       OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='OLD', &
!             FORM='UNFORMATTED', ACCESS='Direct', Recl=18 )
    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='OLD', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       Do i_box = 1, n_boxes_max
          WRITE(261,'(2(x,I8),3(x,E12.5))') N_Dt_previous, i_box, box_lon(i_box), box_lat(i_box), box_lev(i_box)
!          WRITE(261,'(1(x,I8),3(x,E12.5))') i_box, box_lon(i_box), box_lat(i_box), box_lev(i_box)
       End Do
    
    ENDIF


  END SUBROUTINE lagrange_write_std


!---------------------------------------------------------------------

  subroutine lagrange_cleanup()

    if (allocated(box_lon)) deallocate(box_lon)
    if (allocated(box_lat)) deallocate(box_lat)
    if (allocated(box_lev)) deallocate(box_lev)
    if (allocated(box_length)) deallocate(box_length)
    if (allocated(box_depth)) deallocate(box_depth)
    if (allocated(box_width)) deallocate(box_width)

    WRITE(6,'(a)') '--> Lagrange Module Cleanup <--'

  end subroutine lagrange_cleanup


END MODULE Lagrange_Mod
