!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!-----------------------------------------------------------------------------!

MODULE Lagrange_Mod

  USE precision_mod
  USE PhysConstants, ONLY : PI, Re 
  USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: lagrange_init
  PUBLIC :: lagrange_run
  PUBLIC :: lagrange_write_std
  PUBLIC :: lagrange_cleanup


  integer, parameter :: n_boxes_max = 3*6000  ! 24*200*1 : lat*lon*lev
  integer            :: tt
  real(fp), allocatable :: box_lon(:)    
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)
  real(fp), allocatable :: box_depth(:)
  real(fp), allocatable :: box_width(:)
  real(fp), allocatable :: box_length(:)


CONTAINS


!-----------------------------------------------------------------

  SUBROUTINE lagrange_init( am_I_root )

    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    INTEGER                       :: i_box
    INTEGER                       :: ii, jj, kk
    CHARACTER(LEN=255)            :: FILENAME

    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Lagrnage Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'

    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))
    allocate(box_width(n_boxes_max))
    allocate(box_depth(n_boxes_max))
    allocate(box_length(n_boxes_max))


    do ii=1,6000,1  ! change the value of 'n_boxes_max'
        i_box = ii
        box_lon(i_box)  = -145.01e+0_fp                      ! -141.0W degree
        !box_lat(i_box) = -30.05e+0_fp + 0.1e+0_fp * ii     ! -29.95S : 29.95N : 0.1
        box_lat(i_box)  = -30.005e+0_fp + 0.01e+0_fp * ii   ! -29.995S : 299.95N : 0.1
        box_lev(i_box)  = 48.0e+0_fp                        ! 20.0hPa!
    enddo

    do ii=6001,12000,1  ! change the value of 'n_boxes_max'
        i_box = ii
        box_lon(i_box)  = -145.02e+0_fp                           ! -141.0W degree
        !box_lat(i_box) = -30.05e+0_fp + 0.1e+0_fp * ii           ! -29.95S : 29.95N :0.1
        box_lat(i_box)  = -30.005e+0_fp + 0.01e+0_fp * (ii-6000)  ! -29.995S : 299.95N: 0.1
        box_lev(i_box)  = 52.0e+0_fp                              ! 20.0hPa!
    enddo

    do ii=12001,n_boxes_max,1  ! change the value of 'n_boxes_max'
        i_box = ii
        box_lon(i_box)  = -145.03e+0_fp                            ! -141.0W degree
        !box_lat(i_box) = -30.05e+0_fp + 0.1e+0_fp * ii            ! -29.95S : 29.95N : 0.1
        box_lat(i_box)  = -30.005e+0_fp + 0.01e+0_fp * (ii-12000)  ! -29.995S : 299.95N : 0.1
        box_lev(i_box)  = 56.0e+0_fp                               ! 20.0hPa!
    enddo

!    box_lon    = 0.0e+0_fp       ! (/0.0e+0_fp, 5.0e+0_fp, 10.0e+0_fp/)
!    box_lat    = 0.0e+0_fp       ! (/0.0e+0_fp, 4.0e+0_fp, 8.0e+0_fp/)
!    box_lev    = 20.0e+0_fp      ! hPa
    box_width  = 0.0e+0_fp
    box_depth  = 0.0e+0_fp
    box_length = 0.0e+0_fp


!    FILENAME   = 'Lagrange_1hr_box_i_lon_lat_lev.txt'
    FILENAME   = 'Lagrange_1day_box_i_lon_lat_lev.txt'
    tt = 0

    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!    WRITE(261,'(a)') ' --> Lagrange Module Location (i_box,box_lon,box_lat,box_lev) <-- '

    Do i_box = 1, n_boxes_max
       WRITE(261,'(I0.4,3(x,E16.5E4))') i_box,box_lon(i_box),box_lat(i_box),box_lev(i_box)
    End Do


  END SUBROUTINE lagrange_init


!-----------------------------------------------------------------
!=================================================================

  SUBROUTINE lagrange_run(am_I_Root, State_Met, Input_Opt)

    USE Input_Opt_Mod, ONLY : OptInput
    ! USE PhysConstants, ONLY : PI, Re 
    ! Re    : Radius of Earth [m]

    !USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE TIME_MOD,      ONLY : GET_TS_DYN

    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE, XMID, YMID                 
    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON 
    ! USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON 
    ! DLAT( IIPAR, JJPAR, LLPAR ), DLON( IIPAR, JJPAR, LLPAR )
    ! XEDGE  ( IM+1, JM,   L ), YEDGE  ( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

    logical, intent(in) :: am_I_Root
    TYPE(MetState), intent(in) :: State_Met
    !TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in) :: Input_Opt

    REAL :: Dt          ! = 600.0e+0_fp          

    integer :: i_box

    integer :: i_lon
    integer :: i_lat
    integer :: i_lev

    real(fp) :: curr_lon
    real(fp) :: curr_lat
    real(fp) :: curr_pressure
    real(fp) :: X_edge2
    real(fp) :: Y_edge2

    real(fp) :: dbox_lon
    real(fp) :: dbox_lat
    real(fp) :: dbox_lev

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    real(fp) :: curr_u
    real(fp) :: curr_v
    real(fp) :: curr_omeg

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
    X_mid => XMID(:,1,1)   ! IIPAR+1
    Y_mid => YMID(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)

    ! Run Lagrangian advection HERE
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


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2) 
       i_lev = Find_iPLev(curr_pressure,P_edge)


       curr_u    = Interplt_wind(u,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_v    = Interplt_wind(v,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_omeg = Interplt_wind(omeg,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)


       dbox_lon = (Dt*curr_u) / (2.0*PI*Re*cos(curr_lat*PI/180.0)) * 360.0
       dbox_lat = (Dt*curr_v) / (PI*Re) * 180.0
       dbox_lev = Dt * curr_omeg / 100.0     ! Pa => hPa

       box_lon(i_box) = box_lon(i_box) + dbox_lon
       box_lat(i_box) = box_lat(i_box) + dbox_lat
       box_lev(i_box) = box_lev(i_box) + dbox_lev

       !call
       !grow_box(box_length(i_box),box_width(i_box),box_depth(i_box),i_lon,i_lat,i_lev)
    end do

    ! WRITE(6,*)'-Lagrange(i_lev)->',i_lev,P_edge(i_lev),P_edge(i_lev+1)

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
! functions to find location (i,j,k) of boxes 

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
! functions to interpolate wind speed (u,v,omeg) 
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
    if(curr_pressure==P_mid(i_lev))then
    if(curr_lon==X_mid(i_lon))then
    if(curr_lat==Y_mid(i_lat))then

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
      wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk)   + Weight(1,2) * wind(init_lon,init_lat+1,kk)   &
                      + Weight(2,1) * wind(init_lon+1,init_lat,kk) + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
    enddo


    ! second interpolate vertically (Linear)

    wind_lonlat_lev = wind_lonlat(1)    &
                     + (wind_lonlat(2)-wind_lonlat(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev), P_mid(i_lev+1), curr_pressure )

    Interplt_wind = wind_lonlat_lev

    return
  end function


  ! calculation the great-circle distance between two points on the earth surface
  real function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2  ! unit is degree
    !real(fp) :: PI, Re
    
    ! Original equation, better for higher precision (double)
    Distance_Circle = Re * 2.0 * ASIN( (sin( (y1-y2)*PI/180.0 ))**2.0   & 
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

!    FILENAME   = 'Lagrange_1hr_box_i_lon_lat_lev.txt'
    FILENAME   = 'Lagrange_1day_box_i_lon_lat_lev.txt'
    tt = tt +1

!    IF(mod(tt,6)==0)THEN     ! output once every hour
    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)

       OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='OLD', &
             FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       Do i_box = 1, n_boxes_max
          WRITE(261,'(I0.4,3(x,E16.5E4))') i_box, box_lon(i_box), box_lat(i_box),box_lev(i_box)
       End Do
    
    ENDIF

    !=================================================================
    ! LAGRANGE_WRITE_STD begins here!
    !=================================================================

!    RC        =  HCO_SUCCESS
!    Arr1D     => NULL()
!    Arr3D     => NULL()
!    Arr4D     => NULL()
!    timeVec   => NULL()
!
!
!!    nLon     = n_boxes_max
!    nLat     = n_boxes_max
!    nLev     = n_boxes_max
!    nTime    = 1
!
!    NcFile = 'lagrange_test.nc'

    !--------------------------------------------------------
    ! Check if file already exists. 
    ! If so, add new diagnostics to this file
    ! (instead of creating a new one)
    !--------------------------------------------------------

!    INQUIRE( FILE=NcFile, EXIST=IsOldFile ) 

!    IF ( IsOldFile ) THEN
!       CALL Ncop_Wr( fID, NcFile )
!       CALL NC_READ_TIME( fID, nTime, timeunit, timeVec, RC=RC )
!
!       ! new file will have one more time dimension
!       nTime = nTime + 1
!
!    ELSE

       ! Create output file
       ! Pass CREATE_NC4 to make file format netCDF-4 (mps, 3/3/16)
       ! Now create netCDF file with time dimension as UNLIMITED (bmy, 3/8/17)

!       CALL NC_Create( NcFile       = NcFile,                            &
!                       Title        = "lagrange_stdout_test",                   &
!                       nLon         = nLon,                              &
!                       nLat         = nLat,                              &
!                       nLev         = nLev,                           &
!                       nTime        = NF_UNLIMITED,                      &
!                       fId          = fId,                               &
!                       lonId        = lonId,                             &
!                       latId        = latId,                             &
!                       levId        = levId,                             &
!                       timeId       = timeId,                            &
!                       VarCt        = VarCt,                             &
!                       CREATE_NC4   =.TRUE.                             )
!
!    ENDIF

    !-----------------------------------------------------------------
    ! Write variables
    !-----------------------------------------------------------------

       ! Write longitude location of box ("box_lon(:)") to file   
!       CALL NC_Var_Def( fId         = fId,                                &
!                        lonId       = lonId,                              &
!                        latId       = -1,                                 &
!                        levId       = -1,                                 &
!                        timeId      = -1,                                 &
!                        VarName     = 'box_lon',                          &
!                        VarLongName = 'Longitude location of box',        &
!                        VarUnit     = 'degrees_east',                     &
!                        Axis        = 'X',                                &
!                        DataType    = f8,                                 &
!                        VarCt       = VarCt,                              &
!                        Compress    = .TRUE.                             )

!       ALLOCATE( Arr2D( nLon, nTime ) )
!       ALLOCATE( Arr2DOld( nLon, nTime ) )
!       Arr2D(:,nTime) = box_lon(:)
!       Arr2D(:,1:nTime-1) = Arr2DOld(:,:)
!       CALL NC_Var_Write( fId, 'box_lon', Arr2D=Arr1D )
!       DEALLOCATE( Arr2D )
!
!       CALL NC_READ_ARR( fID, TRIM(myName), 1, nlon, 1, nlat, &
!                        1, nlev, 1, ntime-1, ncArr=Arr4DOld, RC=RC )
!       Arr4D(:,:,:,1:ntime-1) = Arr4DOld(:,:,:,:)
     

       ! Write latitude location of box ("box_lat()") to file                                               
!       CALL NC_Var_Def( fId         = fId,                                &
!                        lonId       = -1,                                 &
!                        latId       = latId,                              &
!                        levId       = -1,                                 &
!!                        timeId      = -1,                                 &
!                        VarName     = 'box_lat',                          &
!                        VarLongName = 'Latitude lcoation of box',         &
!                        VarUnit     = 'degrees_north',                    &
!                        Axis        = 'Y',                                &
!                        DataType    = f8,                                 &
!                        VarCt       = VarCt,                              &
!                        Compress    = .TRUE.                             )
!
!       ALLOCATE( Arr1D( nLat ) )
!       Arr1D = box_lat
!       CALL NC_Var_Write( fId, 'box_lat', Arr1D=Arr1D )
!       DEALLOCATE( Arr1D )


    ! Close file
!    CALL NC_CLOSE ( fId )

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
