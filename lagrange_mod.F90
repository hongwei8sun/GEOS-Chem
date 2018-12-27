!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!-----------------------------------------------------------------------------!

MODULE Lagrange_Mod

  USE precision_mod
  USE CMN_SIZE_mod
  USE PhysConstants   ! Physical constants: Re, PI, PI_180

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: init_lagrange
  PUBLIC :: run_lagrange
  PUBLIC :: cleanup_lagrange


  integer, parameter :: n_boxes_max = 5000
  real(fp), allocatable :: box_lon(:)    !first use one point to test
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)
  real(fp), allocatable :: box_depth(:)
  real(fp), allocatable :: box_width(:)
  real(fp), allocatable :: box_length(:)
!
CONTAINS


!\\
! !INTERFACE:
!
  SUBROUTINE init_lagrange( am_I_root )

    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
!    REAL(f4), POINTER :: Ptr2D(:,:)

    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))
    allocate(box_width(n_boxes_max))
    allocate(box_depth(n_boxes_max))
    allocate(box_length(n_boxes_max))

    box_lon    = 0.0e+0_fp
    box_lat    = 0.0e+0_fp
    box_lev    = 0.0e+0_fp
    box_width  = 0.0e+0_fp
    box_depth  = 0.0e+0_fp
    box_length = 0.0e+0_fp

  END SUBROUTINE init_lagrange


!#
  SUBROUTINE run_lagrange(am_I_Root, State_Met, State_Chm, Input_Opt)

    USE Input_Opt_Mod, ONLY : OptInput
    USE PhysConstants, ONLY : PI, Re 

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE, XMID, YMID                  ! choose XEDGE or XMID ??        
    USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON
    ! DLAT( IIPAR, JJPAR, LLPAR ), DLON( IIPAR, JJPAR, LLPAR )
    ! XEDGE  ( IM+1, JM,   L ), YEDGE  ( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

    logical, intent(in) :: am_I_Root
    TYPE(MetState), intent(in) :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in) :: Input_Opt

    REAL :: Dt           = 600.0e+0_fp                   ! TS_CONV is the time step we need here ? 

    integer :: i_box

    integer :: i_lon
    integer :: i_lat
    integer :: i_lev

    real :: curr_lon
    real :: curr_lat
    real :: curr_pressure

    real :: dbox_lon
    real :: dbox_lat
    real :: dbox_lev

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)
    real(fp), pointer :: p_lev(:,:,:)

      ! LONEDG(I) = GET_XEDGE( IIPAR+1, 1, 1 )


    ! Establish pointers
    u => State_Met%U   ! figure out state_met%U is based on lat/lon or modelgrid(i,j)
    v => State_Met%V
    omeg => State_Met%OMEGA
    p_lev => State_Met%PMID  !Pressure (w/r/t moist air) at level centers (hPa)

    ! Run Lagrangian advection HERE
    do i_box = 1,n_boxes_max
       curr_lon = box_lon(i_box)
       curr_lat = box_lat(i_box)
       curr_pressure = box_lev(i_box)         ! 20.0e+0_fp ! hPa
       !curr_pressure = box_pressure(i_box)

       if (curr_lon > 180.0) cycle   ! Input lon should be between (-180,180)
       if (curr_lon < -180.0) cycle
       if (curr_lat > 89.0) cycle    ! Input lat should be between (-89,89)
       if (curr_lat < -89.0) cycle

       i_lon = find_longitude(curr_lon) !  WRITE THIS FUNCTION
       i_lat = find_latitude(curr_lat) !  
       i_lev = find_plev(curr_pressure,i_lon,i_lat,p_lev)

       dbox_lon = (Dt*u(i_lon,i_lat,i_lev)) /(2.0*PI*Re*cos(curr_lat*PI/180.0))*360.0
       dbox_lat = (Dt*v(i_lon,i_lat,i_lev)) / (PI*Re)*180.0
       dbox_lev = Dt*omeg(i_lon,i_lat,i_lev)

       box_lon(i_box) = box_lon(i_box) + dbox_lon
       box_lat(i_box) = box_lat(i_box) + dbox_lat
       box_lev(i_lev) = box_lev(i_box) + dbox_lev

       !call
       !grow_box(box_length(i_box),box_width(i_box),box_depth(i_box),i_lon,i_lat,i_lev)
    end do

    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)
    nullify(omeg)

  END SUBROUTINE run_lagrange


  integer function find_longitude(curr_lon, DLON, IIPAR, XMID)
    implicit none
    real :: curr_lon, DLON(:,:,:)  !注意DLON,XMID的维度
    find_longitude = INT( ANINT( (curr_lon-(XMID(1)))/DLON ) )+1
    ! for 5*4 horizontal resolution run (-180:5:175); ANINT 四舍五入函数
    return
  end function


  integer function find_latitude(curr_lat, DLAT, JJPAR, YMID)
    implicit none
    real :: curr_lat, DLAT
    find_latitude = INT( ANINT( (curr_lat-(YMID(1)))/DLAT ) )+1
    ! since we won't spread aerosols in polar region, ignore 89 and -89 in polar
    ! region. (-89,-86:4:86,89))
    return
  end function


  integer function find_plev(curr_pressure,i_lon,i_lat,p_lev)
    implicit none
    real :: curr_pressure
    real(fp), pointer :: p_lev(:,:,:)
    integer :: i_lon, i_lat
    integer :: locate(1)
    locate = MINLOC(abs( p_lev(i_lon,i_lat,:)-curr_pressure ))
    ! Notice that minimum difference in pressure (hPa) might not the minimum
    ! difference in distance (km)
    find_plev = INT(locate(1))
    return
  end function


!#

  subroutine cleanup_lagrange()
    if (allocated(box_lon)) deallocate(box_lon)
    if (allocated(box_lat)) deallocate(box_lat)
    if (allocated(box_lev)) deallocate(box_lev)
    if (allocated(box_length)) deallocate(box_length)
    if (allocated(box_depth)) deallocate(box_depth)
    if (allocated(box_width)) deallocate(box_width)
  end subroutine cleanup_lagrange

!EOC
END MODULE Lagrange_Mod
