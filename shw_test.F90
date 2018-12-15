module lagrange_mod

  use precision_mod
  use CMN_SIZE_mod
  use PhysConstants   ! Physical constants: Re, PI, PI_180

  implicit none
  private

  public init_lagrange
  public run_lagrange
  public cleanup_lagrange

  integer, parameter :: n_boxes_max = 5000
  real(fp), allocatable :: box_lon(:)    !first use one point to test
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_depth(:)
  real(fp), allocatable :: box_width(:)
  real(fp), allocatable :: box_length(:)

  subroutine init_lagrange(am_I_Root)
    logical, intent(in) :: am_I_Root

    ! Allocate statements
    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_width(n_boxes_max))
    allocate(box_depth(n_boxes_max))
    allocate(box_length(n_boxes_max))

    box_lon    = 0.0e+0_fp
    box_lat    = 0.0e+0_fp
    box_width  = 0.0e+0_fp
    box_depth  = 0.0e+0_fp
    box_length = 0.0e+0_fp

  end subroutine init_lagrange

  subroutine run_lagrange(am_I_Root, state_met, state_chm, input_opt)
    use state_chm_mod, only : chmstate
    use state_met_mod, only : metstate
    use input_opt_mod, only ：OptInput
    
    logical, intent(in) :: am_I_Root
    TYPE(MetState), intent(in) :: state_met
    TYPE(ChmState), intent(inout) :: state_chm
    TYPE(OptInput), intent(in) :: input_opt
    
    REAL(fp) :: REarth       = 6.375e+6_fp
    REAL(fp) :: Pi       = 3.14159265358979323e+0_fp     ! PI    : Double-Precision value of PI          
    REAL(fp) :: Pi_180   = PI / 180e+0_fp                ! PI_180 : Number of radians per degree

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    ! Establish pointers
    u => State_Met%U   ! figure out state_met%U is based on lat/lon or model grid(i,j)
    v => State_Met%V
    omeg => State_Met%OMEGA
    ! Run Lagrangian advection HERE

    do i_box = 1,n_boxes_max
       curr_lon = box_lon(i_box)
       curr_lat = box_lat(i_box)
       !curr_pressure = box_pressure(i_box)

       if (curr_lon > 180.0) cycle   ! Input lon should be between (-180,180)
       if (curr_lon < -180.0) cycle
       if (curr_lat > 89.0) cycle    ! Input lat should be between (-89,89)
       if (curr_lat < -89.0) cycle

       i_lon = find_longitude(curr_lon) !  WRITE THIS FUNCTION
       i_lat = find_latitude(curr_lat) !  
       ! Cheat for now
       curr_pressure = 20.0e+0_fp ! hPa
       i_lev = find_level(curr_pressure,i_lon,i_lat)
       dbox_lon = (time_step*u(i_lon,i_lat,i_lev)) / (2.0*Pi*REarth*cos(curr_lat*Pi/180.0))*360.0
       dbox_lat = (time_step*v(i_lon,i_lat,i_lev)) / (Pi*REarth)*180.0
       box_lon(i_box) = box_lon(i_box) + dbox_lon
       box_lat(i_box) = box_lat(i_box) + dbox_lat

       !call grow_box(box_length(i_box),box_width(i_box),box_depth(i_box),i_lon,i_lat,i_lev)
    end do

    !Dbox_lat = (time_step*v(:,:,:))/(Pi*REarth)*180.0            !model time_step length
    !Dbox_lon = (time_step*u(:,:,:))/(2.0*Pi*REarth*cos(lat*Pi/180.0))*360
!specify the exact location of u(:,:,:)

    !box_lat = box_lat + Dbox_lat
    !box_lon = box_lon + Dbox_lon

    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)

  end subroutine run_lagrange

  real function find_longitude(curr_lon)
    implicit none
    real :: curr_lon
    find_longitude = ( ANINT(curr_lon/5.0)*5.0 - (-180.0) )/5.0 +1   ! for 5*4 horizontal resolution run (-180:5:175); ANINT 四舍五入函数
    return
  end function

  real function find_latitude(curr_lat)
    implicit none
    real :: curr_lat
    find_latitude = ( ANINT((curr_lat+2.0)/4.0)*4.0-2.0 - (-90) )/4.0 +1  ! since we won't spread aerosols in polar region, ignore 89 and -89 in polar region. (-89,-86:4:86,89))
    return
  end function

  subroutine cleanup_lagrange()
    if (allocated(box_lon)) deallocate(box_lon)
    if (allocated(box_lat)) deallocate(box_lat)
    if (allocated(box_length)) deallocate(box_length)
    if (allocated(box_depth)) deallocate(box_depth)
    if (allocated(box_width)) deallocate(box_width)
  end subroutine cleanup_lagrange



end module lagrange_mod
