!-------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!--------------------------------------------------------------------------

! New version
! first simulate plume cross-section in 2D model
! Once the plume cross-section is highly distorted,
! simulate plume cross-section in 1D slab model.

! 1. dissolve plume once hitting tropopause

! 2. add variable 'N_split' to define the number of small parts after 
!    plume split

! 3. Oct 9th, 2020 
!  Update the horizontal stretch code, originally the horizontal plume segment is
!  calculated based on adjacent segments, but this requires all the plume
!  segments are exit in the whole simulation without dissolving. The new code
!  will calculate the horizontal stretch based on the center and one endpoint 
!  (box_lat_edge, box_lon_edge) of the segment, which will not rely on the
!   adjacent segment any more.

! 4. Oct 9th, 2020
! Assume that the box_length would not change because of PV=nRT, 
! only the plume cross-section will be used to represent the volume change because
! of PV=nRT

! 5. Oct 10th, 2020
! Add Third judge in 2D cross-section model.
! At the end of simulation, not only the 1D slab model will dissolve, 2D model
! will also dissolve.

! 6. Oct 10th, 2020
! Inn stead of using adjacent plume segment, the box_alpha is calculated based
! on the center and one endpoint inside one plume segment

! 7. Oct 12nd, 2020
! delete the output for grid area (State_Met_AREA_M2.txt), since this variable
! can be found in the formal output as "AREA".

! 8. Oct 12nd, 2020
! correct the runge-kutta method for advection process (dx/dt=u).
! since the wind speed is 3-h average value, assume wind speed is independent of
! time, only change with locations [u(z,y,x) instead of u(t,z,y,x)]
! RK_Dlev(Ki) = RK_Dt(1) * curr_omeg / 100.0 

! 9. Oct 12nd, 2020
! use runge-kutta method when calcualte box_u(i_box), box_v(i_box),
! box_omeg(i_box)

! 10. Oct 12nd, 2020
! modify the box_alpha and ensure the wind speed along plume cross-section is
! correct. range of box_alpha is [PI, -PI)

! 11. Oct 12nd, 2020
! modify the 1D slab model splitting, correct the box_lon and  box_lon_edge 
! for the new plume after splitting
! to get the correct new box_lon, box_length should times COS(box_alpha) first.

! 12. Oct 13rd, 2020
! correct the calculate for plume splitting box_lon and box_lat, as well as
! box_lon_edge and box_lat_edge.

! Oct 15th, 2020
! add forth judge to 2D cross-section model, let plume dissolve when transport
! from stratosphere into troposphere.

! Nov 1st, 2020
! Instead of setting Pdt=10s, adjust Pdt based on the CFL condition in 2D grids.
!
! Dt2 for 1D grid can only be 60s as minimum. If a less than 60s values is
! needed for Dt2, let 1D slab begin to combine 2 grid cells into 1 grid cell.

! Nov 2, 2020
! Correct the injected aerosol concentration at the beginning
! pay attention to the unit

! Nov 2, 2020
! Add Xscale and Yscale for 2D grid, which can be used to identify the box_theta
! value, when concentratino don't follow gaussian distribution after involve
! settling speed etc.
! 
! when change 2D to 1D, not use degree to judge, but use length (Xscale, Yscale) 
! containing 95% mass to judge

! Nov 3, 2020
! add a lifetime variable Lifetime(n_box) to record the plume lifetim in plume
! model. When lifetime is larger than 3 hours, begin to consider change 2D to 1D
! change Sigma2D to lifetime

! Nov 3, 2020
! redifine LenB = Pdy(Ibox) in Slab_init() function,
! set n_slab_max = 300, which is close to n_y_max;

! Nov 7, 2020
! 1. DlonN, DlatN are corrected: the circumference of the lat circle is 2*pi*r,
! instead of 2*pi*r**2 !!!
!
! 2. check the shadow plume, to ensure the total mass is conseved 
! Mr + Mb = (Mp - Ms) + Mb = (Mp + D_Mp) - (Ms + D_Ms) + [Mb - ( D_Mp - D_Ms)]
! Mr: total mass of the real plume;
! Mp: total mass of the plume model in the code;
! Ms: total mass of the shadow plume model in the code;
! Mb: total mass of the background grid cell containing plume model
! D_Mp: change mass after one time step in plume model
! D_Ms: change mass after one time step in shadow plume model

! Nov 8, 2020
! Instead of using 2nd order chemical reaction to decide the initail lenght of
! 1D grid, use the lenght containing 95% of total mass.

! Nov 11, 2020
! add following code to check cpu time for each section of code:
!       call cpu_time(start)
!       put code to test here
!       call cpu_time(finish)
!       print '("Time = ",f6.3," seconds.")',finish-start

! Nov 11, 2020
! CFL_2d is added to speed up the code in advection in 2D grid model

! Nov 14, 2020
! the concentration will change after the volume change because of idea gas law
! (PV=nRT)

! Nov 14, 2020
! track the whole lifetime of the first plume segment to check:
! (1) shadow plume model performance;
! (2) evolution of the plume is reasonable;

! Nov 16, 2020
! delete shadow plume model, replace by a exrta mass container
!
! reset the injected aerosol concentration in GCM model in lagrangain_initial
! module, to make sure the background concentration is 0 at beginning.
!
! For changing from 2D to 1D in function Interplt_2D() in module Slab_init(), 
! when the interpolate point is out of 2D grid, instead of setting as 0, set 
! as background concentration.
!       IF(Ix0<=1 .or. Ix0>=n_x_max-1)THEN
!        Interplt_2D = 0.0e+0_fp

! Nov 17, 2020
! Add two processes:
! When combine 9 to 1 grid cell in 2D grid model, extra mass will increase.
! When combine 2 to 1 grid cell in 1D grid model, extra mass will increase.

! Nov 19, 2020
! assume plume segment lenght is always along the background wind direction, so
! the box_alpha is calculated based on the background wind speed now.
!
! use Lyapunov component (Ly) to calculate the horizontal stretch of the plume
! segment.
!
! delete the box_lon_edge(i_box), box_lat_edge(i_box), box_lev_edge(i_box)

! Nov 21, 2020
! modify the next_i_lon and next_i_lat to ensure they are in the range of 
! [1, IIPAR/JJPAR].
!
!        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
!        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR
!
!        if(next_i_lat>JJPAR) next_i_lat=JJPAR
!        if(next_i_lat<1) next_i_lat=1

! Nov 27, 2020
!
! update the 1D grid domain, make sure the number of 1D grid cell cam be divided
! by 4, to ensure n_slab_25 is an integer.
!
! add a new variable n_slab_max2=n_slab_max, Concnt1D_bdy(n_slab_max2), to add
! two more grid cell temporarily when calculate the diffusion process inside 1D
! grid.
!       Concnt1D_bdy(1)          = backgrd_concnt(1)
!       Concnt1D_bdy(n_slab_max2) = backgrd_concnt(1)
!       Concnt1D_bdy(2:n_slab_max2-1) = box_concnt_1D(i_box,1:n_slab_max,1)

! Nov 29, 2020
! find the error about total injected mass increasing in plume model:
! box_Ra and box_Rb are commented previously. So only box_length keeps
! increasing. Volume is not conserved.
!      length0           = box_length(i_box)
!      box_length(i_box) = EXP(Ly*Dt) * box_length(i_box)
!      IF(Judge_plume(i_box)==1) THEN
!        box_Ra(i_box) = box_Ra(i_box)*SQRT(length0/box_length(i_box))
!        box_Rb(i_box) = box_Rb(i_box)*SQRT(length0/box_length(i_box))
!
! ERROR: total mass is not conservsed between Plume model and Eulerian model
! Check the third judgement:
! ERROR solved by 
! changing  "IF ( ITS_TIME_FOR_EXIT() ) THEN" 
! to        "IF(DAY==31 .and. HOUR==23)THEN"

! Dec 5, 2020
! add bondary grid cell to 2D grid, similar like what we did for 1D grid.
! By doing so, the boundary of plume would not belong to plume mass, when 
! combine plume grid cells (9 to 1, 2 to 1), these extra mass in plume 
! boundary grid cell will be included.
!
! the dimension order of box_concnt_2D and box_concnt_1D are changed, "n_box"
! dimension are put in the last dimension.

! Dec 6, 2020
! add box_alpha(n_box), and delete box_u, box_v, box_w
!
! delete box_Ptemp(i_box)


! Dec 21, 2020
! ERROR: number of 2D plume model become super larger at 25th day
! Guess the 2D plume model could not change to 1D plume model for some reason,
! check the XY scale and time.
! add a new criteria, once the 2D plume model live for more than 6 hours, force
! the 2D plume model change to 1D plume model.
!
! add curr_Ptemp interplating in plume_module


! when plume is dissolved, delete it from the n_box


! reconsider plume dissolve criteria, which should involve extra mass


!
! dissolve the plume when touch the top of the atmosphere


! check: delete or release all the comment code




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
  PUBLIC :: n_box
  PUBLIC :: n_box_prev
!  PUBLIC :: box_u, box_v, box_omeg
  PUBLIC :: box_alpha
!  PUBLIC :: box_Ptemp

!  PUBLIC :: Plume_I, Plume_J, Plume_L ! (n_box)
  PUBLIC :: Judge_plume ! 0-dissolved; 1-1D slab model; 2-2D model

  PUBLIC :: n_x_max, n_y_max
  PUBLIC :: n_x_mid, n_y_mid
  PUBLIC :: box_concnt_2D ! (n_x_max, n_y_max, N_specie, n_box)

  PUBLIC :: n_slab_max, n_slab_25, n_slab_50, n_slab_75
  PUBLIC :: box_concnt_1D ! (n_slab_max, N_specie, n_box)

!  PUBLIC :: Extra_mass_2D       ! [molec cm-3]
!  PUBLIC :: Extra_mass_1D       ! [molec cm-3]

  PUBLIC :: Total_extra


  real(fp), allocatable :: box_lon(:), box_lat(:), box_lev(:)
!  real(fp), allocatable :: box_lon_edge(:), box_lat_edge(:), box_lev_edge(:)
!  integer, parameter    :: n_box = 60 !6904224     
  integer               :: n_box
   
  integer               :: n_box_prev
!  real(fp), allocatable :: box_u(:), box_v(:), box_omeg(:)
  real(fp), allocatable :: box_alpha(:)
!  real(fp), allocatable :: box_Ptemp(:)


!  integer, allocatable  :: Plume_I(:), Plume_J(:), Plume_L(:) !(n_box)

  integer, parameter    :: n_x_max = 483    ! number of x grids in 2D
  integer, parameter    :: n_y_max = 165    !279     ! number of y grids in 2D
  ! the odd number of n_x_max can ensure a center grid
  integer, parameter    :: n_x_mid = 242
  integer, parameter    :: n_y_mid = 83

  integer, parameter    :: n_x_max2 = n_x_max+2
  integer, parameter    :: n_y_max2 = n_y_max+2

  ! n_slab_max should be divided by 4, to ensure n_slab_25 is an integer.
  integer, parameter    :: n_slab_max = 164 ! close to n_y_max, number of slabs in 1D
  integer, parameter    :: n_slab_max2 = n_slab_max+2
  ! add 2 more slab grid to containing background concentration

  integer               :: n_slab_25, n_slab_50, n_slab_75

  real, allocatable     :: Pdx(:)          ! [m] 2D horizontal resolution
  real, allocatable     :: Pdy(:)          ! [m] 2D vertical resolution
  real, parameter       :: Dx_init = 100
  real, parameter       :: Dy_init = 10

  ! medical concentration of each slab [molec/cm3]
  integer, allocatable  :: Judge_plume(:)
  real(fp), allocatable :: box_concnt_1D(:,:,:)   !(N_slab,N_specie,n_box)
  real(fp), allocatable :: box_concnt_2D(:,:,:,:) 


  real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)
  real(fp), pointer :: P_edge(:)


  integer, parameter    :: N_parcel   = 12 ! 131        
  integer               :: tt     
  ! Aircraft would release 131 aerosol parcels every time step

!  integer               :: i_rec
  integer               :: N_specie


  ! box_radius(n_box) of the cross-section
  real(fp), allocatable :: box_Ra(:) ! long length [m]
  real(fp), allocatable :: box_Rb(:) ! short/vertical side [m]

  ! Theta is the clockwise angle between z-axis (P) and vertical Ra
  real(fp), allocatable :: box_theta(:)     ! degree between vertical and Ra
  real(fp), allocatable :: box_length(:)
  real(fp), allocatable :: Lifetime(:) ! record plume lifetime in plume model



  ! Used for plume_run subroutine:


!  real(fp), allocatable :: Extra_mass_2D(:,:,:,:)
  ! [n_box_max,n_x_max,n_y_max,N_specie]
!  real(fp), allocatable :: Extra_mass_1D(:,:,:) ![]

  real(fp), allocatable :: Total_extra(:)

  integer(fp), allocatable :: Data1D_Int(:)
  real(fp), allocatable :: Data1D(:), Data2D(:,:)
  real(fp), allocatable :: Data3D(:,:,:), Data4D(:,:,:,:)

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

    INTEGER                       :: i_box, i_slab
    INTEGER                       :: ii, jj, kk
    CHARACTER(LEN=255)            :: FILENAME, FILENAME_INIT
    CHARACTER(LEN=255)            :: FILENAME2

    integer :: i_lon, i_lat, i_lev            !1:IIPAR

    REAL(fp) :: lon1, lat1, lon2, lat2
    REAL(fp) :: box_lon_edge, box_lat_edge        


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

    N_specie = 1
    n_box  = N_parcel


    allocate(box_lon(n_box))
    allocate(box_lat(n_box))
    allocate(box_lev(n_box))

!    allocate(box_lon_edge(n_box))
!    allocate(box_lat_edge(n_box))
!    allocate(box_lev_edge(n_box))

!    allocate(box_u(n_box))
!    allocate(box_v(n_box))
!    allocate(box_omeg(n_box))
    allocate(box_alpha(n_box))
!    allocate(box_Ptemp(n_box))

!--------------------------------------------------
! Allocate following variables for plume_run
!--------------------------------------------------
!    allocate(Plume_I(n_box))
!    allocate(Plume_J(n_box))
!    allocate(Plume_L(n_box))

    allocate(Pdx(n_box))
    allocate(Pdy(n_box))

    allocate(box_Ra(n_box))
    allocate(box_Rb(n_box))

    allocate(box_theta(n_box))
    allocate(box_length(n_box))

!    N_specie = State_Chm%nAdvect

    allocate(Judge_plume(n_box))
    allocate(box_concnt_2D(n_x_max,n_y_max,N_specie,n_box))
    allocate(box_concnt_1D(n_slab_max,N_specie,n_box))


!    allocate(Extra_mass_2D(n_box, n_x_max, n_y_max, N_specie))
!    allocate(Extra_mass_1D(n_box, n_slab_max, N_specie))

    allocate(Total_extra(n_box))

    allocate(Lifetime(n_box))


    n_slab_25 = (n_slab_max-2)/4*1
    n_slab_50 = (n_slab_max-2)/4*2
    n_slab_75 = (n_slab_max-2)/4*3


    X_mid  => XMID(:,1,1)   ! IIPAR
    Y_mid  => YMID(1,:,1)
    P_mid  => State_Met%PMID(1,1,:)  ! Pressure at level centers (hPa)

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]


!--------------------------------------------------
! Set the initail value of box location, length
!--------------------------------------------------

    ! (1) run month by month
!    FILENAME_INIT   = '/n/home12/hongwei/hongwei/merra2_2x25_standard_Dec/Lagrange_Init_box_i_lon_lat_lev.txt'
!    OPEN( 361,      FILE=TRIM( FILENAME_INIT   ), STATUS='old', &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!    Do i_box = 1, n_box
!       READ(361,'( 2(x,I8) , 3(x,E12.5) )') n_box_prev, iibox, box_lon(i_box), box_lat(i_box), box_lev(i_box)
!    End Do
    
!--------------------------------------------------

    ! (2)
    do i_box = 1,n_box,1
        box_lon(i_box) = -141.0e+0_fp     ! 0   
        box_lat(i_box) = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) &
                        * (-1.0)**FLOOR(i_box/6000.0) ! -29.95S : 29.95N : 0.1
        box_lev(i_box) = 52.0e+0_fp       ! [hPa] at about 20 km

!        box_lev_edge(i_box) = box_lev(i_box)       ! [hPa] at about 20 km

        
        Pdx(i_box) = Dx_init              ! [m]
        Pdy(i_box) = Dy_init              ! [m]


        box_theta(i_box)   = 0.0e+0_fp    ! [radian]
        Lifetime(i_box)     = 0.0e+0_fp    ! [second]
        Judge_plume(i_box) = 2
    enddo

!--------------------------------------------------------
! Set the initial value of plume character:
! box_theta, box_length,
!  max/min radius for each slab
!--------------------------------------------------------

!    box_theta    = 0.0e+0_fp     ! [radian]

    DO i_box = 1, n_box-1
      box_lon_edge = 0.5 * ( box_lon(i_box) + box_lon(i_box+1) )
      box_lat_edge = 0.5 * ( box_lat(i_box) + box_lat(i_box+1) )
      box_length(i_box)   = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
                                   box_lon_edge, box_lat_edge)  ! [m]
    ENDDO

    i_box = n_box
    box_lon_edge = box_lon(i_box) + 0.5* ( box_lon(i_box) - box_lon(i_box-1) )
    box_lat_edge = box_lat(i_box) + 0.5* ( box_lat(i_box) - box_lat(i_box-1) )
    box_length(i_box) = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
                                   box_lon_edge, box_lat_edge)  ! [m]

    !--------------------------------------------------------
    ! Set the initial concentration of injected aerosols
    ! Here assume the injection rate is 30 kg/km 
    ! (same as 30 g/m) for H2SO4
    !--------------------------------------------------------
    ! State_Chm%nAdvect: the last one is PASV
    DO i_box = 1, n_box
      box_concnt_2D(:,:,N_specie,i_box) = 0.0e+0_fp
      box_concnt_2D(n_x_mid,n_y_mid,N_specie,i_box) = &
                                box_length(i_box)*30.0 &
                              / (Pdx(i_box)*Pdy(i_box)*box_length(i_box))

      ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
      box_concnt_2D(n_x_mid,n_y_mid,N_specie,i_box) = &
        box_concnt_2D(n_x_mid,n_y_mid,N_specie,i_box) &
          / 1.0e+6_fp / 98.0 * AVO

    ENDDO


    tt   = 0


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

!    FILENAME   = 'Lagrange_Init_xyz_' // TRIM(ADJUSTL(YEAR_C)) &
!         // '-' //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) &
!         // '-' // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
!         // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'

!    OPEN( 261,      FILE=TRIM( FILENAME ), STATUS='REPLACE', &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!    Do i_box = 1, n_box
!       WRITE(261,'( 2(x,I8) , 3(x,E12.5) )')n_box_prev, i_box, box_lon(i_box), &
!                                            box_lat(i_box), box_lev(i_box)
!    End Do


!-----------------------------------------------------------------
! Output box concentration
!-----------------------------------------------------------------
!    FILENAME2   = 'Plume_Init_concentration_molec_' // TRIM(ADJUSTL(YEAR_C)) // '-' &
!        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' // &
!        TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) // ':' // &
!        TRIM(ADJUSTL(SECOND_C)) // '.txt'

!    OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!    DO i_box = 1, 1001, 100
!      WRITE(262,*) i_box
!      WRITE(262,*) SUM( box_concnt_2D(i_box,:,:,N_specie) &
!      )*Pdx*Pdy*box_length(i_box)*1e+6_fp
!       ! [molec]
!      WRITE(262,*) box_concnt(i_box,:,:,N_specie)
!    ENDDO

!     WRITE(262,*) box_concnt(1,:,N_specie)

!-----------------------------------------------------------------
! Output State_Met%AD(i_lon,i_lat,i_lev) into State_Met_AD.txt
!-----------------------------------------------------------------
    OPEN( 314,      FILE='State_Met_AD.txt', STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    Do i_lon = 1, IIPAR
    Do i_lat = 1, JJPAR
    Do i_lev = 1, LLPAR
       WRITE(314,'(x,E12.5)') State_Met%AD(i_lon,i_lat,i_lev) ![kg]
    End Do
    End Do
    End Do

!-----------------------------------------------------------------------
! Output State_Met%AREA_M2(i_lon,i_lat,i_lev) into State_Met_AREA_M2.txt
!-----------------------------------------------------------------------
!    OPEN( 315,      FILE='State_Met_AREA_M2.txt', STATUS='REPLACE', &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!    Do i_lon = 1, IIPAR
!    Do i_lat = 1, JJPAR
!    Do i_lev = 1, LLPAR
!       WRITE(315,'(x,E12.5)') GET_AREA_M2(i_lon,i_lat,i_lev) ! [m2]
!    End Do
!    End Do
!    End Do





    ! set initial background concentration of injected aerosol as 0 in GCM
    State_Chm%Species(:,:,:,State_Chm%nAdvect) = 0.0e+0_fp  ! [kg/kg]
    State_Chm%Species(:,:,:,State_Chm%nAdvect-1) = 0.0e+0_fp  ! [kg/kg]



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

    REAL(fp) :: Dt, RK_Dt(5)                  ! = 600.0e+0_fp          

    real(fp), pointer  :: PASV_EU           
    integer            :: nAdv        
    REAL(fp)           :: MW_g

    integer :: i_box

    ! 1:IIPAR
    integer :: i_lon, i_lat, i_lev
    integer :: next_i_lon, next_i_lat, next_i_lev

    integer :: ii, jj, kk

    integer :: Ki

    integer :: n_active_plume


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
    real(fp) :: RK_u(4), RK_v(4), RK_omeg(4)

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)
    real(fp), pointer :: Ptemp(:,:,:)
    real(fp), pointer :: T1(:,:,:)
    real(fp), pointer :: T2(:,:,:)

    real(fp) :: curr_T1, next_T2, ratio
    real(fp) :: V_prev, V_new

    real(fp) :: box_u, box_v, box_omeg
    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: curr_u_PS, curr_v_PS

    real(fp), pointer :: X_edge(:), Y_edge(:)

    real(fp) :: Dx, Dy

    real(fp) :: length0
    real(fp) :: lon1, lon2, lat1, lat2


    real(fp) :: box_lon_edge, box_lat_edge

!    real(fp) :: box_alpha
    real(fp) :: D_wind, D_x, D_y
    real(fp) :: Ly ! Lyapunov exponent [s-1]


    Dt = GET_TS_DYN()

    RK_Dt(1) = 0.0
    RK_Dt(2) = 0.5*Dt
    RK_Dt(3) = 0.5*Dt
    RK_Dt(4) = Dt
    RK_Dt(5) = 0.0

    u     => State_Met%U   ! figure out state_met%U is based on
                           !  lat/lon or modelgrid(i,j)
    v     => State_Met%V   ! V [m s-1]
    omeg  => State_Met%OMEGA     ! Updraft velocity [Pa/s]
    Ptemp => State_Met%THETA     ! Potential temperature [K]
    T1    => State_Met%TMPU1     ! Temperature at start of
                                 !  timestep [K]
    T2    => State_Met%TMPU2     ! Temperature at end of
                                 !  timestep [K]

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!
    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)
     
    
!    if(N_curr<=n_box)then
!      n_box = N_curr
!    else
!      n_box = n_box
!    endif

! add new box every time step


    n_box_prev = n_box
    n_box = n_box + N_parcel

    allocate(Data1D_Int(n_box_prev))
    allocate(Data1D(n_box_prev))
    allocate(Data2D(N_specie, n_box_prev))
    allocate(Data3D(n_slab_max, N_specie, n_box_prev))
    allocate(Data4D(n_x_max, n_y_max, N_specie, n_box_prev))

    Data1D = box_lon(:)
    deallocate(box_lon)
    allocate(box_lon(n_box))
    box_lon(1:n_box_prev) = Data1D

    Data1D = box_lat(:)
    deallocate(box_lat)
    allocate(box_lat(n_box))
    box_lat(1:n_box_prev) = Data1D

    Data1D = box_lev(:)
    deallocate(box_lev)
    allocate(box_lev(n_box))
    box_lev(1:n_box_prev) = Data1D

!    Data1D = box_lon_edge(:)
!    deallocate(box_lon_edge)
!    allocate(box_lon_edge(n_box))
!    box_lon_edge(1:n_box_prev) = Data1D

!    Data1D = box_lat_edge(:)
!    deallocate(box_lat_edge)
!    allocate(box_lat_edge(n_box))
!    box_lat_edge(1:n_box_prev) = Data1D

!    Data1D = box_lev_edge(:)
!    deallocate(box_lev_edge)
!    allocate(box_lev_edge(n_box))
!    box_lev_edge(1:n_box_prev) = Data1D


    Data1D = Pdx(:)
    deallocate(Pdx)
    allocate(Pdx(n_box))
    Pdx(1:n_box_prev) = Data1D

    Data1D = Pdy(:)
    deallocate(Pdy)
    allocate(Pdy(n_box))
    Pdy(1:n_box_prev) = Data1D

    Data1D_Int = Judge_plume(:)
    deallocate(Judge_plume)
    allocate(Judge_plume(n_box))
    Judge_plume(1:n_box_prev) = Data1D_Int


    Data4D = box_concnt_2D(:,:,:,:)
    deallocate(box_concnt_2D)
    allocate(box_concnt_2D(n_x_max,n_y_max,N_specie,n_box))
    box_concnt_2D(:,:,:,1:n_box_prev) = Data4D

!    Data4D = Extra_mass_2D(:,:,:,:)
!    deallocate(Extra_mass_2D)
!    allocate(Extra_mass_2D(n_box,n_x_max,n_y_max,N_specie))
!    Extra_mass_2D(1:n_box_prev,:,:,:) = Data4D


    Data3D = box_concnt_1D(:,:,:)
    deallocate(box_concnt_1D)
    allocate(box_concnt_1D(n_slab_max,N_specie,n_box))
    box_concnt_1D(:,:,1:n_box_prev) = Data3D

!    Data3D = Extra_mass_1D(:,:,:)
!    deallocate(Extra_mass_1D)
!    allocate(Extra_mass_1D(n_box,n_slab_max,N_specie))
!    Extra_mass_1D(1:n_box_prev,:,:) = Data3D


    Data1D = box_Ra(:)
    deallocate(box_Ra)
    allocate(box_Ra(n_box))
    box_Ra(1:n_box_prev) = Data1D

    Data1D = box_Rb(:)
    deallocate(box_Rb)
    allocate(box_Rb(n_box))
    box_Rb(1:n_box_prev) = Data1D

    Data1D = box_theta(:)
    deallocate(box_theta)
    allocate(box_theta(n_box))
    box_theta(1:n_box_prev) = Data1D

    Data1D = box_length(:)
    deallocate(box_length)
    allocate(box_length(n_box))
    box_length(1:n_box_prev) = Data1D



    Data1D = Total_extra(:)
    deallocate(Total_extra)
    allocate(Total_extra(n_box))
    Total_extra(1:n_box_prev) = Data1D


!    Data1D = box_u(:)
!    deallocate(box_u)
!    allocate(box_u(n_box))
!    box_u(1:n_box_prev) = Data1D

!    Data1D = box_v(:)
!    deallocate(box_v)
!    allocate(box_v(n_box))
!    box_v(1:n_box_prev) = Data1D

!    Data1D = box_omeg(:)
!    deallocate(box_omeg)
!    allocate(box_omeg(n_box))
!    box_omeg(1:n_box_prev) = Data1D

    Data1D = box_alpha(:)
    deallocate(box_alpha)
    allocate(box_alpha(n_box))
    box_alpha(1:n_box_prev) = Data1D

!    Data1D = box_Ptemp(:)
!    deallocate(box_Ptemp)
!    allocate(box_Ptemp(n_box))
!    box_Ptemp(1:n_box_prev) = Data1D

!    Data1D_Int = Plume_I(:)
!    deallocate(Plume_I)
!    allocate(Plume_I(n_box))
!    Plume_I(1:n_box_prev) = Data1D_Int

!    Data1D_Int = Plume_J(:)
!    deallocate(Plume_J)
!    allocate(Plume_J(n_box))
!    Plume_J(1:n_box_prev) = Data1D_Int

!    Data1D_Int = Plume_L(:)
!    deallocate(Plume_L)
!    allocate(Plume_L(n_box))
!    Plume_L(1:n_box_prev) = Data1D_Int

    Data1D = Lifetime(:)
    deallocate(Lifetime)
    allocate(Lifetime(n_box))
    Lifetime(1:n_box_prev) = Data1D


    deallocate(Data1D_Int)
    deallocate(Data1D)
    deallocate(Data2D)
    deallocate(Data3D)
    deallocate(Data4D)

    !-----------------------------------------------------------------------
    ! set initial location for new added boxes
    !-----------------------------------------------------------------------

    do i_box = n_box_prev+1, n_box, 1
        box_lon(i_box)   = -141.0e+0_fp   ! 
        box_lat(i_box)   = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) &
                          * (-1.0)**FLOOR(i_box/6000.0) 
        ! -29.995S : 29.995N : 0.01
        box_lev(i_box)   = 52.0e+0_fp       ! about 20 km

!        box_lev_edge(i_box) = box_lev(i_box)

        Pdx(i_box)       = Dx_init
        Pdy(i_box)       = Dy_init

        box_theta(i_box) = 0.0e+0_fp      ! [radian]
        Lifetime(i_box)   = 0.0e+0_fp     ! [second]
        Judge_plume(i_box) = 2

    enddo

    !--------------------------------------------------------
    ! Set the initial value of plume character:
    ! box_theta, box_length, Lifetime
    !--------------------------------------------------------
    DO i_box = n_box_prev+1, n_box-1
      box_lon_edge = 0.5 * ( box_lon(i_box) + box_lon(i_box+1) )
      box_lat_edge = 0.5 * ( box_lat(i_box) + box_lat(i_box+1) )
      box_length(i_box)   = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
                                box_lon_edge, box_lat_edge)  ! [m]
    ENDDO


    i_box = n_box
    box_lon_edge = box_lon(i_box) + 0.5* ( box_lon(i_box) - box_lon(i_box-1) )
    box_lat_edge = box_lat(i_box) + 0.5* ( box_lat(i_box) - box_lat(i_box-1) )
    box_length(i_box) = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
                                box_lon_edge, box_lat_edge)  ! [m]


    !--------------------------------------------------------
    ! Set the initial concentration of injected aerosols
    ! Here assume the injection rate is 30 kg/km for H2SO4
    !--------------------------------------------------------
    ! State_Chm%nAdvect: the last one is PASV
    DO i_box = n_box_prev+1, n_box
      box_concnt_2D(:,:,N_specie,i_box) = 0.0e+0_fp
      box_concnt_2D(n_x_mid,n_y_mid,N_specie,i_box) = &
                                box_length(i_box)*30.0 &
                              / (Pdx(i_box)*Pdy(i_box)*box_length(i_box))

      ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
      box_concnt_2D(n_x_mid,n_y_mid,N_specie,i_box) = &
        box_concnt_2D(n_x_mid,n_y_mid,N_specie,i_box) &
          / 1.0e+6_fp / 98.0 * AVO

    ENDDO



    n_active_plume = 0

    !-----------------------------------------------------------------------
    ! Run Lagrangian trajectory-track HERE
    !-----------------------------------------------------------------------
    do i_box = 1,n_box,1

      IF(Judge_plume(i_box)==0) GOTO 500

      ! record the plume lifetime before dissolving into Eulerian model
      Lifetime(i_box) = Lifetime(i_box) + Dt


      !-----------------------------------------------------------------------
      ! For the center of the plume
      !-----------------------------------------------------------------------

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
!      curr_T1       = box_T1(i_box)        ! K


      ! check this grid box is the neast one or the one located in left-bottom
      ! ???
      i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(curr_pressure,P_edge)


      ! ====================================================================
      ! Add concentraion of PASV into conventional Eulerian GEOS-Chem in 
      ! corresponding with injected parcels in Lagrangian model
      ! For conventional GEOS-Chem for comparison with Lagrangian Model:
      !
      !  AD(I,J,L) = grid box dry air mass [kg]
      !  AIRMW     = dry air molecular wt [g/mol]
      !  MW_G(N)   = specie molecular wt [g/mol]
      !     
      ! the conversion is:
      ! 
      !====================================================================
      if(i_box>n_box_prev)then
         nAdv = State_Chm%nAdvect         ! the last one is PASV_EU
         PASV_EU => State_Chm%Species(i_lon,i_lat,i_lev,nAdv)  ! [kg/kg]

         MW_g = State_Chm%SpcData(nAdv)%Info%emMW_g
         ! Here assume the injection rate is 30 kg/km for H2SO4: 
         PASV_EU = PASV_EU + box_length(i_box)*1.0e-3_fp*30.0 &
                                        /State_Met%AD(i_lon,i_lat,i_lev)

         ! write(6,*)'== test 1 ==>', State_Chm%Spc_Units
      endif


      !--------------------------------------------------------------------
      ! interpolate potential temperature for Plume module:
      !--------------------------------------------------------------------
!      if(abs(curr_lat)>Y_mid(JJPAR))then
!         box_Ptemp(i_box) = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
!                                 i_lev, curr_lon, curr_lat, curr_pressure)
!      else
!         box_Ptemp(i_box) = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
!                                          curr_lon, curr_lat, curr_pressure)
!      endif


      DO Ki = 1,4,1

       !------------------------------------------------------------------
       ! For vertical wind speed:
       ! pay attention for the polar region * * *
       !------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       else
          curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       endif

       RK_omeg(Ki) = curr_omeg

       ! For 
!       dbox_lev = RK_Dt(Ki) * curr_omeg / 100.0     ! Pa => hPa
!       curr_pressure = box_lev(i_box) + dbox_lev

       RK_Dlev(Ki)   = Dt * curr_omeg / 100.0     ! Pa => hPa
       curr_pressure = box_lev(i_box) + RK_Dt(Ki+1) * curr_omeg / 100.0

       if(curr_pressure<P_mid(LLPAR)) &
             curr_pressure = P_mid(LLPAR) !+ ( P_mid(LLPAR) - curr_pressure )
       if(curr_pressure>P_mid(1)) &
             curr_pressure = P_mid(1) !- ( curr_pressure - P_mid(1) )


       !------------------------------------------------------------------
       ! For the region where lat<72, use Regualr Longitude-Latitude Mesh:
       !------------------------------------------------------------------
       if(abs(curr_lat)<=72.0)then

         curr_u = Interplt_wind_RLL(u, i_lon, i_lat, i_lev, curr_lon, &
                                                  curr_lat, curr_pressure)
         curr_v = Interplt_wind_RLL(v, i_lon, i_lat, i_lev, curr_lon, &
                                                  curr_lat, curr_pressure)

         RK_u(Ki) = curr_u
         RK_v(Ki) = curr_v

         dbox_lon  = (RK_Dt(Ki+1)*curr_u) &
                        / (2.0*PI*Re*COS(box_lat(i_box)*PI/180.0)) * 360.0
         dbox_lat  = (RK_Dt(Ki+1)*curr_v) / (PI*Re) * 180.0

         curr_lon  = box_lon(i_box) + dbox_lon
         curr_lat  = box_lat(i_box) + dbox_lat
 

         RK_Dlon(Ki) = (Dt*curr_u) &
                      / (2.0*PI*Re*COS(box_lat(i_box)*PI/180.0)) * 360.0
         RK_Dlat(Ki) = (Dt*curr_v) / (PI*Re) * 180.0

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

         RK_u(Ki) = curr_u_PS
         RK_v(Ki) = curr_v_PS

         dbox_x_PS = RK_Dt(Ki+1)*curr_u_PS
         dbox_y_PS = RK_Dt(Ki+1)*curr_v_PS


         RK_Dx_PS = Dt*curr_u_PS
         RK_Dy_PS = Dt*curr_v_PS


         !------------------------------------------------------------------
         ! change from (lon,lat) in RLL to (x,y) in PS: 
         !------------------------------------------------------------------
         if(box_lat(i_box)<0)then
           box_x_PS = -1.0* Re* COS(box_lon(i_box)*PI/180.0) &
                                / TAN(box_lat(i_box)*PI/180.0)
           box_y_PS = -1.0* Re* SIN(box_lon(i_box)*PI/180.0) &
                                / TAN(box_lat(i_box)*PI/180.0)
         else
           box_x_PS = Re* COS(box_lon(i_box)*PI/180.0) &
                        / TAN(box_lat(i_box)*PI/180.0)
           box_y_PS = Re* SIN(box_lon(i_box)*PI/180.0) &
                        / TAN(box_lat(i_box)*PI/180.0)
         endif

         RK_x_PS  = box_x_PS + RK_Dx_PS
         RK_y_PS  = box_y_PS + RK_Dy_PS

         box_x_PS  = box_x_PS + dbox_x_PS
         box_y_PS  = box_y_PS + dbox_y_PS


         !------------------------------------------------------------------
         ! change from (x,y) in PS to (lon,lat) in RLL
         !------------------------------------------------------------------
         if(box_x_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI 
         endif
         if(box_x_PS<0.0 .and. box_y_PS<=0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI -180.0
         endif
         if(box_x_PS<0.0 .and. box_y_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI +180.0
         endif
           
         if(curr_lat<0.0)then
           curr_lat= -1* ATAN( Re/SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         else
           curr_lat= ATAN( Re / SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         endif

         !------------------------------------------------------------------
         ! For 4th order Runge Kutta
         !------------------------------------------------------------------
         if(RK_x_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS<=0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI -180.0
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI +180.0
         endif

         if(RK_lat<0.0)then
           RK_lat = -1 * ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         else
           RK_lat = ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         endif

         RK_Dlon(Ki) = RK_lon - box_lon(i_box)
         RK_Dlat(Ki) = RK_lat - box_lat(i_box)

       endif ! if(abs(curr_lat)>72.0)then

      ENDDO ! Ki = 1,4,1

      box_lon(i_box) = box_lon(i_box) + &
                (RK_Dlon(1)+2.0*RK_Dlon(2)+2.0*RK_Dlon(3)+RK_Dlon(4))/6.0
      box_lat(i_box) = box_lat(i_box) + &
                (RK_Dlat(1)+2.0*RK_Dlat(2)+2.0*RK_Dlat(3)+RK_Dlat(4))/6.0
      box_lev(i_box) = box_lev(i_box) + &
                (RK_Dlev(1)+2.0*RK_Dlev(2)+2.0*RK_Dlev(3)+RK_Dlev(4))/6.0



      box_u    = ( RK_u(1) + 2.0*RK_u(2) &
                         + 2.0*RK_u(3) + RK_u(4) ) / 6.0
      box_v    = ( RK_v(1) + 2.0*RK_v(2) &
                         + 2.0*RK_v(3) + RK_v(4) ) / 6.0
      box_omeg = ( RK_omeg(1) + 2.0*RK_omeg(2) &
                         + 2.0*RK_omeg(3) + RK_omeg(4) ) / 6.0



!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '0-1D', i_box, &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,        &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp *       &
!                                        SUM(box_concnt_1D(i_box,:,1)),   &
!        Total_extra(i_box)




      !--------------------------------------------------------------------
      ! interpolate temperature for plume volumn change (PV=nRT):
      !--------------------------------------------------------------------
      if(abs(curr_lat)>Y_mid(JJPAR))then
         curr_T1 = Interplt_wind_RLL_polar(T1, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
      else
         curr_T1 = Interplt_wind_RLL(T1, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
      endif



      i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev(i_box),P_edge)



      if(abs(box_lat(i_box))>Y_mid(JJPAR))then
         next_T2 = Interplt_wind_RLL_polar(T2, i_lon, i_lat, i_lev, &
                           box_lon(i_box), box_lat(i_box), box_lev(i_box))
      else
         next_T2 = Interplt_wind_RLL(T2, i_lon, i_lat, i_lev, &
                           box_lon(i_box), box_lat(i_box), box_lev(i_box))
      endif


      ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
      ratio = ( (next_T2/box_lev(i_box))/(curr_T1/curr_pressure) )**(1/2)

      ! assume the volume change mainly apply to the cross-section, 
      ! the box_length would not change 


      IF(Judge_plume(i_box)==1) THEN
        V_prev = box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp

        box_Ra(i_box) = box_Ra(i_box) *ratio
        box_Rb(i_box) = box_Rb(i_box) *ratio

        V_new = box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp

        box_concnt_1D(:,:,i_box) = box_concnt_1D(:,:,i_box)*V_prev/V_new
!        Extra_mass_1D(i_box,:,:) = Extra_mass_1D(i_box,:,:)*V_prev/V_new
      ELSE IF(Judge_plume(i_box)==2) THEN
        V_prev = Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp

        Pdx(i_box) = Pdx(i_box) *ratio
        Pdy(i_box) = Pdy(i_box) *ratio

        V_new = Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp

        box_concnt_2D(:,:,:,i_box) = box_concnt_2D(:,:,:,i_box)*V_prev/V_new
!        Extra_mass_2D(i_box,:,:,:) = Extra_mass_2D(i_box,:,:,:)*V_prev/V_new
      ENDIF




!      IF(i_box==500.and.Judge_plume(i_box)==2) WRITE(6,*) 0, '2D', i_box, &
!         SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))*V_new, &
!         Total_extra(i_box)
!

!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '1-1D', i_box,&
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,        &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp *       &
!                                        SUM(box_concnt_1D(i_box,:,1)),   &
!        Total_extra(i_box)


!
!
!      IF(i_box==500.and.Judge_plume(i_box)==2) WRITE(6,*) '2D', V_new, &
!         box_concnt_2D(i_box,n_x_mid,n_y_mid,1), &
!                                box_concnt_2D(i_box,n_x_max,n_y_max,1)
!
!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '1D', &
!         box_concnt_1D(i_box,n_slab_50:n_slab_50+3,1)

!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '1D', &
!         box_concnt_1D(i_box,n_slab_max-3:n_slab_max,1)



      !------------------------------------------------------------------
      ! calcualte the box_alpha [0,2*PI)
      !------------------------------------------------------------------
      if((box_u**2+box_v**2)==0)then
          box_alpha(i_box) = 0.0
      else
        IF(box_v>=0)THEN
          box_alpha(i_box) = ACOS( box_u/SQRT(box_u**2+box_v**2) ) 
        ELSE
          box_alpha(i_box) = 2*PI - ACOS( box_u/SQRT(box_u**2+box_v**2) )
        ENDIF
      endif

      !------------------------------------------------------------------
      ! calcualte the Lyaponov exponent (Ly), unit: s-1
      !------------------------------------------------------------------
      IF(box_alpha(i_box)>=1.75*PI)THEN
        next_i_lon = i_lon + 1
        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR

        next_i_lat = i_lat
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1) next_i_lat=1


        D_wind = ABS(u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)) &
                +ABS(v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev))
        D_x    = Dx/360.0 * 2*PI *Re*COS(box_lat(i_box)/180*PI)
        Ly     = D_wind/D_x

      ELSEIF(box_alpha(i_box)>=1.25*PI)THEN
        next_i_lon = i_lon
        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR

        next_i_lat = i_lat - 1
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1) next_i_lat=1

        D_wind = ABS(u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)) &
                +ABS(v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev))
        D_y    = Dy/360.0 * 2*PI*Re
        Ly     = D_wind/D_y

      ELSEIF(box_alpha(i_box)>=0.75*PI)THEN
        next_i_lon = i_lon - 1
        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR

        next_i_lat = i_lat
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1) next_i_lat=1

        D_wind = ABS(u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)) &
                +ABS(v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev))
        D_x    = Dx/360.0 * 2*PI *Re*COS(box_lat(i_box)/180*PI) 
        Ly     = D_wind/D_x

      ELSEIF(box_alpha(i_box)>=0.25*PI)THEN
        next_i_lon = i_lon
        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR

        next_i_lat = i_lat + 1
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1) next_i_lat=1

        D_wind = ABS(u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)) &
                +ABS(v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev))
        D_y    = Dy/360.0 * 2*PI*Re
        Ly     = D_wind/D_y

      ELSE
        next_i_lon = i_lon + 1 
        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR

        next_i_lat = i_lat
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1) next_i_lat=1

        D_wind = ABS(u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)) &
                +ABS(v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev))
        D_x    = Dx/360.0 * 2*PI *Re*COS(box_lat(i_box)/180*PI)
        Ly     = D_wind/D_x

      ENDIF

      !------------------------------------------------------------------
      ! Horizontal stretch:
      ! Adjust the length/radius of box based on Lyaponov exponent (Ly)
      !------------------------------------------------------------------

!      IF(i_box==48)THEN
!      WRITE(6,*) '*** Horizontal strecth: ', box_lat(i_box), box_length(i_box)
!      WRITE(6,*) box_alpha, D_wind/Ly, D_wind, Ly
!      WRITE(6,*) i_lon, next_i_lon, i_lat, next_i_lat, Pdx(i_box), Pdy(i_box)
!
!      ENDIF


      length0           = box_length(i_box)
      box_length(i_box) = EXP(Ly*Dt) * box_length(i_box)


      IF(Judge_plume(i_box)==1) THEN
        box_Ra(i_box) = box_Ra(i_box)*SQRT(length0/box_length(i_box))
        box_Rb(i_box) = box_Rb(i_box)*SQRT(length0/box_length(i_box))
      ELSE IF(Judge_plume(i_box)==2) THEN
        Pdx(i_box) = Pdx(i_box)*SQRT(length0/box_length(i_box))
        Pdy(i_box) = Pdy(i_box)*SQRT(length0/box_length(i_box))

      ENDIF



      n_active_plume = n_active_plume + 1


!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '1-1D', i_box,&
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,        &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp *       &
!                                        SUM(box_concnt_1D(i_box,:,1)),   &
!        Total_extra(i_box)




500     CONTINUE

    end do  !do i_box = 1,n_box

!    WRITE(6,*) '--- n_active_plume, n_box:',  n_active_plume, n_box

    !------------------------------------------------------------------
    ! Horizontal stretch:
    ! Adjust the length/radius of the box based on new location
    !------------------------------------------------------------------

!    DO i_box = 1, n_box
! 
!      length0           = box_length(i_box)
!      box_length(i_box) = 2 * Distance_Circle( box_lon(i_box), box_lat(i_box), &
!                                   box_lon_edge(i_box), box_lat_edge(i_box) ) ! [m]
!
!      IF(Judge_plume(i_box)==1) THEN
!        box_Ra(i_box) = box_Ra(i_box)*SQRT(length0/box_length(i_box))
!        box_Rb(i_box) = box_Rb(i_box)*SQRT(length0/box_length(i_box))
!      ELSE IF(Judge_plume(i_box)==2) THEN
!         Pdx(i_box) = Pdx(i_box)*SQRT(length0/box_length(i_box))
!         Pdy(i_box) = Pdy(i_box)*SQRT(length0/box_length(i_box))
!
!      ENDIF
!
!    ENDDO 


    !------------------------------------------------------------------
    ! Everything is done, clean up pointers
    !------------------------------------------------------------------
    nullify(u)
    nullify(v)
    nullify(omeg)

    nullify(Ptemp)
    nullify(T1)
    nullify(T2)

    nullify(X_edge)
    nullify(Y_edge)

    nullify(PASV_EU)


  END SUBROUTINE lagrange_run


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real(fp) function Interplt_wind_RLL(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: wind_lonlat(2), wind_lonlat_lev

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

       Weight(i,j) = 1.0/distance(i,j) / SUM( 1.0/distance(:,:) )

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

  real(fp) function Interplt_wind_RLL_polar(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
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

  real(fp) function Interplt_uv_PS_polar(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

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
      curr_x = -1.0* Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = -1.0* Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    else
      curr_x = Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
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
        x_PS(i)= Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
        y_PS(i)= Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
      else
        x_PS(i)= -1.0* Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
        y_PS(i)= -1.0* Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
      endif


       do k=1,2
          kk = k + init_lev - 1
          IF(i_uv==1)THEN ! i_ux==1 for u
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_PS(i,k) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                        + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF

          IF(i_uv==0)THEN ! for v
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                       - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          else
            uv_PS(i,k) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                           + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
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
          uv_polars(ii) = -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                         / SIN(Y_mid(jj)*PI/180.0) &
                         + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                         / (SIN(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_polars(ii) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
        enddo

            uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF

          IF(i_uv==0)THEN ! for v
          do ii = 1,IIPAR
             if(Y_mid(jj)>0)then
               uv_polars(ii) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                             - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
             else
               uv_polars(ii) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                                 + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
             endif
          enddo
             uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF
    enddo



    ! calculate the distance between particle and grid point
    do i = 1,3
       distance_PS(i)= SQRT( (x_PS(i)-curr_x)**2.0 + (y_PS(i)-curr_y)**2.0 )
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

  real(fp) function Interplt_uv_PS(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

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
      curr_x = -1.0* Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = -1.0* Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    else
      curr_x = Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
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

        x_PS(i,j) = Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)  
        y_PS(i,j) = Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
          
        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then ! i_ux==1 for u
             uv_PS(i,j,k)= -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                         / SIN(Y_mid(jj)*PI/180.0) &
                                + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2) )
           endif

           if(i_uv==0)then ! for v
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      else

      x_PS(i,j)= -1.0* Re* COS(X_mid(ii)*PI/180.0) /TAN(Y_mid(jj)*PI/180.0)
      y_PS(i,j)= -1.0* Re* SIN(X_mid(ii)*PI/180.0) /TAN(Y_mid(jj)*PI/180.0)

        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif

           if(i_uv==0)then
             uv_PS(i,j,k) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                              + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      endif

    enddo
    enddo

    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
       distance_PS(i,j) = SQRT( (x_PS(i,j)-curr_x)**2.0 &
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

    USE TIME_MOD              ! For computing date & time 
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_EXIT

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR

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

    REAL(fp) :: Dt          ! = 600.0e+0_fp          
    REAL(fp) :: Pdt
    REAL(fp) :: Dt2

    integer  :: n_box_max
    integer  :: N_split

    integer  :: i_box, i_specie
    integer  :: ii_box
    integer  :: i_x, i_y, i, j
    integer  :: i_slab

    integer  :: i_lon, i_lat, i_lev

    integer  :: t1s

    real(fp) :: curr_lon, curr_lat, curr_pressure
!    real(fp) :: next_lon, next_lat, next_pressure

    real(fp) :: X_edge2, Y_edge2

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    real(fp), pointer :: Ptemp(:,:,:)

    real(fp), pointer :: P_BXHEIGHT(:,:,:)

!    real(fp)          :: curr_u, curr_v, curr_omeg


    real(fp)          :: curr_Ptemp, Ptemp_shear  ! potential temperature
    real(fp)          :: U_shear, V_shear, UV_shear

    real(fp)          :: Dx, Dy
    real(fp)          :: wind_s_shear
    real(fp)          :: theta_previous

    ! angle between travel direction and lon: beta belongs to [-PI/2,3*PI/2)

    ! same as angle between cross-section (right when 
    ! toward travel direction) and lon: box_alpha belongs to [0, 2*PI)
!    real(fp)          :: box_alpha ! beta! [radians]

    real(fp), pointer :: X_edge(:)
    real(fp), pointer :: Y_edge(:)


    real(fp)  :: Outer(n_slab_max), Inner(n_slab_max)


    real(fp)  :: eddy_v, eddy_h
    real(fp)  :: eddy_A, eddy_B
    real(fp)  :: Cv, Ch, Omega_N, N_BV

    real(fp)  :: grid_volumn
    real(fp)  :: V_grid_1D, V_grid_2D
    real(fp)  :: mass_plume, mass_plume_new
    real(fp)  :: D_mass_plume, D_mass_extra ! mass change in plume and extra
                                            ! model every time step [molec]

    real(fp)  :: Ratio_radius, Ra_Rb

    real(fp)  :: L_b, L_O, Ee

    real(fp)  :: CFL
    real(fp)  :: CFL_2d(n_x_max,n_y_max) ! for 2D advection

    real(fp)  :: Pu(n_x_max,n_y_max) ! for 2D advection

    real(fp)  :: Pc(n_x_max,n_y_max),Pc2(n_x_max,n_y_max) !, Ec(n_x_max,n_y_max)
    real(fp)  :: Pc_bdy(n_x_max2,n_y_max2)
!    real(fp)  :: D_concnt(n_x_max,n_y_max)
    real(fp)  :: Xscale, Yscale, frac_mass
    real(fp)  :: C2d_prev(n_x_max,n_y_max), C2d_prev_extra(n_x_max,n_y_max)

    real(fp)  :: Cslab(n_slab_max) !, Extra_Cslab(n_slab_max)

    real(fp)  :: Concnt2D_bdy(n_x_max2, n_y_max2)
    real(fp)  :: Concnt1D_bdy(n_slab_max2)

    real(fp)  :: C1d_prev(n_slab_max), C1d_prev_extra(n_slab_max)

    real(fp)  :: DlonN, DlatN ! split 1D from 1 to 5 segments

    real(fp)  :: backgrd_concnt(N_specie)

    real(fp)  :: Rate_Mix, Rate_Eul, Eul_concnt


    real(fp)  :: start, finish


!    CHARACTER(LEN=255)     :: FILENAME2
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg


    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR

    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()


    RC        =  GC_SUCCESS
    ErrMsg    =  ''

!    FILENAME2   = 'Plume_theta_max_min_radius.txt'

!    n_box_max = n_box

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
    ! Convert specie to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90')
       RETURN
    ENDIF


!    if(N_curr<=n_box)then
!      n_box = N_curr
!    else
!      n_box = n_box
!    endif

    !======================================================================
    ! Set up the new injected plumes/box
    ! identify the location of the plume box, 
    ! put the corresponding concentration of backgound grid into the initial 
    ! concentration inside plume in unit of [molec/cm3].
    !======================================================================
    IF(n_box_prev==N_parcel) THEN
    DO i_box = 1, N_parcel, 1       ! Only for initail injected box

       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)   ! hPa


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       if(i_lon>IIPAR) i_lon=i_lon-IIPAR
       if(i_lon<1) i_lon=i_lon+IIPAR

       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
       if(i_lat>JJPAR) i_lat=JJPAR
       if(i_lat<1) i_lat=1

       i_lev = Find_iPLev(curr_pressure,P_edge)


       do i_specie = 1, N_specie
         box_concnt_2D(:,:,i_specie,i_box) = box_concnt_2D(:,:,i_specie,i_box) &
                + State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)

         ! [molec]
!         Extra_mass_2D(i_box,:,:,i_specie) = &
!                        State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)

         ! [molec]
!         Total_extra(i_box) = &
!                        Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
!                      * SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))

         Total_extra(i_box) = &
              Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
           * State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) &
           * (n_x_max-2)*(n_y_max-2)
       enddo



    ENDDO
    ENDIF !IF(n_box_prev==N_parcel) THEN



    DO i_box = n_box_prev+1, n_box, 1       ! Only for new injected box

       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)   ! hPa


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       if(i_lon>IIPAR) i_lon=i_lon-IIPAR
       if(i_lon<1) i_lon=i_lon+IIPAR

       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
       if(i_lat>JJPAR) i_lat=JJPAR
       if(i_lat<1) i_lat=1

       i_lev = Find_iPLev(curr_pressure,P_edge)


       do i_specie = 1, N_specie
         box_concnt_2D(:,:,i_specie,i_box) = box_concnt_2D(:,:,i_specie,i_box) &
                + State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
        
         ![molec]
!         Extra_mass_2D(i_box,:,:,i_specie) = &
!                 State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)

         ! [molec]
!         Total_extra(i_box) = &
!                        Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
!                      * SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))


         Total_extra(i_box) = &
              Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
           * State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) &
           * (n_x_max-2)*(n_y_max-2)



       enddo 

    ENDDO ! DO i_box = n_box_prev+1, n_box_max, 1


    !=====================================================================
    ! Run Plume distortion & dilution HERE
    !=====================================================================
    n_box_max = n_box
    do i_box = 1,n_box_max

         
       IF(Judge_plume(i_box)==0) GOTO 200

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

       ! make sure the location is not out of range
!       do while (box_lat_edge(i_box) > Y_edge(JJPAR+1))
!          box_lat_edge(i_box) = Y_edge(JJPAR+1) &
!                         - ( box_lat_edge(i_box)-Y_edge(JJPAR+1) )
!       end do
!
!       do while (box_lat_edge(i_box) < Y_edge(1))
!          box_lat_edge(i_box) = Y_edge(1) + ( box_lat_edge(i_box)-Y_edge(1) )
!       end do

!       do while (box_lon_edge(i_box) > X_edge(IIPAR+1))
!          box_lon_edge(i_box) = box_lon_edge(i_box) - 360.0
!       end do

!       do while (box_lon_edge(i_box) < X_edge(1))
!          box_lon_edge(i_box) = box_lon_edge(i_box) + 360.0
!       end do


       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)        ! hPa

!       curr_u        = box_u(i_box)
!       curr_v        = box_v(i_box)
!       curr_omeg     = box_omeg(i_box)
!       curr_Ptemp    = box_Ptemp(i_box)


       !------------------------------------------------------------------
       ! calcualte the box_alpha
       ! Next adjacent point
       !------------------------------------------------------------------
!       next_lon      = box_lon_edge(i_box)
!       next_lat      = box_lat_edge(i_box)
!       next_pressure = box_lev_edge(i_box)        ! hPa

       ! Judge the sign of box_alpha
!       if((next_lon-curr_lon)>0.0)then
!         beta = ATAN( (next_lat - curr_lat) / (next_lon - curr_lon) )
!         ! (-PI/2,PI/2)
!       else if((next_lon-curr_lon)<0.0)then
!         beta = PI + ATAN( (next_lat - curr_lat) / (next_lon - curr_lon) )
!         ! (PI/2, 3*PI/2)
!       else if((next_lat-curr_lat)>0.0)then
!         beta = 0.5*PI
!       else if((next_lat-curr_lat)<0.0)then
!         beta = -0.5*PI
!       endif
!       ! beta is [-PI/2,3*PI/2), box_alpha is [-PI, PI)
!       box_alpha = beta - PI/2

!      if((box_u**2+box_v**2)==0)then
!          box_alpha = 0.0
!      else
!        IF(box_v(i_box)>=0)THEN
!          box_alpha = ACOS( box_u/SQRT(box_u**2+box_v**2) )
!        ELSE
!          box_alpha = 2*PI-ACOS( box_u/SQRT(box_u**2+box_v**2) )
!        ENDIF
!      endif


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       if(i_lon>IIPAR) i_lon=i_lon-IIPAR
       if(i_lon<1) i_lon=i_lon+IIPAR

       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
       if(i_lat>JJPAR) i_lat=JJPAR
       if(i_lat<1) i_lat=1

       i_lev = Find_iPLev(curr_pressure,P_edge)


!       backgrd_concnt(i_box,:) = State_Chm%Species(i_lon,i_lat,i_lev,:)
       backgrd_concnt(1) = State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
       ! [molec/cm3]

!       Plume_I(i_box) = i_lon
!       Plume_J(i_box) = i_lat
!       Plume_L(i_box) = i_lev


       !====================================================================
       ! calculate the wind shear along plume corss-section
       ! calculate the diffusivity in horizontal and vertical direction
       !====================================================================
       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha(i_box), i_lon, &
                        i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! *** attention *** ???
       wind_s_shear = 0.002


       ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
       Cv = 0.2
       Omega_N = 0.1
       Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, i_lon, i_lat, &
                             i_lev,curr_lon, curr_lat, curr_pressure)


       !--------------------------------------------------------------------
       ! interpolate potential temperature for Plume module:
       !--------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
         curr_Ptemp = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
       else
         curr_Ptemp = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
       endif


       N_BV = SQRT(Ptemp_shear*g0/curr_Ptemp)
       IF(N_BV<=0.001) N_BV = 0.001

       ! diffusivity unit: [m2/s]
       eddy_v = Cv * Omega_N**2 / N_BV
       eddy_v = 1.0 ! ???

       ! Calculate horizontal eddy diffusivity:
!       Ch = 0.1
!       U_shear = Vertical_shear(u, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon,
!       curr_lat, curr_pressure)
!       V_shear = Vertical_shear(v, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon,
!       curr_lat, curr_pressure)
!       UV_shear = SQRT( U_shear**2 + V_shear**2 )
!       eddy_h(i_slab) = Ch*UV_shear*(Init_radius+(i_slab-1)*D_radius)**2
       eddy_h = 10.0  ! ???

! test for horizontal diffusivity
!       L_b = 2.0*PI*SQRT(curr_u**2+curr_v**2)/N_BV
!       Ee = 0.5 * Omega_N**2 * wind_s_shear**2
!       L_O = 2.0*PI*SQRT( Ee/(N_BV**3) )


       !====================================================================
       ! Now choose next step based on Judge_plume
       ! 0: the plume already dissolved into background grids, skip this box
       ! 1: use 1D slab model to simulate the plume corss-section
       ! 2: use 2D model to simulate the plume cross-section
       !====================================================================

       IF(Judge_plume(i_box)==2) THEN

!       call cpu_time(start)


       ! Define the wind field based on wind shear
       DO i_y = 1, n_y_max
         Pu(:,i_y) = (i_y-n_y_mid)*Pdy(i_box) * wind_s_shear
         ! [m s-1]
       ENDDO


       !-------------------------------------------------------------------
       ! Calculate the advection-diffusion in 2D grids
       !-------------------------------------------------------------------

       Pdt = Dt*0.5*0.5*0.5

       ! Find the best Pdt to meet CFL condition:
 700     CONTINUE

       Pdt = 0.5*Pdt

       CFL = Pdt*Pu(2,n_y_max)/Pdx(i_box)
       IF(ABS(CFL)>1) GOTO 700
       IF(ABS(2*eddy_h*Pdt/(Pdx(i_box)**2))>1) GOTO 700
       IF(ABS(2*eddy_v*Pdt/(Pdy(i_box)**2))>1) GOTO 700



       DO t1s = 1, NINT(Dt/Pdt)


         !--------------------------------------------------------------------
         ! update the boundary grids in 2D model based on background 
         ! contentration
         ! background should be updated every time step
         ! shw !!!
         !--------------------------------------------------------------------


!         box_concnt_2D(i_box,1,:,1)       = backgrd_concnt(1)
!         box_concnt_2D(i_box,n_x_max,:,1) = backgrd_concnt(1)
!         box_concnt_2D(i_box,:,1,1)       = backgrd_concnt(1)
!         box_concnt_2D(i_box,:,n_y_max,1) = backgrd_concnt(1)


!         Extra_mass_2D(i_box,1,:,1)       = backgrd_concnt(1)
!         Extra_mass_2D(i_box,n_x_max,:,1) = backgrd_concnt(1)
!         Extra_mass_2D(i_box,:,1,1)       = backgrd_concnt(1)
!         Extra_mass_2D(i_box,:,n_y_max,1) = backgrd_concnt(1)

        
!         IF(i_box==1.and.t1s==1) &
!                  WRITE(6,*)'Judge_plume(i_box), Pdx(i_box), Pdy(i_box):', &
!                                      Judge_plume(i_box), Pdx(i_box), Pdy(i_box)


         !--------------------------------------------------------------
         ! if Pdx become smaller than half Dx_initm,
         ! combine 9 grids into 1 grids. 
         ! Update: box_concnt_2D(), Pdx(), Pdy(),  Extra_mass_2D()
         !--------------------------------------------------------------
         IF(Pdx(i_box)<=0.5*Dx_init)THEN

600     CONTINUE

         
           Pc = box_concnt_2D(:,:,1,i_box) ! [molec cm-3]
!           Ec = Extra_mass_2D(i_box,:,:,1) ! [molec cm-3]

           box_concnt_2D(:,:,1,i_box)       = backgrd_concnt(1)
!           Extra_mass_2D(i_box,:,:,1)       = backgrd_concnt(1)


           DO i = 1, n_x_max/3, 1
           DO j = 1, n_y_max/3, 1

             i_x = (i-1)*3+1
             i_y = (j-1)*3+1

             box_concnt_2D(i+n_x_max/3,j+n_y_max/3,1,i_box) = &
                                      SUM(Pc(i_x:i_x+2,i_y:i_y+2))/9

!             Extra_mass_2D(i_box,i+n_x_max/3,j+n_y_max/3,1) = &
!                                      SUM(Ec(i_x:i_x+2,i_y:i_y+2))/9
! ??? Total_extra(i_box) = 
           ENDDO
           ENDDO

           Pdx(i_box) = Pdx(i_box)*3
           Pdy(i_box) = Pdy(i_box)*3


           ! uodate Total_extra(i_box)
           V_grid_2D       = Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp
           mass_plume      = V_grid_2D/9 *SUM(Pc(:,:))

           mass_plume_new = V_grid_2D &
                *SUM( box_concnt_2D(:,:,1,i_box))

           D_mass_plume = mass_plume_new - mass_plume

           IF(D_mass_plume>0)THEN
             WRITE(6,*)'--- 9 to 1 grid cell:', i_box, mass_plume, mass_plume_new
             Total_extra(i_box) = Total_extra(i_box) + D_mass_plume
           ELSE
             ! WRITE(6,*)'--- ERROR in 9 to 1 grid cell in 2D grid ---'
           ENDIF


!           WRITE(6,*)'--- Check combining 9 to 1 grid cells in 2D:'
!           WRITE(6,*) Pdx(i_box), Pdy(i_box)
!                      SUM( box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)  &
!                          -Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)) &
!                      *Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp,  &
!                      SUM(Pc-Ec)*Pdx(i_box)*Pdy(i_box)/9                   &
!                                              *box_length(i_box)*1.0e+6_fp



           IF(Pdx(i_box)<=0.5*Dx_init) WRITE(6,*) &
               "ERROR 0: more combination: ", i_box, Pdx(i_box), Pdy(i_box)

           IF(Pdx(i_box)<=0.5*Dx_init) GOTO 600


         ENDIF ! IF(Pdx(i_box)<0.5*Dx_init)THEN


!         CFL = Pdt*Pu(2,n_y_max)/Pdx(i_box)
!         IF(ABS(CFL)>1) WRITE(6,*) '* ERROR 1 *', Pu(2,n_y_max), Pdx(i_box)
!
!         IF(ABS(2*eddy_h*Pdt/(Pdx(i_box)**2))>1) &
!                WRITE(6,*) '*** ERROR 2 ***', eddy_h, Pdt, Pdx(i_box)
!
!         IF(ABS(2*eddy_v*Pdt/(Pdy(i_box)**2))>1) &
!                WRITE(6,*) '*** ERROR 3 ***', eddy_v, Pdt, Pdy(i_box)



         !--------------------------------------------------------------
         ! calculate [vertical wind shear advection] and [diffusion]
         !--------------------------------------------------------------

         ! for interaction with background grid cell
         C2d_prev           = box_concnt_2D(:,:,1,i_box)
!         C2d_prev_extra     = Extra_mass_2D(i_box,:,:,1)

         V_grid_2D       = Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp


!         IF(i_box==1.and.t1s==1) WRITE(6,*) 1, i_box, Lifetime(i_box), &
!           V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)),
!           Total_extra(i_box)


         Concnt2D_bdy(:,:) = backgrd_concnt(1)
         Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1) = box_concnt_2D(:,:,1,i_box)


         ! advection

!         Pc = box_concnt_2D(i_box,:,:,1) ! [molec cm-3]
!         Ec = Extra_mass_2D(i_box,:,:,1) ! [molec cm-3]
         Pc_bdy = Concnt2D_bdy

         CFL_2d = Pdt*Pu(:,:)/Pdx(i_box) ! [n_x_max, n_y_max]

         Concnt2D_bdy(2:n_x_max2-1, 2:n_y_max2-1) = &
                Pc_bdy(2:n_x_max2-1,2:n_y_max2-1) &
                - 0.5 *CFL_2d(:,:) *( Pc_bdy(3:n_x_max2,2:n_y_max2-1) &
                             -Pc_bdy(1:n_x_max2-2,2:n_y_max2-1) ) &
                + 0.5 *CFL_2d(:,:)**2 *(   Pc_bdy(3:n_x_max2,2:n_y_max2-1) &
                                -2*Pc_bdy(2:n_x_max2-1,2:n_y_max2-1) &
                                  +Pc_bdy(1:n_x_max2-2,2:n_y_max2-1) )


!         Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1) = &
!                Ec(2:n_x_max-1,2:n_y_max-1) &
!                - 0.5 *CFL_2d(:,:) *( Ec(3:n_x_max,2:n_y_max-1) &
!                             -Ec(1:n_x_max-2,2:n_y_max-1) ) &
!                + 0.5 *CFL_2d(:,:)**2 *(   Ec(3:n_x_max,2:n_y_max-1) &
!                                -2*Ec(2:n_x_max-1,2:n_y_max-1) &
!                                  +Ec(1:n_x_max-2,2:n_y_max-1) )


!         DO i_y = 2, n_y_max-1, 1
!           CFL = Pdt*Pu(2,i_y)/Pdx(i_box) ! [1]
!
!           box_concnt_2D(i_box,2:n_x_max-1,i_y,1) = Pc(2:n_x_max-1,i_y) &
!            - 0.5 *CFL *(Pc(3:n_x_max,i_y)-Pc(1:n_x_max-2,i_y))         &
!            + 0.5 *CFL**2 *(Pc(3:n_x_max,i_y)-2*Pc(2:n_x_max-1,i_y)     &
!                                               +Pc(1:n_x_max-2,i_y))
!
!
!           Extra_mass_2D(i_box,2:n_x_max-1,i_y,1) = Ec(2:n_x_max-1,i_y) &
!            - 0.5 *CFL *(Ec(3:n_x_max,i_y)-Ec(1:n_x_max-2,i_y))         &
!            + 0.5 *CFL**2 *(Ec(3:n_x_max,i_y)-2*Ec(2:n_x_max-1,i_y)     &
!                                               +Ec(1:n_x_max-2,i_y))
!
!         ENDDO


!         IF(i_box==2095.and.t1s==1) WRITE(6,*) 2, i_box, &
!           V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)), &
!           V_grid_2D*SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))


         ! diffusion
         Pc_bdy = Concnt2D_bdy

         Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1) =               &
                Pc_bdy(2:n_x_max2-1,2:n_y_max2-1)                             &
           + Pdt*( eddy_h*(   Pc_bdy(1:n_x_max2-2,2:n_y_max2-1)               &
                           -2*Pc_bdy(2:n_x_max2-1,2:n_y_max2-1)               &
                             +Pc_bdy(3:n_x_max2,2:n_y_max2-1)  )              &
                                                /(Pdx(i_box)**2)        &
                 + eddy_v*(   Pc_bdy(2:n_x_max2-1,1:n_y_max2-2)               &
                           -2*Pc_bdy(2:n_x_max2-1,2:n_y_max2-1)               &
                             +Pc_bdy(2:n_x_max2-1,3:n_y_max2)  )              &
                                                /(Pdy(i_box)**2) )


!         Ec = Extra_mass_2D(i_box,:,:,1) ! [molec cm-3]

!         Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1) = &
!                Ec(2:n_x_max-1,2:n_y_max-1)               &
!           + Pdt*(eddy_h*(   Ec(1:n_x_max-2,2:n_y_max-1)  &
!                          -2*Ec(2:n_x_max-1,2:n_y_max-1)  &
!                            +Ec(3:n_x_max,2:n_y_max-1)  )/(Pdx(i_box)**2) &
!                 +eddy_v*(   Ec(2:n_x_max-1,1:n_y_max-2)  &
!                          -2*Ec(2:n_x_max-1,2:n_y_max-1)  &
!                            +Ec(2:n_x_max-1,3:n_y_max)  )/(Pdy(i_box)**2))


!         IF(i_box==2095.and.t1s==1) WRITE(6,*) 3, i_box, &
!           V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)), &
!           V_grid_2D*SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))

         box_concnt_2D(:,:,1,i_box) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)

         !================================================================
         ! Update the concentration in the background grid cell
         ! after the interaction with 2D plume
         !================================================================
         grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
         V_grid_2D       = Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp

!         D_mass_extra = V_grid_2D &
!                *SUM( Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)      &
!                    - C2d_prev_extra(2:n_x_max-1,2:n_y_max-1) )
         ! [molec]


         ! the boundary always represents the background concentration
         ! [molec]
         mass_plume      = V_grid_2D*SUM(C2d_prev(:,:))

         mass_plume_new = V_grid_2D &
                *SUM( box_concnt_2D(:,:,1,i_box))


         D_mass_plume = mass_plume_new - mass_plume




         IF(D_mass_plume<0)THEN
           ! generally D_mass_plume is negative
           D_mass_plume = (1-Total_extra(i_box)/mass_plume)*D_mass_plume


           Total_extra(i_box) = Total_extra(i_box) &
                           + Total_extra(i_box)/mass_plume *D_mass_plume
         ENDIF


         backgrd_concnt(1) = ( backgrd_concnt(1)*grid_volumn &
                                          - D_mass_plume) /grid_volumn

         State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = &
                                                backgrd_concnt(1)


!         box_concnt_2D(i_box,1,:,1)       = backgrd_concnt(1)
!         box_concnt_2D(i_box,n_x_max,:,1) = backgrd_concnt(1)
!         box_concnt_2D(i_box,:,1,1)       = backgrd_concnt(1)
!         box_concnt_2D(i_box,:,n_y_max,1) = backgrd_concnt(1)


!         IF(i_box==1.and.t1s==1)THEN
!           WRITE(6,*)'--- Check advection-diffusion in 2D:',            &
!                     SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1) &
!                       -Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)) &
!                     *Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp,&
!                    D_mass_plume, D_mass_extra
!         ENDIF



         !====================================================================
         ! Change from 2D to 1D, 
         ! once the tilting degree is bigger than 88 deg (88/180*3.14)
         !====================================================================
         

         IF(Lifetime(i_box)>3.0*3600.0)THEN

!         D_concnt  = box_concnt_2D(i_box,:,:,1) - Extra_mass_2D(i_box,:,:,1)
         frac_mass = 0.95

!         Xscale = Get_XYscale(D_concnt, i_box, frac_mass, 2)
!         Yscale = Get_XYscale(D_concnt, i_box, frac_mass, 1)

         Pc2 = box_concnt_2D(:,:,1,i_box)
         Xscale = Get_XYscale(Pc2, i_box, frac_mass, 2)
         Yscale = Get_XYscale(Pc2, i_box, frac_mass, 1)

         box_theta(i_box) = ATAN( Xscale/Yscale )


!         IF(i_box==1.and.t1s==1)THEN
!           WRITE(6,*)'i_box, box_theta(i_box), 87.0/180.0*PI: ', &
!                               i_box, box_theta(i_box), 87.0/180.0*PI
!         ENDIF

         IF( box_theta(i_box)>(87.0/180.0*PI) &
                        .and. Lifetime(i_box)>6.0*3600.0 ) THEN
           WRITE(6,*) "*** ERROR ***"
           WRITE(6,*) i_box, Xscale, Yscale, box_theta(i_box)
         ENDIF


         IF( box_theta(i_box)>(87.0/180.0*PI) &
                        .or. Lifetime(i_box)>6.0*3600.0 ) THEN !!! shw ???


           Pc = box_concnt_2D(:,:,1,i_box)
!           Ec = Extra_mass_2D(i_box,:,:,1)

           ! calculate accutate box_theta
!           box_theta(i_box) = &
!             Find_theta(box_theta(i_box), Pc, Yscale, i_box)

!           box_theta(i_box) = 1.521235 !!! shw

           IF(box_theta(i_box)>1.56) WRITE(6,*)&
            '*** ERROR *** i_box, box_theta(i_box):', i_box, box_theta(i_box),&
                                                                Lifetime(i_box)

           ! assign length/width and concentration to 1D slab model 
           ! return initial box_Ra(i_box), box_Rb(i_box), box_concnt_1D(i_box)
           CALL Slab_init(box_theta(i_box), backgrd_concnt(1), i_box, Pc, Yscale)


           Judge_plume(i_box)=1
           V_grid_1D = box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp

           mass_plume_new = V_grid_1D &
                            * SUM(box_concnt_1D(1:n_slab_max,1,i_box))

           mass_plume = V_grid_2D &
                        * SUM(box_concnt_2D(:,:,1,i_box))

           !!! shw
!           Total_extra(i_box) = &
!                        V_grid_1D * SUM(Extra_mass_1D(i_box,2:n_slab_max-1,1))


           !!! update the background concentration:
           D_mass_plume = mass_plume_new - mass_plume


!           D_mass_extra = V_grid_1D * SUM(Extra_mass_1D(i_box,2:n_slab_max-1,1)) &
!                 - V_grid_2D * SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))

           IF(D_mass_plume<0)THEN
             ! generally D_mass_plume is negative
             D_mass_plume = (1-Total_extra(i_box)/mass_plume)*D_mass_plume


             Total_extra(i_box) = Total_extra(i_box) &
                             + Total_extra(i_box)/mass_plume *D_mass_plume
           ENDIF


!           IF(i_box==2095.and.t1s==1) WRITE(6,*) 4, i_box, &
!                 V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)), &
!                 V_grid_2D*SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))



           backgrd_concnt(1) = ( backgrd_concnt(1)*grid_volumn &
                                                - D_mass_plume ) /grid_volumn

           State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = &
                                                    backgrd_concnt(1)



           GOTO 400 ! begin 1D slab model

         ENDIF ! IF(box_theta(i_box)>98/180*3.14)
         ENDIF ! IF(Lifetime(i_box)>3*3600.0)THEN

         !===================================================================
         ! Third judge:
         ! If this is the final time step, the model is going to end
         ! dissolve all the plume into GCM in the last time step
         !===================================================================
!         IF ( ITS_TIME_FOR_EXIT() ) THEN
         IF(DAY==31 .and. HOUR==23)THEN

!           backgrd_concnt(1) = &
!              ( SUM(  box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)     &
!                    - Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)  )  &
!                      *V_grid_2D   &
!                  + backgrd_concnt(1)*grid_volumn ) / grid_volumn

           backgrd_concnt(1) = &
            ( SUM( box_concnt_2D(:,:,1,i_box) )  &
               *V_grid_2D - Total_extra(i_box)  &
                  + backgrd_concnt(1)*grid_volumn ) / grid_volumn


           State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) =    &
                                                       backgrd_concnt(1)

           Judge_plume(i_box) = 0
           GOTO 200 ! skip this box, go to next box

         ENDIF

         !===================================================================
         ! Fourth judge:
         ! If the plume touch the tropopause, dissolve the plume
         !===================================================================
         IF ( box_lev(i_box)>State_Met%TROPP(i_lon,i_lat) ) THEN

!           backgrd_concnt(1) = &
!              ( SUM(  box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)     &
!                    - Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)  )  &
!                      *V_grid_2D   &
!                  + backgrd_concnt(1)*grid_volumn ) / grid_volumn


           backgrd_concnt(1) = &
            ( SUM( box_concnt_2D(:,:,1,i_box) )  &
               *V_grid_2D - Total_extra(i_box)  &
                  + backgrd_concnt(1)*grid_volumn ) / grid_volumn


           State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) =    &
                                                       backgrd_concnt(1)

           Judge_plume(i_box) = 0
           GOTO 200 ! skip this box, go to next box

         ENDIF



!       IF(i_box==2095.and.t1s==1) WRITE(6,*) 5, i_box, &
!         V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)), &
!         V_grid_2D*SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))




       ENDDO ! DO t1s = 1, NINT(Dt/Pdt)
        
        
!       IF(i_box==2095) WRITE(6,*) 5.5, i_box, &
!         V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)), &
!         V_grid_2D*SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))



!       call cpu_time(finish)
!       WRITE(6,*)'Time (finish-start) for 2D:', i_box, finish-start


       ENDIF ! IF(Judge_plume(i_box)==2) THEN


       ! *******************************************************************
       !====================================================================
       ! begin 1D slab model
       !====================================================================
400     CONTINUE

       IF(Judge_plume(i_box)==1) THEN


!       call cpu_time(start)


       ! ============================================
       ! Decide the time step based on CFL condition
       ! 2*k*Dt/Dr<1 (or Dr-2*k*Dt>0) for diffusion
       !==================================================================

       ! ignore the diffusivity in long radius direction
       eddy_B = eddy_v*SIN(abs(box_theta(i_box))) &
                          + eddy_h*COS(box_theta(i_box)) ! b

       Dt2 = Dt

       IF((2*eddy_B*Dt2/(box_Rb(i_box)**2))>0.5) THEN
          Dt2 = Dt*0.1
       ENDIF


 300     CONTINUE


       IF((2*eddy_B*Dt2/(box_Rb(i_box)**2))>0.5) THEN

       

!!! shw put this as a function


         !-------------------------------------------------------------------
         ! Decrease half of slab by combining 2 slab into 1 slab
         !-------------------------------------------------------------------
         Cslab       = box_concnt_1D(:,N_specie,i_box)
!         Extra_Cslab = Extra_mass_1D(i_box,:,N_specie)


         DO i_slab = 1, n_slab_50
           box_concnt_1D(i_slab+N_slab_25,N_specie,i_box) = &
             (  Cslab(i_slab*2-1) + Cslab(i_slab*2) ) /2

!           Extra_mass_1D(i_box,i_slab+N_slab_25,N_specie) = &
!             (  Extra_Cslab(i_slab*2-1) + Extra_Cslab(i_slab*2) ) /2

         ENDDO

         box_Rb(i_box) = box_Rb(i_box)*2

!         DO i_specie = 1,N_specie
!           Extra_mass(i_box,i_specie) = &
!                      SUM(Extra_mass(i_box,i_specie))
!            ![molec]
!         ENDDO


         !-------------------------------------------------------------------
         ! Meanwhile, add new slab into plume to compensate lost slab
         !-------------------------------------------------------------------
         DO i_slab = 1,n_slab_25
         do i_specie = 1,N_specie
           box_concnt_1D(i_slab,i_specie,i_box) = &
                   State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)

!           Extra_mass_1D(i_box,i_slab,i_specie) = &
!                   State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
         enddo
         ENDDO


         DO i_slab = n_slab_75+1, n_slab_max
         do i_specie = 1,N_specie
           box_concnt_1D(i_slab,i_specie,i_box) = &
                   State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)

!           Extra_mass_1D(i_box,i_slab,i_specie) = &
!                   State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
         enddo
         ENDDO



         ! uodate Total_extra(i_box)
         V_grid_1D       = box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp
         mass_plume      = V_grid_1D/2 *SUM(Cslab(1:n_slab_max))

         mass_plume_new = V_grid_1D &
              *SUM( box_concnt_1D(1:n_slab_max,1,i_box) )

         D_mass_plume = mass_plume_new - mass_plume

         IF(D_mass_plume>0)THEN
!           WRITE(6,*)'--- 2 to 1 grid cell:', i_box, mass_plume, mass_plume_new
           Total_extra(i_box) = Total_extra(i_box) + D_mass_plume
         ELSE 
           ! WRITE(6,*)'--- ERROR in 2 to 1 grid cell in 1D grid ---', i_box
         ENDIF


       ENDIF ! IF( 2*eddy_B*Dt2/(box_Rb(i_box)**2)>1 )

    
       IF(2*eddy_B*Dt2/(box_Rb(i_box)**2)>0.5 )THEN
         GOTO 300
       ENDIF


       !--------------------------------------------------------------------
       ! Begin the loop to calculate the advection & diffusion in slab model
       !--------------------------------------------------------------------
!       DO i_specie = 1,N_specie
       DO i_specie = 1, 1
       do t1s=1,NINT(Dt/Dt2)


       !--------------------------------------------------------------------
       ! update the boundary grids in 2D model based on background 
       ! contentration
       !--------------------------------------------------------------------
! ???      box_concnt_1D(i_box,1,1)          = backgrd_concnt(1)
! ???      box_concnt_1D(i_box,n_slab_max,1) = backgrd_concnt(1)

!       Extra_mass_1D(i_box,1,1)          = backgrd_concnt(1)
!       Extra_mass_1D(i_box,n_slab_max,1) = backgrd_concnt(1)


       V_grid_1D       = box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp
       C1d_prev       = box_concnt_1D(:,1,i_box)
!       C1d_prev_extra = Extra_mass_1D(i_box,:,1)

       !--------------------------------------------------------------------
       ! For deformation of cross-section caused by wind shear 
       ! (A.D.Naiman et al., 2010):
       ! This should be moved to the outer loop !!! shw
       !--------------------------------------------------------------------

       theta_previous   = box_theta(i_box)
       box_theta(i_box) = ATAN( TAN(box_theta(i_box)) + wind_s_shear*Dt2 )

       ! make sure use box_theta or TAN(box_theta)  ??? 
       box_Ra(i_box) = box_Ra(i_box) &
                                  * (TAN(box_theta(i_box))**2+1)**0.5 &
                                  / (TAN(theta_previous)**2+1)**0.5
       box_Rb(i_box) = box_Rb(i_box) &
                                  * (TAN(box_theta(i_box))**2+1)**(-0.5) &
                                  / (TAN(theta_previous)**2+1)**(-0.5)

       !--------------------------------------------------------------------
       ! For the concentration change caused by eddy diffusion:
       ! box_concnt(n_box,N_slab,N_specie)
       !--------------------------------------------------------------------

       ! ignore the diffusivity in long radius direction
       eddy_B = eddy_v*SIN(abs(box_theta(i_box))) &
                          + eddy_h*COS(box_theta(i_box)) ! b



       !---------------------------------------------------------------------
       ! Calculate the diffusion in slab grids
       !---------------------------------------------------------------------

       Cslab       = box_concnt_1D(:,N_specie,i_box)
!       Extra_Cslab = Extra_mass_1D(i_box,:,N_specie)

       Concnt1D_bdy(1)          = backgrd_concnt(1)
       Concnt1D_bdy(n_slab_max2) = backgrd_concnt(1)
       Concnt1D_bdy(2:n_slab_max2-1) = box_concnt_1D(1:n_slab_max,1,i_box)




       IF(2.0*eddy_B*Dt2/(box_Rb(i_box)**2)>1) WRITE(6,*)'*** CFL ERROR ***'

!       box_concnt_1D(i_box,2:n_slab_max-1,1) = Cslab(2:n_slab_max-1)     &
!              + Dt2*eddy_B *(Cslab(3:n_slab_max)-2*Cslab(2:n_slab_max-1) &
!                            +Cslab(1:n_slab_max-2)) /box_Rb(i_box)**2

!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '2-1D', i_box,&
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,        &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp *       &
!                                        SUM(box_concnt_1D(i_box,:,1)),   &
!        Total_extra(i_box), backgrd_concnt(1)


!       IF(i_box==222.and.t1s==1) WRITE(6,*) '1-1D', &
!         box_concnt_1D(i_box,n_slab_50:n_slab_50+3,1)
!       IF(i_box==222.and.t1s==1) WRITE(6,*) '1-1D', &
!         box_concnt_1D(i_box,n_slab_max-3:n_slab_max,1),Concnt1D_bdy(n_slab_max2)
!       IF(i_box==222.and.t1s==1) WRITE(6,*) Concnt1D_bdy(n_slab_max2), &
!                                          V_grid_1D*SUM(box_concnt_1D(i_box,:,1))


       box_concnt_1D(1:n_slab_max,1,i_box) = Concnt1D_bdy(2:n_slab_max2-1)     &
        + Dt2*eddy_B*(Concnt1D_bdy(3:n_slab_max2)-2*Concnt1D_bdy(2:n_slab_max2-1) &
                     +Concnt1D_bdy(1:n_slab_max2-2)) /box_Rb(i_box)**2


!       IF(i_box==222.and.t1s==1) WRITE(6,*) '1-1D', &
!         box_concnt_1D(i_box,n_slab_50:n_slab_50+3,1)
!       IF(i_box==222.and.t1s==1) WRITE(6,*) '1-1D', &
!         box_concnt_1D(i_box,n_slab_max-3:n_slab_max,1),Concnt1D_bdy(n_slab_max2)
!       IF(i_box==222.and.t1s==1) WRITE(6,*) Concnt1D_bdy(n_slab_max2), &
!                                          V_grid_1D*SUM(box_concnt_1D(i_box,:,1))

!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '3-1D', i_box,&
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,        &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp *       &
!                                        SUM(box_concnt_1D(i_box,:,1)),   &
!        Total_extra(i_box), backgrd_concnt(1)



!       Extra_mass_1D(i_box,2:n_slab_max-1,1) = Extra_Cslab(2:n_slab_max-1)      &
!         + Dt2*eddy_B *(Extra_Cslab(3:n_slab_max)-2*Extra_Cslab(2:n_slab_max-1) &
!                       +Extra_Cslab(1:n_slab_max-2)) /box_Rb(i_box)**2


       !====================================================================
       ! Update the concentration in the background grid cell
       ! after the interaction with 1D plume
       !====================================================================

       grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ![cm3]


       mass_plume = V_grid_1D * SUM( C1d_prev(1:n_slab_max) )
       mass_plume_new = V_grid_1D * SUM(box_concnt_1D(1:n_slab_max,1,i_box))

       ! [molec]
       D_mass_plume = mass_plume_new - mass_plume


       IF(D_mass_plume<0)THEN
       ! generally D_mass_plume is negative

       Total_extra(i_box) = Total_extra(i_box) &
              + Total_extra(i_box)/mass_plume *D_mass_plume

       D_mass_plume = (1-Total_extra(i_box)/mass_plume)*D_mass_plume
       ENDIF


       backgrd_concnt(1) = ( backgrd_concnt(1)*grid_volumn &
                                              - D_mass_plume) / grid_volumn

       State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = &
                                                backgrd_concnt(1)


! ???       box_concnt_1D(i_box,1,1)          = backgrd_concnt(1)
! ???       box_concnt_1D(i_box,n_slab_max,1) = backgrd_concnt(1)

!       IF(i_box==1)THEN
!         WRITE(6,*)'--- Check advection-diffusion in 1D:',            &
!                   SUM(box_concnt_1D(i_box,2:n_slab_max-1,1) &
!                     -Extra_mass_1D(i_box,2:n_slab_max-1,1)) &
!                   *box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,&
!                  D_mass_plume, D_mass_extra
!       ENDIF


!      IF(i_box==222.and.Judge_plume(i_box)==1) WRITE(6,*) '4-1D', i_box,&
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp,        &
!        box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp *       &
!                                        SUM(box_concnt_1D(i_box,:,1)),   &
!        Total_extra(i_box), backgrd_concnt(1)


       !===================================================================
       ! Second judge:
       ! According to a 2-order process rate, if there is no big difference
       ! between PIG and Eulerian model, then dissolve the PIG into
       ! background grid cell.                                          shw
       ! test this once big time step
       !===================================================================

       Rate_Mix = V_grid_1D * SUM( box_concnt_1D(1:n_slab_max,1,i_box)**2 ) &
                + backgrd_concnt(1)**2 * grid_volumn

       Eul_concnt = ( SUM(box_concnt_1D(1:n_slab_max,1,i_box)) *V_grid_1D &
                    + backgrd_concnt(1)*grid_volumn ) / grid_volumn
       Rate_Eul = Eul_concnt**2*grid_volumn




!       IF(i_box==500.and.t1s==1) WRITE(6,*) '*** Second judge:', &
!         Rate_Mix, Rate_Eul, backgrd_concnt(1)




       IF(ABS(Rate_Mix-Rate_Eul)/Rate_Eul<0.3)THEN !!! shw ???


!         backgrd_concnt(1) = &
!            ( SUM(  box_concnt_1D(i_box, 2:n_slab_max-1,1)             &
!                  - Extra_mass_1D(i_box, 2:n_slab_max-1, 1)  )         &
!              * V_grid_1D &
!                + backgrd_concnt(1)*grid_volumn ) / grid_volumn

         backgrd_concnt(1) = &
           ( SUM(box_concnt_1D(1:n_slab_max,1,i_box)) * V_grid_1D   &
                - Total_extra(i_box)  &
                + backgrd_concnt(1)*grid_volumn ) / grid_volumn

         State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) = backgrd_concnt(1)

         Judge_plume(i_box) = 0


         GOTO 200 ! skip this box, go to next box

       ENDIF ! IF(ABS(Rate_Mix-Rate_Eul)/Rate_Eul<0.01)THEN


       !===================================================================
       ! Third judge:
       ! If this is the final time step, the model is going to end
       ! dissolve all the plume into GCM in the last time step
       !===================================================================
!       IF ( ITS_TIME_FOR_EXIT() ) THEN
        IF(DAY==31 .and. HOUR==23)THEN

!         backgrd_concnt(1) = &
!            ( SUM(  box_concnt_1D(i_box, 2:n_slab_max-1,1)             &
!                  - Extra_mass_1D(i_box, 2:n_slab_max-1, 1)  )         &
!              * V_grid_1D &
!                + backgrd_concnt(1)*grid_volumn ) / grid_volumn

         backgrd_concnt(1) = &
           ( SUM(box_concnt_1D(1:n_slab_max,1,i_box)) * V_grid_1D   &
                - Total_extra(i_box)  &
                + backgrd_concnt(1)*grid_volumn ) / grid_volumn


         State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) =    &
                                                     backgrd_concnt(1)

         Judge_plume(i_box) = 0


         GOTO 200 ! skip this box, go to next box

       ENDIF


       !===================================================================
       ! Fourth judge:
       ! If the plume touch the tropopause, dissolve the plume
       !===================================================================
       IF ( box_lev(i_box)>State_Met%TROPP(i_lon,i_lat) ) THEN

!         backgrd_concnt(1) = &
!            ( SUM(  box_concnt_1D(i_box, 2:n_slab_max-1,1)             &
!                  - Extra_mass_1D(i_box, 2:n_slab_max-1, 1)  )         &
!              * V_grid_1D   &
!                + backgrd_concnt(1)*grid_volumn ) / grid_volumn

         backgrd_concnt(1) = &
           ( SUM(box_concnt_1D(1:n_slab_max,1,i_box)) * V_grid_1D   &
                - Total_extra(i_box)  &
                + backgrd_concnt(1)*grid_volumn ) / grid_volumn


         State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1) =    &
                                                     backgrd_concnt(1)

         Judge_plume(i_box) = 0

         WRITE(6,*) '*** Dissolved in 4th Judge ***', i_box

         GOTO 200 ! skip this box, go to next box

       ENDIF


       enddo ! do t1s=1,NINT(Dt/Dt2)
       ENDDO ! do i_Species = 1,nSpecies



       !===================================================================
       ! For filament structure:
       ! If slab length(Ra) is bigger than 2* horizontal resolution (2*Dx),
       ! split slab into five smaller segment (Ra/5)
       !===================================================================
       N_split = 3

!       IF(i_box==4) WRITE(6,*)'--- n_box, Ra ---', n_box_max, box_Ra(i_box), Dx

!!! shw
       IF( box_Ra(i_box) > Dx*110.0*1000.0 ) THEN


         ! extend the total number of box to include the new box
         n_box_prev = n_box
         n_box = n_box + (N_split-1)

         allocate(Data1D_Int(n_box_prev))
         allocate(Data1D(n_box_prev))
         allocate(Data2D(N_specie, n_box_prev))
         allocate(Data3D(n_slab_max, N_specie, n_box_prev))
         allocate(Data4D(n_x_max, n_y_max, N_specie, n_box_prev))

         Data1D = box_lon(:)
         deallocate(box_lon)
         allocate(box_lon(n_box))
         box_lon(1:n_box_prev) = Data1D

         Data1D = box_lat(:)
         deallocate(box_lat)
         allocate(box_lat(n_box))  
         box_lat(1:n_box_prev) = Data1D
          
         Data1D = box_lev(:)
         deallocate(box_lev)
         allocate(box_lev(n_box))
         box_lev(1:n_box_prev) = Data1D

!         Data1D = box_lon_edge(:)
!         deallocate(box_lon_edge)
!         allocate(box_lon_edge(n_box))
!         box_lon_edge(1:n_box_prev) = Data1D

!         Data1D = box_lat_edge(:)
!         deallocate(box_lat_edge)
!         allocate(box_lat_edge(n_box))
!         box_lat_edge(1:n_box_prev) = Data1D

!         Data1D = box_lev_edge(:)
!         deallocate(box_lev_edge)
!         allocate(box_lev_edge(n_box))
!         box_lev_edge(1:n_box_prev) = Data1D


         Data1D = Pdx(:)
         deallocate(Pdx)
         allocate(Pdx(n_box))
         Pdx(1:n_box_prev) = Data1D

         Data1D = Pdy(:)
         deallocate(Pdy)
         allocate(Pdy(n_box))
         Pdy(1:n_box_prev) = Data1D


         Data1D_Int = Judge_plume(:)
         deallocate(Judge_plume)
         allocate(Judge_plume(n_box))
         Judge_plume(1:n_box_prev) = Data1D_Int


         Data4D = box_concnt_2D(:,:,:,:)
         deallocate(box_concnt_2D)
         allocate(box_concnt_2D(n_x_max,n_y_max,N_specie,n_box))
         box_concnt_2D(:,:,:,1:n_box_prev) = Data4D

!         Data4D = Extra_mass_2D(:,:,:,:)
!         deallocate(Extra_mass_2D)
!         allocate(Extra_mass_2D(n_box,n_x_max,n_y_max,N_specie))
!         Extra_mass_2D(1:n_box_prev,:,:,:) = Data4D


         Data3D = box_concnt_1D(:,:,:)
         deallocate(box_concnt_1D)
         allocate(box_concnt_1D(n_slab_max,N_specie,n_box))
         box_concnt_1D(:,:,1:n_box_prev) = Data3D

!         Data3D = Extra_mass_1D(:,:,:)
!         deallocate(Extra_mass_1D)
!         allocate(Extra_mass_1D(n_box,n_slab_max,N_specie))
!         Extra_mass_1D(1:n_box_prev,:,:) = Data3D


         Data1D = box_Ra(:)
         deallocate(box_Ra)
         allocate(box_Ra(n_box))
         box_Ra(1:n_box_prev) = Data1D

         Data1D = box_Rb(:)
         deallocate(box_Rb)
         allocate(box_Rb(n_box))
         box_Rb(1:n_box_prev) = Data1D

         Data1D = box_theta(:)
         deallocate(box_theta)
         allocate(box_theta(n_box))
         box_theta(1:n_box_prev) = Data1D

         Data1D = box_length(:)
         deallocate(box_length)
         allocate(box_length(n_box))
         box_length(1:n_box_prev) = Data1D


!         Data1D = box_u(:)
!         deallocate(box_u)
!         allocate(box_u(n_box))
!         box_u(1:n_box_prev) = Data1D
!
!         Data1D = box_v(:)
!         deallocate(box_v)
!         allocate(box_v(n_box))
!         box_v(1:n_box_prev) = Data1D
!
!         Data1D = box_omeg(:)
!         deallocate(box_omeg)
!         allocate(box_omeg(n_box))
!         box_omeg(1:n_box_prev) = Data1D

         Data1D = box_alpha(:)
         deallocate(box_alpha)
         allocate(box_alpha(n_box))
         box_alpha(1:n_box_prev) = Data1D

!         Data1D = box_Ptemp(:)
!         deallocate(box_Ptemp)
!         allocate(box_Ptemp(n_box))
!         box_Ptemp(1:n_box_prev) = Data1D

         Data1D = Total_extra(:)
         deallocate(Total_extra)
         allocate(Total_extra(n_box))
         Total_extra(1:n_box_prev) = Data1D


!         Data1D_Int = Plume_I(:)
!         deallocate(Plume_I)
!         allocate(Plume_I(n_box))
!         Plume_I(1:n_box_prev) = Data1D_Int

!         Data1D_Int = Plume_J(:)
!         deallocate(Plume_J)
!         allocate(Plume_J(n_box))
!         Plume_J(1:n_box_prev) = Data1D_Int

!         Data1D_Int = Plume_L(:)
!         deallocate(Plume_L)
!         allocate(Plume_L(n_box))
!         Plume_L(1:n_box_prev) = Data1D_Int

         Data1D = Lifetime(:)
         deallocate(Lifetime)
         allocate(Lifetime(n_box))
         Lifetime(1:n_box_prev) = Data1D


         deallocate(Data1D_Int)
         deallocate(Data1D)
         deallocate(Data2D)
         deallocate(Data3D)
         deallocate(Data4D)

         !-----------------------------------------------------------------------
         ! set initial location for new added boxes from splitting
         !-----------------------------------------------------------------------

         ! Based on Great_circle distance
!         DlonN = 2.0 * 180/PI *ASIN( ABS( SIN(box_Ra(i_box)*COS(box_alpha)/N_split /Re *0.5) &
!                                     / COS(box_lat(i_box)/180*PI)**2 ) )
        
         DlonN = ABS(box_Ra(i_box)*COS(box_alpha(i_box))/N_split) / &
                        (2*PI *(Re*COS(box_lat(i_box)/180*PI)) ) * 360

         DlatN = ABS(box_Ra(i_box)*SIN(box_alpha(i_box))/N_split) / (2*PI*Re) * 360


         DO i = 1, (N_split-1)/2, 1
           box_lon(n_box_prev+i) = box_lon(i_box) - i*DlonN
!           box_lon_edge(n_box_prev+i) = box_lon_edge(i_box) - i*DlonN

           box_lat(n_box_prev+i) = box_lat(i_box) - i*DlatN
!           box_lat_edge(n_box_prev+i) = box_lat_edge(i_box) - i*DlatN
         ENDDO

         DO i = 1, (N_split-1)/2, 1
           box_lon(n_box_prev+(N_split-1)/2+i) = box_lon(i_box) + i*DlonN
!           box_lon_edge(n_box_prev+(N_split-1)/2+i) = box_lon_edge(i_box) + i*DlonN

           box_lat(n_box_prev+(N_split-1)/2+i) = box_lat(i_box) + i*DlatN
!           box_lat_edge(n_box_prev+(N_split-1)/2+i) = box_lat_edge(i_box) + i*DlatN
         ENDDO 

!         WRITE(6,*) '--- box_Ra(i_box), box_alpha, DlonN, DlatN, 2.0*Dx:'
!         WRITE(6,*) box_Ra(i_box), box_alpha,  DlonN, DlatN, 2.0*Dx



         do ii_box = n_box_prev+1, n_box, 1

           box_lev(ii_box)     = box_lev(i_box)
!           box_lev_edge(ii_box)     = box_lev_edge(i_box)

           ! not need to set Pdx and Pdy

           box_Ra(ii_box)      = box_Ra(i_box)/N_split
           box_Rb(ii_box)      = box_Rb(i_box)

           box_theta(ii_box)   = box_theta(i_box)
           Lifetime(ii_box)     = Lifetime(i_box)
           Judge_plume(ii_box) = 1

           box_length(ii_box)  = box_length(i_box)

!           box_u(ii_box)       = box_u(i_box)
!           box_v(ii_box)       = box_v(i_box)
!           box_omeg(ii_box)    = box_omeg(i_box)

           box_alpha(ii_box)   = box_alpha(i_box)

!           box_Ptemp(ii_box)   = box_Ptemp(i_box)


           box_concnt_1D(:,:,ii_box) = box_concnt_1D(:,:,i_box)
!           Extra_mass_1D(ii_box,:,:) = Extra_mass_1D(i_box,:,:)
           Total_extra(ii_box) = Total_extra(i_box)/N_split

           ! No use
           Pdx(ii_box)       = 0.0
           Pdy(ii_box)       = 0.0

         enddo

         box_Ra(i_box) = box_Ra(i_box)/N_split
         Total_extra(i_box) = Total_extra(i_box)/N_split

       ENDIF ! IF( box_Ra(i_box) > 2*Dx )


!       call cpu_time(finish)
!       WRITE(6,*)'Time (finish-start) for 1D:', i_box, finish-start


       ENDIF ! IF(Judge_plume(i_box)==1) THEN


 200 CONTINUE


       V_grid_2D       = Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp

!       IF(i_box==2095) WRITE(6,*) 6, i_box, &
!         V_grid_2D*SUM(box_concnt_2D(i_box,2:n_x_max-1,2:n_y_max-1,1)), &
!         V_grid_2D*SUM(Extra_mass_2D(i_box,2:n_x_max-1,2:n_y_max-1,1))


    enddo ! do i_box = 1,n_box_max

    !=======================================================================
    ! Convert specie back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'plume_mod.F90' )
       RETURN
    ENDIF







! IF plume is dissolved, delete the plume to release memor
! Judge_plume(i_box) = 0

!     n_box_prev = n_box
!     n_box = n_box - N_plume0
!
!     allocate(Data1D_Int(n_box_prev))
!     allocate(Data1D(n_box_prev))
!     allocate(Data2D(N_specie, n_box_prev))
!     allocate(Data3D(n_slab_max, N_specie, n_box_prev))
!     allocate(Data4D(n_x_max, n_y_max, N_specie, n_box_prev))










    ! N_curr is used to add particles, here add 5 parcles every time step


    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)
    nullify(omeg)

    nullify(Ptemp)
    nullify(P_BXHEIGHT)

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

  real(fp) function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2      ! unit is degree
    real(fp)     :: xx1, yy1, xx2, yy2  ! unit is radian
    !real(fp) :: PI, Re

    xx1 = x1/180.0*PI
    yy1 = y1/180.0*PI
    xx2 = x2/180.0*PI
    yy2 = y2/180.0*PI

!     Distance_Circle = Re * 
!       ATAN( SQRT( (COS(y2)*SIN(x2-x1))**2 + &
!              ( COS(y1)*SIN(y2)-SIN(y1)*OCS(y2)COS((x2-x1)) )**2 ) &
!            / (SIN(y1)*SIN(y2)+COS(y1)*COS(y2)*COS(x2-x1)) )

    ! output distance is in unit of [m]
    Distance_Circle = Re * 2.0 * ASIN(SQRT( (SIN((yy1-yy2)*0.5))**2.0 &
                     +COS(yy1)*COS(yy2)*(SIN((xx1-xx2)*0.5))**2.0 ))
    return
  end function


  !-------------------------------------------------------------------  
  ! calculate the wind_s (inside a plume sross-section) shear along pressure
  ! direction
  !-------------------------------------------------------------------  

  real(fp) function Wind_shear_s(u, v, P_BXHEIGHT, plume_alpha, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: plume_alpha
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

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2

      kk = k + init_lev - 1

      if(init_lon==0)then
          u_lonlat(k) =  Weight(1,1) * u(IIPAR,init_lat,kk) &
                       + Weight(1,2) * u(IIPAR,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(IIPAR,init_lat,kk) &
                       + Weight(1,2) * v(IIPAR,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      else if(init_lon==IIPAR)then
          u_lonlat(k) =  Weight(1,1) * u(init_lon,init_lat,kk) &
                       + Weight(1,2) * u(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * u(1,init_lat,kk)   &
                       + Weight(2,2) * u(1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(init_lon,init_lat,kk) &
                       + Weight(1,2) * v(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * v(1,init_lat,kk)   &
                       + Weight(2,2) * v(1,init_lat+1,kk)
      else
          u_lonlat(k) =  Weight(1,1) * u(init_lon,init_lat,kk) &
                       + Weight(1,2) * u(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(init_lon,init_lat,kk) &
                       + Weight(1,2) * v(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      endif

      wind_s(k) = u_lonlat(k)*COS(plume_alpha) + v_lonlat(k)*SIN(plume_alpha)

    enddo


    ! second vertical shear of wind_s

    ! This code should be changed !!!
   ! Because it is the pressure center in [hPa] instead of height center in
    ! [m]
    ! Delt_height    = 0.5 * ( P_BXHEIGHT(init_lon,init_lat,init_lev) +
    ! P_BXHEIGHT(init_lon,init_lat,init_lev+1) )
    if(init_lon==0)then
      Delt_height = Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev),    &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                 + Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev+1),   &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else if(init_lon==IIPAR)then
      Delt_height = Pa2meter( P_BXHEIGHT(1,init_lat,init_lev),        &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(1,init_lat,init_lev+1),       &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else
      Delt_height = Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),&
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    endif


    ! find the z height of each pressure level in GEOS-Chem

    Wind_shear_s = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function

!-------------------------------------------------------------------
! transform the pressure level [Pa] to the height level [m]
!-------------------------------------------------------------------

  real(fp) function Pa2meter(Box_height, P1, P2, Judge)
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

  real(fp) function Vertical_shear(var, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
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

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2
      kk            = k + init_lev - 1
      
      IF(init_lon==0)THEN
        var_lonlat(k) =  Weight(1,1) *var(IIPAR,i_lat,kk)   &
                       + Weight(1,2) *var(IIPAR,i_lat+1,kk)   &
                       + Weight(2,1) *var(1,i_lat,kk) &
                       + Weight(2,2) *var(1,i_lat+1,kk)
      ELSE IF(init_lon==IIPAR)THEN
        var_lonlat(k) =  Weight(1,1) *var(IIPAR,i_lat,kk)   &
                       + Weight(1,2) *var(IIPAR,i_lat+1,kk)   &
                       + Weight(2,1) *var(1,i_lat,kk) &
                       + Weight(2,2) *var(1,i_lat+1,kk)
      ELSE
        var_lonlat(k) =  Weight(1,1) *var(init_lon,i_lat,kk)   &
                       + Weight(1,2) *var(init_lon,i_lat+1,kk)   &
                       + Weight(2,1) *var(init_lon+1,i_lat,kk) &
                       + Weight(2,2) *var(init_lon+1,i_lat+1,kk)
      ENDIF
    enddo


    ! second vertical shear of wind
    if(init_lon==0)then
      Delt_height = Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev),    &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev+1),   &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else if(init_lon==IIPAR)then
      Delt_height = Pa2meter( P_BXHEIGHT(1,init_lat,init_lev),        &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(1,init_lat,init_lev+1),       &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else
      Delt_height = Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),&
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    endif
    ! find the z height of each pressure level in GEOS-Chem

    Vertical_shear = ( var_lonlat(2) - var_lonlat(1) ) / Delt_height

    return
  end function


!===================================================================
!
!===================================================================
  REAL(fp) FUNCTION Get_XYscale(concnt1_2D, i_box, frac, axis)

    IMPLICIT NONE

    REAL(fp) :: concnt1_2D(n_x_max, n_y_max), frac
    INTEGER  :: i_box, axis

    INTEGER  :: i, j
    REAL(fp) :: temp, C_sum, mass_total

    REAL(fp) :: D_len
    INTEGER  :: N_max, N_frac

    REAL(fp), allocatable :: concnt1_2D_sum(:)


      IF(axis==2)THEN ! sum y
        N_max = n_x_max
        D_len = Pdx(i_box)
      ELSE IF(axis==1)THEN ! sum x
        N_max = n_y_max
        D_len = Pdy(i_box)
      ENDIF

    allocate(concnt1_2D_sum(N_max))

      concnt1_2D_sum = SUM(concnt1_2D, DIM=axis)


      ! sort the array from high value to low:
      DO i = N_max-1, 1, -1
      DO j = 1, i
        IF(concnt1_2D_sum(j)<concnt1_2D_sum(j+1))THEN
          temp = concnt1_2D_sum(j)
          concnt1_2D_sum(j) = concnt1_2D_sum(j+1)
          concnt1_2D_sum(j+1) = temp
        ENDIF
      ENDDO
      ENDDO


      mass_total = SUM(concnt1_2D_sum)
      C_sum = 0.0

      ! find the XYscale containing frac of total mass
      DO i = 1, N_max, 1
        C_sum = C_sum + concnt1_2D_sum(i)
        if(C_sum>mass_total*frac)then
          N_frac = i
          EXIT
        endif
      ENDDO


      Get_XYscale = N_frac * D_len

    RETURN

  END FUNCTION



!  REAL(fp) FUNCTION Find_theta(theta1, concnt1_2D, Height, Ibox)
!    ! n_x_max, n_y_max, Pdx, Pdy are global variables
!
!    IMPLICIT NONE
!
!    REAL(fp)    :: theta1, concnt1_2D(n_x_max, n_y_max), Height    
!    REAL(fp)    :: R
!    REAL(fp)    :: test(100), X(100), Y(100), Conc(100)
!
!    REAL(fp)    :: angle(1)
!
!    INTEGER     :: i, Ibox
!
!      R = Height/COS(theta1) ! 2.5 is a random number
!
!
!      DO i= 1, 100
!        test(i) = PI/2 - (PI/2-theta1)*2/100*i
!        Y(i)    = R *COS( test(i) )
!        X(i)    = R *SIN( test(i) )
!        Conc(i) = Interplt_2D(X(i), Y(i), concnt1_2D, Ibox, 1)
!      ENDDO
!
!
!      angle      = test( MAXLOC(Conc(:)) )
!      Find_theta = angle(1)
!
!    return 
!
!  END FUNCTION



  SUBROUTINE Slab_init(theta1, bkgrd_concnt, Ibox, Pc_2D, Height1)
 
    IMPLICIT NONE

    REAL(fp)    :: theta1, bkgrd_concnt, Height1
    REAL(fp)    :: Pc_2D(n_x_max,n_y_max) !, Ec_2D(n_x_max,n_y_max)
    INTEGER     :: Ibox


    INTEGER, parameter     :: Nb = n_slab_max
    INTEGER, parameter     :: Na = 128

    INTEGER     :: Nb_mid, Na_mid

    REAL(fp)    :: X2d(Na,Nb), Y2d(Na,Nb), C2d(Na,Nb) !, Extra_C2d(Na,Nb)

    REAL(fp)    :: LenB, LenA
    REAL(fp)    :: Adx, Ady, Bdx, Bdy
    REAL(fp)    :: Prod, M, Lb, La

    REAL(fp)    :: C_slab(Nb) !, Extra_slab(Nb)

    INTEGER     :: i, j

      Nb_mid = INT(Nb/2)
      Na_mid = INT(Na/2)


      LenB = Pdy(Ibox) ! 8.0 *Height1 / Nb /SIN(theta1)
      LenA = 1.3 *Height1 *TAN(theta1) /Na /SIN(theta1)


      ! interval in long radius
      Adx = LenA*SIN(theta1)
      Ady = LenA*COS(theta1)

      ! interval in short radius
      Bdy = LenB*SIN(theta1)
      Bdx = LenB*COS(theta1)


      ! find the location of 1D grid in 2D XY grids
      DO i=1,Na,1
        X2d(i,Nb_mid) = -Adx*Na_mid + Adx*(i-0.5)
        Y2d(i,Nb_mid) = -Ady*Na_mid + Ady*(i-0.5)
      ENDDO


      X2d(:,Nb_mid) = X2d(:,Nb_mid) + 0.5*Bdx
      Y2d(:,Nb_mid) = Y2d(:,Nb_mid) - 0.5*Bdy

      DO j=Nb_mid+1, Nb, 1
        X2d(:,j) = X2d(:,j-1) - Bdx
        Y2d(:,j) = Y2d(:,j-1) + Bdy
      ENDDO


      DO j=Nb_mid-1, 1, -1
        X2d(:,j) = X2d(:,j+1) + Bdx
        Y2d(:,j) = Y2d(:,j+1) - Bdy
      ENDDO


      DO i=1,Na,1
      DO j=1,Nb,1
        C2d(i,j)       = Interplt_2D(X2d(i,j), Y2d(i,j), Pc_2D, bkgrd_concnt, Ibox, 2)
!        Extra_C2d(i,j) = Interplt_2D(X2d(i,j), Y2d(i,j), Ec_2D, Ibox, 2)
      ENDDO
      ENDDO


      ! define the length/width of slab based on a seconde order chamical reaction

!      Prod = SUM( C2d(Na_mid,:)**2 ) *(LenA*LenB)
!      M    = SUM( C2d(Na_mid,:) ) *(LenA*LenB)
!      Lb   = LenB
!      La   = M**2 /Prod / Lb

!      DO i=1,Nb,1
!        C_slab(i)     = SUM(C2d(:,i)) *(LenA*LenB)/ (La*Lb)
!        Extra_slab(i) = SUM(Extra_C2d(:,i)) *(LenA*LenB)/ (La*Lb)
!      ENDDO

        ! define the length/width of slab based on 95% total mass
        Lb = LenB
        La = Height1 *TAN(theta1)

        DO i=1,Nb,1
          C_slab(i)     = SUM(C2d(:,i)) *(LenA*LenB)/ (La*Lb)
!          Extra_slab(i) = SUM(Extra_C2d(:,i)) *(LenA*LenB)/ (La*Lb)
        ENDDO


      box_Ra(Ibox) = La
      box_Rb(Ibox) = Lb
      box_concnt_1D(:,1,Ibox) = C_slab(:)
!      Extra_mass_1D(Ibox,:,1) = Extra_slab(:)

  END SUBROUTINE



  REAL(fp) FUNCTION Interplt_2D(x0, y0, C_2D, bkgrd_concnt, Ibox, ids)
    ! n_x_max, n_y_max, Pdx, Pdy are global variables

    IMPLICIT NONE

    INTEGER     :: Ibox

    REAL(fp)    :: x0, y0, C_2D(n_x_max,n_y_max), bkgrd_concnt
    REAL(fp)    :: X1d(n_x_max), Y1d(n_y_max)
    REAL(fp)    :: C1, C2

    integer     :: ids ! 1 for Find_theta(); 2 for Slab_init()
    integer     :: Ix0, Iy0
    integer     :: i, j

      ! define the coordinate system
      DO i=1, n_x_max
        X1d(i) = Pdx(Ibox)*(i-n_x_mid)
      ENDDO
      DO j=1, n_y_max
        Y1d(j) = Pdy(Ibox)*(j-n_y_mid)
      ENDDO


      Ix0 = floor( (x0-X1d(1)) / Pdx(Ibox) ) + 1
      Iy0 = floor( (y0-Y1d(1)) / Pdy(Ibox) ) + 1

      IF(Ix0<1 .or. Ix0>=n_x_max)THEN
        Interplt_2D = bkgrd_concnt
!        IF(ids==1)WRITE(6,*)"*** ERROR: Find_theta() out bounds ***"

      ELSE IF((Iy0<1 .or. Iy0>=n_y_max))THEN
        Interplt_2D = bkgrd_concnt
!        IF(ids==1)WRITE(6,*)"*** ERROR: Find_theta() out bounds***"

      ELSE
        C1 = Interplt_linear(x0, X1d(Ix0), X1d(Ix0+1), &
                C_2D(Ix0,Iy0),C_2D(Ix0+1,Iy0))

        C2 = Interplt_linear(x0, X1d(Ix0), X1d(Ix0+1), &
                C_2D(Ix0,Iy0+1), C_2D(Ix0+1,Iy0+1))

        Interplt_2D = Interplt_linear(y0, Y1d(Iy0), &
                                          Y1d(Iy0+1), C1, C2)

      ENDIF

    return

  END FUNCTION


  REAL(fp) FUNCTION Interplt_linear(xx, x1, x2, Cx1, Cx2)

    IMPLICIT NONE

    REAL(fp)    :: xx, x1, x2, Cx1, Cx2


    Interplt_linear = ( Cx1*(x2-xx) + Cx2*(xx-x1) ) /(x2-x1)

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
    INTEGER                   :: i_box, i_slab
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

    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)


    FILENAME   = 'Lagrange_xyz_' // TRIM(ADJUSTL(YEAR_C)) // '-'   &
        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'


    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )


    DO i_box = 1,n_box,1
      WRITE(261, *) i_box, Lifetime(i_box)/3600/24, Judge_plume(i_box)
    ENDDO ! DO i_box = 1,n_box,1

    ENDIF ! IF(mod(tt,1440)==0)THEN


!-------------------------------------------------------------------
! Output plume location
!-------------------------------------------------------------------
!    IF(mod(tt,6)==0)THEN     ! output once every hour
    IF(mod(tt,14400000)==0)THEN   ! output once every day (24 hours)

      FILENAME2   = 'Plume_concentration_molec_' // TRIM(ADJUSTL(YEAR_C)) // &
         '-' //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) //      &
         '-' //TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) //    &
         ':' //TRIM(ADJUSTL(SECOND_C)) // '.txt'

      OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
            FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      WRITE(262,*)'total plume number:', n_box


      DO i_box = 1,n_box,1
!      IF(i_box<3)THEN !!! shw

        IF(Judge_plume(i_box)==0) THEN
          WRITE(262,*) i_box, ': dissolved'
        ENDIF

        IF(Judge_plume(i_box)==2) THEN
          WRITE(262,*) i_box, ': 2D model, total mass [molec]:'
          WRITE(262,*) Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
                * SUM(box_concnt_2D(:,:,N_specie,i_box))
        ENDIF

        IF(Judge_plume(i_box)==1) THEN
          WRITE(262,*) i_box, ': 1D model'
          WRITE(262,*) SUM(box_concnt_1D(1:n_slab_max,N_specie,i_box)) &
                   *box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp
        ENDIF

        WRITE(262,*)' '

!     ENDIF ! IF(i_box<3)THEN
     ENDDO ! DO i_box = 1,n_box,1


!       DO i_box = 1, n_box
!          WRITE(262,*) i_box, SUM(box_concnt(i_box,:,N_specie)*V_slab(i_box,:))
!       ENDDO

!        WRITE(262,*) box_concnt(1,:,N_specie)

    ENDIF ! IF(mod(tt,144)==0)THEN
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
    if (allocated(box_Ra))       deallocate(box_Ra)
    if (allocated(box_Rb))       deallocate(box_Rb)
    if (allocated(box_theta))    deallocate(box_theta)

    WRITE(6,'(a)') '--> Lagrange and Plume Module Cleanup <--'


  end subroutine lagrange_cleanup


END MODULE Lagrange_Mod
