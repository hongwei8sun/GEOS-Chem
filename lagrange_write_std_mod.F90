!EOC
!------------------------------------------------------------------------------
!                  lagrange module output                                     !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lagrange_write_std_mod.F90
!
!
MODULE LAGRANGE_WRITE_STD_MOD

  USE precision_mod
  USE CMN_SIZE_mod
  USE PhysConstants   ! Physical constants: Re, PI, PI_180

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: LAGRANGE_WRITE_STD


  ! Fill value used in HEMCO  netCDF files.
  REAL(fp), PARAMETER :: FillValue = 1.e-31_fp

CONTAINS


!EOC
!------------------------------------------------------------------------------
!                  lagrange module output                                     !
!------------------------------------------------------------------------------
!BOP


  SUBROUTINE LAGRANGE_write_std( am_I_Root )
!
! !USES:
!
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
    USE PhysConstants, ONLY : PI, Re

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE, XMID, YMID       
    ! choose XEDGE or XMID ??        
    USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON

    ! Parameters for netCDF routines
    include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    LOGICAL,                    INTENT(IN   ) :: am_I_Root   ! root CPU?

!    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
!
! !LOCAL VARIABLES:
!    INTEGER                   :: I
    INTEGER                   :: YYYY, MM, DD, h, m, s
    REAL(f8), POINTER         :: Arr1D(:)
    INTEGER,  POINTER         :: Int1D(:)
    REAL(f4), POINTER         :: Arr3D(:,:,:)
    REAL(f4), POINTER         :: Arr4D(:,:,:,:)
    CHARACTER(LEN=255)        :: NcFile
    CHARACTER(LEN=255)        :: Pfx, title, Reference, Contact 
    INTEGER                   :: fId, lonId, latId, levId, TimeId
    INTEGER                   :: VarCt
    INTEGER                   :: nLon, nLat, nLev, nLevTmp, nTime 
!    INTEGER                   :: L

    CHARACTER(LEN=255), PARAMETER :: LOC = 'LAGRANGE_WRITE_STD (lagrange_write_std_mod.F90)' 

    real(fp), pointer :: X_mid(:)
    real(fp), pointer :: Y_mid(:)

    !=================================================================
    ! LAGRANGE_WRITE_STD begins here!
    !=================================================================

       X_mid => XMID(:,1,1)
       Y_mid => YMID(1,:,1)

       nLon     = IIPAR
       nLat     = JJPAR
       nLev     = LLPAR
       nTime    = 1       

       ! Create output file
       ! Pass CREATE_NC4 to make file format netCDF-4 (mps, 3/3/16)
       ! Now create netCDF file with time dimension as UNLIMITED (bmy, 3/8/17)
       CALL NC_Create( NcFile       = NcFile,                            &
                       Title        = "lagrange_test",                             &
                       nLon         = nLon,                              &
                       nLat         = nLat,                              &
                       nLev         = nLevTmp,                           &
                       nTime        = NF_UNLIMITED,                      &
                       fId          = fId,                               &
                       lonId        = lonId,                             &
                       latId        = latId,                             &
                       levId        = levId,                             &
                       timeId       = timeId,                            &
                       VarCt        = VarCt,                             &
                       CREATE_NC4   =.TRUE.                             )


    !-----------------------------------------------------------------
    ! Write grid dimensions (incl. time) 
    !-----------------------------------------------------------------

       ! Write longitude axis variable ("lon") to file   
       CALL NC_Var_Def( fId         = fId,                                &
                        lonId       = lonId,                              &
                        latId       = -1,                                 &
                        levId       = -1,                                 &
                        timeId      = -1,                                 &
                        VarName     = 'lon',                              &
                        VarLongName = 'Longitude',                        &
                        VarUnit     = 'degrees_east',                     &
                        Axis        = 'X',                                &
                        DataType    = f8,                                 &
                        VarCt       = VarCt,                              & 
                        Compress    = .TRUE.                             )
       CALL NC_Var_Write( fId, 'lon', Arr1D=X_MID )
       
       ! Write longitude axis variable ("lat") to file                                               
       CALL NC_Var_Def( fId         = fId,                                &                          
                        lonId       = latId,                              &                          
                        latId       = -1,                                 &                          
                        levId       = -1,                                 &                          
                        timeId      = -1,                                 &                          
                        VarName     = 'lat',                              &                          
                        VarLongName = 'Latitude',                        &                          
                        VarUnit     = 'degrees_east',                     &                          
                        Axis        = 'Y',                                &                          
                        DataType    = f8,                                 &                          
                        VarCt       = VarCt,                              &                          
                        Compress    = .TRUE.                             )                           
       CALL NC_Var_Write( fId, 'lat', Arr1D=Y_MID )
       
          ! Deallocate arrays
   
 !-----------------------------------------------------------------
    ! Cleanup
    !-----------------------------------------------------------------

    ! Close file
    CALL NC_CLOSE ( fId )


    ! Return 

  END SUBROUTINE LAGRANGE_write_std


END MODULE LAGRANGE_WRITE_STD_MOD

