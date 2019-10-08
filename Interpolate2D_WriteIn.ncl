;---------------------------------------------
; regrid_1.ncl
;
; Concepts illustrated:
;   - Interpolating from one grid to another using bilinear interpolation
;   - Drawing color-filled contours over a cylindrical equidistant map
;   - Paneling two plots vertically on a page
;   - Adding a common labelbar to paneled plots
;   - Adding units attributes to lat/lon arrays
;   - Turning off the individual contour labelbar
;   - Manually creating lat/lon coordinate arrays 
;---------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

month      = (/"01","02","03","04","05","06","07","08","09","10","11","12"/)
day 	   = (/"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"/)
Ndays      = (/31,28,31,30,31,30,31,31,30,31,30,31/)

do imon = 0,11,1
do iday = 0,Ndays(imon)-1,1
print((/imon/))
print((/iday/))

;---Read in netCDF file
  f_new = addfile("/n/home12/hongwei/hongwei/data/goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T3NVASM.5.12.4/2015/"+month(imon)+"/MERRA2_400.tavg3_3d_asm_Nv.2015"+month(imon)+day(iday)+".nc4","r")
  f_old = addfile("/n/home12/hongwei/hongwei/data/GEOS-Chem/MERRA2_A3dyn/2015/"+month(imon)+"/MERRA2.2015"+month(imon)+day(iday)+".A3dyn.4x5.nc4","w")

;---Read in data
  uwnd1 = f_new->U(:,::-1,:,:); remove cyclic point
  lon1  = f_new->lon(:)
  lat1  = f_new->lat(:)

  uwnd2 = f_old->U(:,:,:,:); remove cyclic point
  lon2  = f_old->lon(:)
  lat2  = f_old->lat(:)

;---------------------------------------------
; interpolate to new grid
;---------------------------------------------

  uwnd1_new = linint2(lon1,lat1,uwnd1,True,lon2,lat2,0)
  copy_VarMeta(uwnd2,uwnd1_new)
;  printVarSummary(uwnd1_new&lev)

; replace the U, V in old file with interpolated U, V from new file

  f_old->U = uwnd2               ; overwrite

end do
end do

end
