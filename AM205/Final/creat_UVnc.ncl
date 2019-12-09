 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

 f0 = addfile("./OLD/MERRA2.20141231.A3dyn.2x25.nc4", "r")

 lat   = f0->lat
 lon   = f0->lon
 time0 = f0->time

 uwnd = f0->U(:,38,:,:)
 vwnd = f0->V(:,38,:,:)
 
 uwnd_old = uwnd
 vwnd_old = vwnd

 days	= (/"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"/)

do i = 0,30 
 print(i+"test")
 f = addfile("./OLD/MERRA2.201501"+days(i)+".A3dyn.2x25.nc4", "r")

 uwnd = f->U(:,38,:,:)
 vwnd = f->V(:,38,:,:)

 uwnd&time = uwnd&time + (i+1)*1440

 uwnd_new = array_append_record (uwnd_old, uwnd, 0)
 vwnd_new = array_append_record (vwnd_old, vwnd, 0)

 delete(uwnd_old)
 delete(vwnd_old)
 uwnd_old = uwnd_new
 vwnd_old = vwnd_new

 delete(uwnd_new)
 delete(vwnd_new)

end do


 f1	= addfile("./OLD/MERRA2.20150201.A3dyn.2x25.nc4","r")

 uwnd = f1->U(:,38,:,:)
 vwnd = f1->V(:,38,:,:)

 uwnd&time = uwnd&time + 32*1440

 uwnd_new = array_append_record (uwnd_old, uwnd, 0)
 vwnd_new = array_append_record (vwnd_old, vwnd, 0)

 delete(uwnd_old)
 delete(vwnd_old)
 uwnd_old = uwnd_new
 vwnd_old = vwnd_new

 delete(uwnd_new)
 delete(vwnd_new)

 time         = uwnd_old&time
 Nt           = dimsizes(time)
 tNew         = ispan(time(0),time(Nt-1),10)
 copy_VarMeta(time0(0),tNew(0))
 tNew@delta_t = "0000-00-00 00:10:00"
 tNew@time_increment = "001000"

print(tNew(135:dimsizes(tNew)-135)-1440)

;===============================================================
; initial time is 20141231-01:30, time step is 10min.
;==

 U_new = linint1_n_Wrap(time, uwnd_old, False, tNew, 0, 0)
 V_new = linint1_n_Wrap(time, vwnd_old, False, tNew, 0, 0)

 system("/bin/rm -f UV_wind.nc")    ; remove any pre-existing file
 ncdf = addfile("UV_wind.nc" ,"c")  ; open output netCDF file

 ncdf->time  = tNew(135:dimsizes(tNew)-135)
 ncdf->lat   = lat
 ncdf->lon   = lon

 ncdf->uwnd = U_new(135:dimsizes(tNew)-135,:,:)
 ncdf->vwnd = V_new(135:dimsizes(tNew)-135,:,:)


