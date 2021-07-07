load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin

  print("Started: "+systemfunc("date"))
; Begin user modification area
; IOP1 runs from 0000 UTC 2014/02/15 through 0000 UTC 2014/03/26, which is too
; long for the CCPP SCM, so select a start date in format YYYY, MM, DD 
; (Note: YYYY must be 2014 for IOP1, and start date needs to be the same for 
; both DEF and SCM files.)
  YYYY = "2014"
  MM   = "03"
  DD   = "05"

; Add directory and file info
  arm_dir = "/glade/work/damico/"
  iops    = "mao180varanaecmwfM1.c1.20140201.000000.cdf"

; Add the IOP file
  goam = addfile(arm_dir+iops,"r")

; Subset the times, the current setup has a four day simulation starting at 
; 0900 UTC (~sunrise) and ending at sunrise of the final day
  start_date = 0.
  sel_date   = 0.

  start_date@units = "seconds since 2014-02-15 00:00:00"
  sel_date@units   = "seconds since "+YYYY+"-"+MM+"-"+DD+" 09:00:00"

  sel_date_on_start = cd_convert(sel_date,start_date@units)

  end_date = 0.

  if (MM .eq. "02" .and. toint(DD) .gt. 23) then
    end_date@units = "seconds since "+YYYY+"-0"+(toint(MM)+1)+"-0"+ \
                     ((toint(DD)+4)-28)+" 09:00:00"
  else if (toint(DD) .lt. 6) then
    end_date@units = "seconds since "+YYYY+"-"+MM+"-0"+(toint(DD)+4)+ \
                     " 09:00:00"
  else if (MM .eq. "03" .and. toint(DD) .ge. 22) then
    end_date@units = "seconds since "+YYYY+"-03-26 00:00:00"
  else
    end_date@units = "seconds since "+YYYY+"-"+MM+"-"+(toint(DD)+4)+ \
                     " 09:00:00"
  end if
  end if
  end if

  end_date_on_start = cd_convert(end_date,start_date@units)

  time_total = goam->time_offset

  tm_s = closest_val(sel_date_on_start,time_total)
  tm_f = closest_val(end_date_on_start,time_total)

; Define the initial time
  t0           = todouble(0.)
  t0@long_name = "Initial time"
  t0@units     = sel_date@units
  t0@calendar  = "Gregorian"

  nt0 = dimsizes(t0)

; Add latitude and longitude
  latM = goam->lat
  lonM = goam->lon

  lat = todouble(latM) ; netCDF variables should all be double
  lon = todouble(lonM * -1.) ; Convert longitude to degrees east

  lat@long_name = "latitude"
  lat@units     = "degrees_north"
  lon@long_name = "longitude"
  lon@units     = "degrees_east"

  ny = dimsizes(lat)
  nx = dimsizes(lon)

; Levels
  lev_hPa = goam->lev

  lev_pressure           = todouble(lev_hPa * 100.) ; Convert to Pa
  lev_pressure@long_name = "pressure for variable pressure"
  lev_pressure@units     = "Pa"

  np = dimsizes(lev_pressure)

  pressure = new((/ nt0,np,ny,nx /),double)

  pressure(0,:,0,0)  = lev_pressure
  pressure@long_name = "Pressure"
  pressure@units     = "Pa"
  pressure!1         = "lev_pressure"

; Surface pressure
  psurf = new((/ nt0,ny,nx /),double)
  psurf = goam->p_srf_aver(tm_s)

  ps           = psurf * 100.
  ps@long_name = "Surface pressure"
  ps@units     = "Pa"

; Initial profiles (u & v)
  lev_u           = lev_pressure
  lev_u@long_name = "pressure for variable u"

  u = new((/ nt0,np,ny,nx /),double)

  u(0,:,0,0)  = goam->u(tm_s,:)
  u@long_name = "Zonal wind"
  u@units     = "m s-1"
  u!1         = "lev_u"

  lev_v           = lev_pressure
  lev_v@long_name = "pressure for variable v"

  v = new((/ nt0,np,ny,nx /),double)

  v(0,:,0,0)  = goam->v(tm_s,:)
  v@long_name = "Meridional wind"
  v@units     = "m s-1"
  v!1         = "lev_v"

; Initial temperature profile
  lev_temp           = lev_pressure
  lev_temp@long_name = "pressure for variable temp"

  temp = new((/ nt0,np,ny,nx /),double)

  temp(0,:,0,0)  = goam->T(tm_s,:)
  temp@long_name = "Temperature"
  temp@units     = "K"
  temp!1         = "lev_temp"

; Initial profile of rv
  lev_rv           = lev_pressure
  lev_rv@long_name = "pressure for variable rv"

  rvap = new((/ nt0,np,ny,nx /),double)

  rvap(0,:,0,0) = goam->q(tm_s,:)
  rv            = rvap / 1000. ; Convert to kg kg-1
  rv@long_name  = "Water vapor mixing ratio"
  rv@units      = "kg kg-1"
  rv!1          = "lev_rv"

; Forcing profiles (start with omega)
  time_omega           = time_total(tm_s:tm_f)
  time_omega@long_name = "Forcing time for omega"
  time_omega@units     = t0@units
  time_omega@calendar  = t0@calendar

  nt = dimsizes(time_omega)

  lev_omega           = lev_pressure
  lev_omega@long_name = "pressure for variable omega"

  omega_mbh = new((/ nt,np,ny,nx /),double)

  omega_mbh(:,:,0,0) = goam->omega(tm_s:tm_f,:)
  omega              = omega_mbh / 60. / 60. * 100. ; Convert from hPa hr-1 to
                                                    ; Pa s-1
  omega@long_name    = "Vertical pressure velovity"
  omega@units        = "Pa s-1"
  omega!0            = "time_omega"
  omega!1            = "lev_omega"

; Temperature advection (horizontal forcing)
  time_temp_adv           = time_omega
  time_temp_adv@long_name = "Forcing time for temp_adv"

  lev_temp_adv           = lev_pressure
  lev_temp_adv@long_name = "pressure for temp_adv"

  temp_Khr = new((/ nt,np,ny,nx /),double)

  temp_Khr(:,:,0,0)  = goam->T_adv_h(tm_s:tm_f,:)
  temp_adv           = temp_Khr / 60. / 60. ; Convert to K s-1
  temp_adv@long_name = "Temperature large-scale advection"
  temp_adv@units     = "K s-1"
  temp_adv!0         = "time_temp_adv"
  temp_adv!1         = "lev_temp_adv"

; rv advection
  time_rv_adv           = time_omega
  time_rv_adv@long_name = "Forcing time for rv_adv"

  lev_rv_adv           = lev_pressure
  lev_rv_adv@long_name = "pressure for rv_adv"

  rv_gkghr = new((/ nt,np,ny,nx /),double)

  rv_gkghr(:,:,0,0) = goam->q_adv_h(tm_s:tm_f,:)
  rv_adv            = rv_gkghr / 60. / 60. / 1000. ; Convert to kg kg-1 s-1
  rv_adv@long_name  = "Water vapor mixing ratio large-scale advection"
  rv_adv@units      = "kg kg-1 s-1"
  rv_adv!0          = "time_rv_adv"
  rv_adv!1          = "lev_rv_adv"

;; u & v nudging profiles (no u & v horizontal advection available)
;  time_u_nudging           = time_omega
;  time_u_nudging@long_name = "Forcing time for u_nudging"
;
;  lev_u_nudging           = lev_pressure
;  lev_u_nudging@long_name = "pressure for u_nudging"
;
;  u_nudging = new((/ nt,np,ny,nx /),double)
;
;  u_nudging(:,:,0,0)  = goam->u(tm_s:tm_f,:)
;  u_nudging@long_name = "Zonal wind profile for nudging"
;  u_nudging@units     = u@units
;  u_nudging!0         = "time_u_nudging"
;  u_nudging!1         = "lev_u_nudging"
;
;  time_v_nudging           = time_omega
;  time_v_nudging@long_name = "Forcing time for v_nudging"
;
;  lev_v_nudging           = lev_pressure
;  lev_v_nudging@long_name = "pressure for v_nudging"
;
;  v_nudging = new((/ nt,np,ny,nx /),double)
;
;  v_nudging(:,:,0,0)  = goam->v(tm_s:tm_f,:)
;  v_nudging@long_name = "Meridional wind profile for nudging"
;  v_nudging@units     = v@units
;  v_nudging!0         = "time_v_nudging"
;  v_nudging!1         = "lev_v_nudging"
;
; Surface temperature, LHF, and SHF
  time_ts           = time_omega
  time_ts@long_name = "Forcing time for ts"

  tsc = new((/ nt,ny,nx /),double)

  tsc(:,0,0)   = goam->T_skin(tm_s:tm_f)
  ts           = tsc + 273.15 ; Convert to K
  ts@long_name = "Surface temperature"
  ts@units     = "K"
  ts!0         = "time_ts"

  time_sfc_sens_flx           = time_omega
  time_sfc_sens_flx@long_name = "Forcing time for sfc_sens_flx"

  sfc_sens_flx = new((/ nt,ny,nx /),double)

  sfc_sens_flx(:,0,0)    = goam->SH(tm_s:tm_f)
  sfc_sens_flx@long_name = "Surface sensible heat flux (positive upward)"
  sfc_sens_flx@units     = "W m-2"
  sfc_sens_flx!0         = "time_sfc_sens_flx"

  time_sfc_lat_flx           = time_omega
  time_sfc_lat_flx@long_name = "Forcing time for sfc_lat_flx"

  sfc_lat_flx = new((/ nt,ny,nx /),double)

  sfc_lat_flx(:,0,0)    = goam->LH(tm_s:tm_f)
  sfc_lat_flx@long_name = "Surface latent heat flux (positive upward)"
  sfc_lat_flx@units     = "W m-2"
  sfc_lat_flx!0         = "time_sfc_lat_flx"

; END OF ADDING VARIABLES, begin writing netCDF
;  setfileoption("nc","FileStructure","Advanced")
;  setfileoption("nc","Format","NetCDF4") ;--> NOTE: DEPHY prefers netCDF3

  fout = "GOAMAZON_IOP1_DEF_driver.nc"
  system("/bin/rm -fv "+fout)
  fo   = addfile(fout,"c")

  setfileoption(fo,"DefineMode",True)

  fAtt                    = True
  fAtt@case               = "GOAMAZON/IOP1"
  fAtt@title              = "Forcing and initial conditions for GOAMAZON/IOP1"
  fAtt@reference          = "Tang et al. (2016,Atmos. Chem. Phys.), "+ \
                            "Martin et al. (2017,BAMS)"
  fAtt@author             = "D. DAmico (DEF File), Xie and Tang (Forcing "+ \
                            "and IC File)"
  fAtt@comment            = ""
  fAtt@startDate          = YYYY+MM+DD+"0900"
  if (MM .eq. "02" .and. toint(DD) .gt. 23) then
    fAtt@endDate          = YYYY+"0"+(toint(MM)+1)+"0"+((toint(DD)+4)-28)+ \
                            "0900"
  else if (toint(DD) .lt. 6) then
    fAtt@endDate          = YYYY+MM+"0"+(toint(DD)+4)+"0900"
  else if (MM .eq. "03" .and. toint(DD) .ge. 22) then
    fAtt@endDate          = YYYY+"03260000"
  else
    fAtt@endDate          = YYYY+MM+(toint(DD)+4)+"0900"
  end if
  end if
  end if
  fAtt@adv_temp           = 1
  fAtt@adv_theta          = 0
  fAtt@adv_thetal         = 0
  fAtt@rad_temp           = "adv"
  fAtt@rad_theta          = 0
  fAtt@rad_thetal         = 0
  fAtt@adv_qv             = 0
  fAtt@adv_qt             = 0
  fAtt@adv_rv             = 1
  fAtt@adv_rt             = 0
  fAtt@forc_w             = 0
  fAtt@forc_omega         = 1
  fAtt@forc_geo           = 0
  fAtt@nudging_u          = 0 ;10800
  fAtt@nudging_v          = 0 ;10800
  fAtt@nudging_temp       = 0
  fAtt@nudging_theta      = 0
  fAtt@nudging_thetal     = 0
  fAtt@nudging_qv         = 0
  fAtt@nudging_qt         = 0
  fAtt@nudging_rv         = 0
  fAtt@nudging_rt         = 0
;  fAtt@p_nudging_u        = 100785.5
;  fAtt@p_nudging_v        = 100785.5
  fAtt@zorog              = goam->alt
  fAtt@z0                 = 0.02 ; See Carneiro et al. 2021 (might want to run 
                                 ; with a variety of z0 values)
  fAtt@surfaceType        = "land"
  fAtt@surfaceForcing     = "surfaceFlux"
  fAtt@surfaceForcingWind = "z0"

  fileattdef(fo,fAtt)

;; Define the dimensions (w/ u & v nudging)
;  dimNames = (/ "t0","lat","lon","lev_pressure","lev_u","lev_v","lev_temp", \
;                "lev_rv","time_omega","lev_omega","time_temp_adv", \
;                "lev_temp_adv","time_rv_adv","lev_rv_adv","time_u_nudging", \
;                "lev_u_nudging","time_v_nudging","lev_v_nudging","time_ts", \
;                "time_sfc_sens_flx","time_sfc_lat_flx" /)
;  dimSizes = (/ nt0,ny,nx,np,np,np,np,np,nt,np,nt,np,nt,np,nt,np,nt,np \
;                nt,nt,nt /)
;  dimUnlim = (/ False,False,False,False,False,False,False,False,False, \
;                False,False,False,False,False,False,False,False,False, \
;                False,False,False /)
; Define the dimensions (w/o u & v nudging)
  dimNames = (/ "t0","lat","lon","lev_pressure","lev_u","lev_v","lev_temp", \
                "lev_rv","time_omega","lev_omega","time_temp_adv", \
                "lev_temp_adv","time_rv_adv","lev_rv_adv","time_ts", \
                "time_sfc_sens_flx","time_sfc_lat_flx" /)
  dimSizes = (/ nt0,ny,nx,np,np,np,np,np,nt,np,nt,np,nt,np, \
                nt,nt,nt /)
  dimUnlim = (/ False,False,False,False,False,False,False,False,False, \
                False,False,False,False,False,False,False,False /)

  filedimdef(fo,dimNames,dimSizes,dimUnlim)

; Add the dimensions, and the variables that use them (attempting to mix them
; to match the DEF netCDF file for GABLS1)
  filevardef(fo,"t0",typeof(t0),"t0")
  filevarattdef(fo,"t0",t0)
  fo->t0 = (/ t0 /)

  filevardef(fo,"lat",typeof(lat),"lat")
  filevarattdef(fo,"lat",lat)
  fo->lat = (/ lat /)

  filevardef(fo,"lon",typeof(lon),"lon")
  filevarattdef(fo,"lon",lon)
  fo->lon = (/ lon /)

  filevardef(fo,"ps",typeof(ps),(/ "t0","lat","lon" /))
  filevarattdef(fo,"ps",ps)
  fo->ps = (/ ps /)

  filevardef(fo,"lev_pressure",typeof(lev_pressure),"lev_pressure")
  filevarattdef(fo,"lev_pressure",lev_pressure)
  fo->lev_pressure = (/ lev_pressure /)

  filevardef(fo,"pressure",typeof(pressure),(/ "t0","lev_pressure","lat", \
                                               "lon" /))
  filevarattdef(fo,"pressure",pressure)
  fo->pressure = (/ pressure /)

  filevardef(fo,"lev_u",typeof(lev_u),"lev_u")
  filevarattdef(fo,"lev_u",lev_u)
  fo->lev_u = (/ lev_u /)

  filevardef(fo,"u",typeof(u),(/ "t0","lev_u","lat","lon" /))
  filevarattdef(fo,"u",u)
  fo->u = (/ u /)

  filevardef(fo,"lev_v",typeof(lev_v),"lev_v")
  filevarattdef(fo,"lev_v",lev_v)
  fo->lev_v = (/ lev_v /)

  filevardef(fo,"v",typeof(v),(/ "t0","lev_v","lat","lon" /))
  filevarattdef(fo,"v",v)
  fo->v = (/ v /)

  filevardef(fo,"lev_temp",typeof(lev_temp),"lev_temp")
  filevarattdef(fo,"lev_temp",lev_temp)
  fo->lev_temp = (/ lev_temp /)

  filevardef(fo,"temp",typeof(temp),(/ "t0","lev_temp","lat","lon" /))
  filevarattdef(fo,"temp",temp)
  fo->temp = (/ temp /)

  filevardef(fo,"lev_rv",typeof(lev_rv),"lev_rv")
  filevarattdef(fo,"lev_rv",lev_rv)
  fo->lev_rv = (/ lev_rv /)

  filevardef(fo,"rv",typeof(rv),(/ "t0","lev_rv","lat","lon" /))
  filevarattdef(fo,"rv",rv)
  fo->rv = (/ rv /)

  filevardef(fo,"time_omega",typeof(time_omega),"time_omega")
  filevarattdef(fo,"time_omega",time_omega)
  fo->time_omega = (/ time_omega /)

  filevardef(fo,"lev_omega",typeof(lev_omega),"lev_omega")
  filevarattdef(fo,"lev_omega",lev_omega)
  fo->lev_omega = (/ lev_omega /)

  filevardef(fo,"omega",typeof(omega),(/ "time_omega","lev_omega","lat", \
                                         "lon" /))
  filevarattdef(fo,"omega",omega)
  fo->omega = (/ omega /)

  filevardef(fo,"time_temp_adv",typeof(time_temp_adv),"time_temp_adv")
  filevarattdef(fo,"time_temp_adv",time_temp_adv)
  fo->time_temp_adv = (/ time_temp_adv /)

  filevardef(fo,"lev_temp_adv",typeof(lev_temp_adv),"lev_temp_adv")
  filevarattdef(fo,"lev_temp_adv",lev_temp_adv)
  fo->lev_temp_adv = (/ lev_temp_adv /)

  filevardef(fo,"temp_adv",typeof(temp_adv),(/ "time_temp_adv","lev_temp_adv",\
                                               "lat","lon" /))
  filevarattdef(fo,"temp_adv",temp_adv)
  fo->temp_adv = (/ temp_adv /)

  filevardef(fo,"time_rv_adv",typeof(time_rv_adv),"time_rv_adv")
  filevarattdef(fo,"time_rv_adv",time_rv_adv)
  fo->time_rv_adv = (/ time_rv_adv /)

  filevardef(fo,"lev_rv_adv",typeof(lev_rv_adv),"lev_rv_adv")
  filevarattdef(fo,"lev_rv_adv",lev_rv_adv)
  fo->lev_rv_adv = (/ lev_rv_adv /)

  filevardef(fo,"rv_adv",typeof(rv_adv),(/ "time_rv_adv","lev_rv_adv", \
                                           "lat","lon" /))
  filevarattdef(fo,"rv_adv",rv_adv)
  fo->rv_adv = (/ rv_adv /)

;  filevardef(fo,"time_u_nudging",typeof(time_u_nudging),"time_u_nudging")
;  filevarattdef(fo,"time_u_nudging",time_u_nudging)
;  fo->time_u_nudging = (/ time_u_nudging /)

;  filevardef(fo,"lev_u_nudging",typeof(lev_u_nudging),"lev_u_nudging")
;  filevarattdef(fo,"lev_u_nudging",lev_u_nudging)
;  fo->lev_u_nudging = (/ lev_u_nudging /)

;  filevardef(fo,"u_nudging",typeof(u_nudging),(/ "time_u_nudging", \
;                                                 "lev_u_nudging","lat", \
;                                                 "lon" /))
;  filevarattdef(fo,"u_nudging",u_nudging)
;  fo->u_nudging = (/ u_nudging /)

;  filevardef(fo,"time_v_nudging",typeof(time_v_nudging),"time_v_nudging")
;  filevarattdef(fo,"time_v_nudging",time_v_nudging)
;  fo->time_v_nudging = (/ time_v_nudging /)

;  filevardef(fo,"lev_v_nudging",typeof(lev_v_nudging),"lev_v_nudging")
;  filevarattdef(fo,"lev_v_nudging",lev_v_nudging)
;  fo->lev_v_nudging = (/ lev_v_nudging /)

;  filevardef(fo,"v_nudging",typeof(v_nudging),(/ "time_v_nudging", \
;                                                 "lev_v_nudging","lat", \
;                                                 "lon" /))
;  filevarattdef(fo,"v_nudging",u_nudging)
;  fo->v_nudging = (/ v_nudging /)

  filevardef(fo,"time_ts",typeof(time_ts),"time_ts")
  filevarattdef(fo,"time_ts",time_ts)
  fo->time_ts = (/ time_ts /)

  filevardef(fo,"ts",typeof(ts),(/ "time_ts","lat","lon" /))
  filevarattdef(fo,"ts",ts)
  fo->ts = (/ ts /)

  filevardef(fo,"time_sfc_sens_flx",typeof(time_sfc_sens_flx), \
            "time_sfc_sens_flx")
  filevarattdef(fo,"time_sfc_sens_flx",time_sfc_sens_flx)
  fo->time_sfc_sens_flx = (/ time_sfc_sens_flx /)

  filevardef(fo,"sfc_sens_flx",typeof(sfc_sens_flx), (/ "time_sfc_sens_flx", \
                                                        "lat","lon" /))
  filevarattdef(fo,"sfc_sens_flx",sfc_sens_flx)
  fo->sfc_sens_flx = (/ sfc_sens_flx /)

  filevardef(fo,"time_sfc_lat_flx",typeof(time_sfc_lat_flx), \
             "time_sfc_lat_flx")
  filevarattdef(fo,"time_sfc_lat_flx",time_sfc_lat_flx)
  fo->time_sfc_lat_flx = (/ time_sfc_lat_flx /)

  filevardef(fo,"sfc_lat_flx",typeof(sfc_lat_flx), (/ "time_sfc_lat_flx", \
                                                      "lat","lon" /))
  filevarattdef(fo,"sfc_lat_flx",sfc_lat_flx)
  fo->sfc_lat_flx = (/ sfc_lat_flx /)

  print("Created and wrote new file "+fout)
  print("Complete: "+systemfunc("date"))

end

