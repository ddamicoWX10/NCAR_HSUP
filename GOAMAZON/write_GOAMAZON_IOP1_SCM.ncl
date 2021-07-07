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

; Add directory and file information
  arm_dir = "/glade/work/damico/"
  iops    = "mao180varanaecmwfM1.c1.20140201.000000.cdf"

; Constants
  R_dry = 287.    ; ideal gas constant for dry air (J kg^-1 K^-1)
  R_vap = 461.5   ; gas constant for water vapor (J kg^-1 K^-1)
  c_p   = 1004.   ; specific heat at const pressure (J kg^-1 K^-1)
  l_v   = 2.5e6   ; latent heat of vaporization (J kg^-1)
  g     = 9.80665 ; acceleration due tp gravity (m s^-2)
  P00   = 100000. ; initial pressure (Pa)

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

  hgt_int = ispan(0,6000,10)

  height           = new((/ nt0,dimsizes(hgt_int),ny,nx /),double)
  height(0,:,0,0)  = todouble(hgt_int)
  height@long_name = "Height above ground"
  height@units     = "m"

  np = dimsizes(height(0,:,0,0))

; Pressure levels
  pall = goam->lev

  p_in           = todouble(pall * 100.)
  p_in@long_name = "pressure"
  p_in@units     = "Pa"

  npi = dimsizes(p_in)

  temp0           = goam->T(tm_s,:)
  temp0@long_name = "Temperature"
  temp0@units     = "K"

  rvap0         = goam->q(tm_s,:)
  rv0           = rvap0 / 1000.
  rv0@long_name = "Water vapor mixing ratio"
  rv0@units     = "kg kg-1"

  tv = temp0
  tv = temp0 * ((1 + (rv0 / (R_dry/R_vap))) / (1 + rv0))

  hght    = new((/ npi /),double)
  hght(0) = 0.
  do i=1,npi-1,1
    hght(i) = hght(i-1) + (((R_dry * ((tv(i) + tv(i-1)) * 0.5)) / g) * \
              log(p_in(i-1) / p_in(i)))
  end do

  lev           = linint1(hght,p_in,False,height(0,:,0,0),0)
  lev@long_name = "pressure"
  lev@units     = "Pa"

  pressure          = new((/ nt0,np,ny,nx /),double)
  pressure(0,:,0,0) = lev

; Surface pressure
  psurf = new((/ nt0,ny,nx /),double)
  psurf = goam->p_srf_aver(tm_s)

  ps           = psurf * 100.
  ps@long_name = "Surface pressure"
  ps@units     = "Pa"

; Initial profiles (u & v)
  ui = goam->u(tm_s,:)

  uz = linint1(hght,ui,False,height(0,:,0,0),0)

  u           = new((/ nt0,np,ny,nx /),double)
  u(0,:,0,0)  = uz
  u@long_name = "Zonal wind"
  u@units     = "m s-1"

  vi = goam->v(tm_s,:)

  vz = linint1(hght,vi,False,height(0,:,0,0),0)

  v = new((/ nt0,np,ny,nx /),double)

  v(0,:,0,0)  = vz
  v@long_name = "Meridional wind"
  v@units     = "m s-1"

; Initial temperature profile
  tempz = linint1(hght,temp0,False,height(0,:,0,0),0)

  temp = new((/ nt0,np,ny,nx /),double)

  temp(0,:,0,0)  = tempz
  temp@long_name = "Temperature"
  temp@units     = "K"

; From temperature, calculate potential temperature (needed in DEPHY format)
  theta           = temp * ((P00 / pressure) ^ (R_dry / c_p))
  theta@long_name = "Potential temperature"
  theta@units     = "K"

; Initial profile of rv
  rvz = linint1(hght,rv0,False,height(0,:,0,0),0)

  rv = new((/ nt0,np,ny,nx /),double)

  rv           = rvz / 1000.
  rv@long_name = "Water vapor mixing ratio"
  rv@units     = "kg kg-1"

; Calculate water vapor specific humidity (qv) from rv
  qv           = rv / (1. + rv)
  qv@long_name = "Specific humidity"
  qv@units     = "kg kg-1"

; Total water specific humidity (qt) and mixing ratio (rt). Going to be equal
; to qv and rv respectively because qi,ri,ql, and rl are undefined in dataset.
  qt           = qv
  qt@long_name = "Total water content"
  rt           = rv
  rt@long_name = "Total water mixing ratio"

; Initial profiles of ql,qi,rl,ri, and tke (each will be all zeroes, as they
; are undefined in the dataset)
  ql           = new((/ nt0,np,ny,nx /),double)
  ql           = 0.
  ql@long_name = "Liquid water content"
  ql@units     = "kg kg-1"

  qi           = new((/ nt0,np,ny,nx /),double)
  qi           = 0.
  qi@long_name = "Ice water content"
  qi@units     = "kg kg-1"

  rl           = new((/ nt0,np,ny,nx /),double)
  rl           = 0.
  rl@long_name = "Liquid water mixing ratio"
  rl@units     = "kg kg-1"

  ri           = new((/ nt0,np,ny,nx /),double)
  ri           = 0.
  ri@long_name = "Ice water mixing ratio"
  ri@units     = "kg kg-1"

  tke           = new((/ nt0,np,ny,nx /),double)
  tke           = 0.
  tke@long_name = "Turbulent kinetic energy"
  tke@units     = "m2 s-2"

; Time for forcing
  time           = time_total(tm_s:tm_f)
  time@long_name = "Forcing time"
  time@units     = t0@units

  nt = dimsizes(time)

; Surface pressure for forcing
  ps_f              = new((/ nt,ny,nx /),double)
  ps_forc           = new((/ nt,ny,nx /),double)
  ps_f(:,0,0)       = goam->p_srf_aver(tm_s:tm_f)
  ps_forc           = ps_f * 1.
  ps_forc@long_name = "Surface pressure for forcing"

; Height for forcing
  height_forc = new((/ nt,np,ny,nx /),double)
  do i=0,np-1,1
    height_forc(:,i,:,:) = height(:,i,:,:)
  end do
  height_forc@long_name = "Height for forcing"

; Pressure for forcing
  pressure_forc = new((/ nt,np,ny,nx /),double)
  do i=0,np-1,1
    pressure_forc(:,i,:,:) = lev(i)
  end do
  pressure_forc@long_name = "Pressure for forcing"

; Forcing profiles (start with omega)
  omega_mbh          = new((/ nt,npi,ny,nx /),double)
  omega_mbh(:,:,0,0) = goam->omega(tm_s:tm_f,:)

  omega_z = new((/ nt,np,ny,nx /),double)
  do i=0,nt-1,1
    omega_z(i,:,0,0) = linint1(hght,omega_mbh(i,:,0,0),False,height(0,:,0,0), \
                               0)
  end do

  omega              = omega_z / 60. / 60. * 100.
  omega@long_name    = "Vertical pressure velocity"
  omega@units        = "Pa s-1"

; Temperature advection (horizontal forcing)
  temp_Khr          = new((/ nt,npi,ny,nx /),double)
  temp_Khr(:,:,0,0) = goam->T_adv_h(tm_s:tm_f,:)

  temp_advZ = new((/ nt,np,ny,nx /),double)
  do i=0,nt-1,1
    temp_advZ(i,:,0,0) = linint1(hght,temp_Khr(i,:,0,0),False, \
                                 height(0,:,0,0),0)
  end do

  temp_adv           = temp_advZ / 60. / 60.
  temp_adv@long_name = "Temperature large-scale advection"
  temp_adv@units     = "K s-1"

; rv advection
  rv_gkghr          = new((/ nt,npi,ny,nx /),double)
  rv_gkghr(:,:,0,0) = goam->q_adv_h(tm_s:tm_f,:)

  rv_advZ = new((/ nt,np,ny,nx /),double)
  do i=0,nt-1,1
    rv_advZ(i,:,0,0) = linint1(hght,rv_gkghr(i,:,0,0),False,height(0,:,0,0), \
                               0)
  end do

  rv_adv            = rv_advZ / 60. / 60. / 1000.
  rv_adv@long_name  = "Water vapor mixing ratio large-scale advection"
  rv_adv@units      = "kg kg-1 s-1"

;; u & v nudging profiles (no u & v horizontal advection available)
;  u_nudgingZ          = new((/ nt,npi,ny,nx /),double)
;  u_nudgingZ(:,:,0,0) = goam->u(tm_s:tm_f,:)

;  u_nudging = new((/ nt,np,ny,nx /),double)
;  do i=0,nt-1,1
;    u_nudging(i,:,0,0) = linint1(hght,u_nudgingZ(i,:,0,0),False, \
;                                 height(0,:,0,0),0)
;  end do

;  u_nudging@long_name = "Zonal wind profile for nudging"
;  u_nudging@units     = u@units

;  v_nudgingZ          = new((/ nt,npi,ny,nx /),double)
;  v_nudgingZ(:,:,0,0) = goam->v(tm_s:tm_f,:)

;  v_nudging = new((/ nt,np,ny,nx /),double)
;  do i=0,nt-1,1
;    v_nudging(i,:,0,0) = linint1(hght,v_nudgingZ(i,:,0,0),False, \
;                                 height(0,:,0,0),0)
;  end do

;  v_nudging@long_name = "Meridional wind profile for nudging"
;  v_nudging@units     = v@units

; Surface temperature, LHF, and SHF
  tsc = new((/ nt,ny,nx /),double)

  tsc(:,0,0)   = goam->T_skin(tm_s:tm_f)
  ts           = tsc + 273.15
  ts@long_name = "Surface temperature"
  ts@units     = "K"

  sfc_shf      = new((/ nt,ny,nx /),double)
  sfc_sens_flx = new((/ nt,ny,nx /),double)

  sfc_shf(:,0,0)         = goam->SH(tm_s:tm_f)
  sfc_sens_flx           = sfc_shf * 1.
  sfc_sens_flx@long_name = "Surface sensible heat flux (positive upward)"
  sfc_sens_flx@units     = "W m-2"

  sfc_lhf     = new((/ nt,ny,nx /),double)
  sfc_lat_flx = new((/ nt,ny,nx /),double)

  sfc_lhf(:,0,0)        = goam->LH(tm_s:tm_f)
  sfc_lat_flx           = sfc_lhf * 1.
  sfc_lat_flx@long_name = "Surface latent heat flux (positive upward)"
  sfc_lat_flx@units     = "W m-2"

; END OF ADDING VARIABLES, begin writing netCDF
;  setfileoption("nc","FileStructure","Advanced")
;  setfileoption("nc","Format","NetCDF4") ;--> NOTE: DEPHY prefers netCDF3

  fout = "GOAMAZON_IOP1_SCM_driver.nc"
  system("/bin/rm -fv "+fout)
  fo   = addfile(fout,"c")

  setfileoption(fo,"DefineMode",True)

  fAtt                    = True
  fAtt@case               = "GOAMAZON/IOP1"
  fAtt@title              = "Forcing and initial conditions for GOAMAZON/IOP1"
  fAtt@reference          = "Tang et al. (2016,Atmos. Chem. Phys.), "+ \
                            "Martin et al. (2017,BAMS)"
  fAtt@author             = "D. DAmico (SCM File), Xie and Tang (Forcing "+ \
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

; Define the dimensions
  dimNames = (/ "t0","lat","lon","lev","time" /)
  dimSizes = (/  nt0, ny,   nx,   np,   nt    /)
  dimUnlim = (/ False,False,False,False,False /)

  filedimdef(fo,dimNames,dimSizes,dimUnlim)

; Add the dimensions, and the variables that use them (attempting to mix them
; to match the SCM netCDF file for GABLS1)
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

  filevardef(fo,"lev",typeof(lev),"lev")
  filevarattdef(fo,"lev",lev)
  fo->lev = (/ lev /)

  filevardef(fo,"pressure",typeof(pressure),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"pressure",pressure)
  fo->pressure = (/ pressure /)

  filevardef(fo,"height",typeof(height),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"height",height)
  fo->height = (/ height /)

  filevardef(fo,"u",typeof(u),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"u",u)
  fo->u = (/ u /)

  filevardef(fo,"v",typeof(v),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"v",v)
  fo->v = (/ v /)

  filevardef(fo,"temp",typeof(temp),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"temp",temp)
  fo->temp = (/ temp /)

  filevardef(fo,"theta",typeof(theta),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"theta",theta)
  fo->theta = (/ theta /)

  filevardef(fo,"qv",typeof(qv),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"qv",qv)
  fo->qv = (/ qv /)

  filevardef(fo,"qt",typeof(qt),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"qt",qt)
  fo->qt = (/ qt /)

  filevardef(fo,"rv",typeof(rv),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"rv",rv)
  fo->rv = (/ rv /)

  filevardef(fo,"rt",typeof(rt),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"rt",rt)
  fo->rt = (/ rt /)

  filevardef(fo,"rl",typeof(rl),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"rl",rl)
  fo->rl = (/ rl /)

  filevardef(fo,"ri",typeof(ri),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"ri",ri)
  fo->ri = (/ ri /)

  filevardef(fo,"ql",typeof(ql),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"ql",ql)
  fo->ql = (/ ql /)

  filevardef(fo,"qi",typeof(qi),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"qi",qi)
  fo->qi = (/ qi /)

  filevardef(fo,"tke",typeof(tke),(/ "t0","lev","lat","lon" /))
  filevarattdef(fo,"tke",tke)
  fo->tke = (/ tke /)

  filevardef(fo,"time",typeof(time),"time")
  filevarattdef(fo,"time",time)
  fo->time = (/ time /)

  filevardef(fo,"ps_forc",typeof(ps_forc),(/ "time","lat","lon" /))
  filevarattdef(fo,"ps_forc",ps_forc)
  fo->ps_forc = (/ ps_forc /)

  filevardef(fo,"pressure_forc",typeof(pressure_forc),(/ "time","lev", \
                                                         "lat","lon" /))
  filevarattdef(fo,"pressure_forc",pressure_forc)
  fo->pressure_forc = (/ pressure_forc /)

  filevardef(fo,"height_forc",typeof(height_forc),(/ "time","lev","lat", \
                                                     "lon" /))
  filevarattdef(fo,"height_forc",height_forc)
  fo->height_forc = (/ height_forc /)

  filevardef(fo,"omega",typeof(omega),(/ "time","lev","lat","lon" /))
  filevarattdef(fo,"omega",omega)
  fo->omega = (/ omega /)

  filevardef(fo,"temp_adv",typeof(temp_adv),(/ "time","lev","lat","lon" /))
  filevarattdef(fo,"temp_adv",temp_adv)
  fo->temp_adv = (/ temp_adv /)

  filevardef(fo,"rv_adv",typeof(rv_adv),(/ "time","lev","lat","lon" /))
  filevarattdef(fo,"rv_adv",rv_adv)
  fo->rv_adv = (/ rv_adv /)

;  filevardef(fo,"u_nudging",typeof(u_nudging),(/ "time","lev","lat", \
;                                                 "lon" /))
;  filevarattdef(fo,"u_nudging",u_nudging)
;  fo->u_nudging = (/ u_nudging /)

;  filevardef(fo,"v_nudging",typeof(v_nudging),(/ "time","lev","lat", \
;                                                 "lon" /))
;  filevarattdef(fo,"v_nudging",v_nudging)
;  fo->v_nudging = (/ v_nudging /)

  filevardef(fo,"ts",typeof(ts),(/ "time","lat","lon" /))
  filevarattdef(fo,"ts",ts)
  fo->ts = (/ ts /)

  filevardef(fo,"sfc_sens_flx",typeof(sfc_sens_flx), (/ "time","lat","lon" /))
  filevarattdef(fo,"sfc_sens_flx",sfc_sens_flx)
  fo->sfc_sens_flx = (/ sfc_sens_flx /)

  filevardef(fo,"sfc_lat_flx",typeof(sfc_lat_flx), (/ "time","lat","lon" /))
  filevarattdef(fo,"sfc_lat_flx",sfc_lat_flx)
  fo->sfc_lat_flx = (/ sfc_lat_flx /)

  print("Created and wrote new file "+fout)
  print("Complete: "+systemfunc("date"))

end

