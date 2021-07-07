load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; User modification area, the forcing and directories vary between
; LASSO simulations
  gtype = "pdf"

  pbl_list = (/ "satmedmf","ACM","YSU" /)

  dir_all = "/glade/u/home/damico/SCMv4p1_update/gmtb-scm/scm/bin/"+ \
            "output_wangara_day33_SCM_GFS_v15p2"
;            "output_wangara_day33_noForce_SCM_GFS_v15p2"

  filename = "output.nc"

; End user modification area

  ns = dimsizes(pbl_list)

; Constants
  R_dry = 287.    ; ideal gas constant for dry air (J kg^-1 K^-1)
  R_vap = 461.5   ; gas constant for water vapor (J kg^-1 K^-1)
  c_p   = 1004.   ; specific heat at const pressure (J kg^-1 K^-1)
  c_pv  = 1996.   ; spec heat of water vapor (J kg^-1 K^-1)
  l_v   = 2.5e6   ; latent heat of vaporization (J kg^-1)
  g     = 9.80665 ; acceleration due tp gravity (m s^-2)
  P00   = 100000. ; initial pressure (Pa)
  eps   = (R_dry / R_vap)
  pi    = 4.*atan(1.)

  print("Started: "+systemfunc("date"))

; Add SCM files
  fEDMF = addfile(dir_all+"/"+filename,"r")
  fACM  = addfile(dir_all+"_"+pbl_list(1)+"/"+filename,"r")
  fYSU  = addfile(dir_all+"_"+pbl_list(2)+"/"+filename,"r")

; Time
  tscm = fEDMF->time
  YYYY = fEDMF->init_year
  MM   = fEDMF->init_month
  DD   = fEDMF->init_day
  HH   = fEDMF->init_hour

  hscm       = tscm
  hscm       = tscm / 60. / 60.
  if (MM .lt. 10 .and. DD .lt. 10) then
    hscm@units = "hours since "+YYYY+"-0"+MM+"-0"+DD+"_"+HH+":00:00"
  else if (MM .lt. 10)
    hscm@units = "hours since "+YYYY+"-0"+MM+"-"+DD+"_"+HH+":00:00"
  else
    hscm@units = "hours since "+YYYY+"-"+MM+"-"+DD+"_"+HH+":00:00"
  end if
  end if

; SCM variables
  presctl_Pa = fEDMF->pres
  tmpctl_K   = fEDMF->T
  qvctl_kgkg = fEDMF->qv
  qcctl_kgkg = fEDMF->qc
  qcctl_conv = fEDMF->qc_conv
  qictl_kgkg = fEDMF->qi
  tsfcctl_K  = fEDMF->T_s
  psfcctl_Pa = fEDMF->pres_s
  prest1_Pa  = fACM->pres
  tmpt1_K    = fACM->T
  qvt1_kgkg  = fACM->qv
  qct1_kgkg  = fACM->qc
  qct1_conv  = fACM->qc_conv
  qit1_kgkg  = fACM->qi
  tsfct1_K   = fACM->T_s
  psfct1_Pa  = fACM->pres_s
  prest2_Pa  = fYSU->pres
  tmpt2_K    = fYSU->T
  qvt2_kgkg  = fYSU->qv
  qct2_kgkg  = fYSU->qc
  qct2_conv  = fYSU->qc_conv
  qit2_kgkg  = fYSU->qi
  tsfct2_K   = fYSU->T_s
  psfct2_Pa  = fYSU->pres_s

; Find SCM height (not an SCM output variable, use hypsometric equation)
; Start with whole atmosphere virtual temperature
  tvctl = tmpctl_K
  tvctl = tmpctl_K * ((1 + (qvctl_kgkg / (R_dry/R_vap))) / \
                     (1 + qvctl_kgkg))

  tvctl@description = "virtual temperature"

  tvt1 = tmpt1_K
  tvt1 = tmpt1_K * ((1 + (qvt1_kgkg / (R_dry/R_vap))) / \
                   (1 + qvt1_kgkg))

  tvt1@description = "virtual temperature"

  tvt2 = tmpt2_K
  tvt2 = tmpt2_K * ((1 + (qvt2_kgkg / (R_dry/R_vap))) / \
                   (1 + qvt2_kgkg))

  tvt2@description = "virtual temperature"

; Surface virtual temperature
  tvsctl = tsfcctl_K
  tvsctl = tsfcctl_K *((1 + (qvctl_kgkg(:,0,:) / (R_dry/R_vap))) / \
                      (1 + qvctl_kgkg(:,0,:)))

  tvsctl@description = "surface virtual temperature"

  tvst1 = tsfct1_K
  tvst1 = tsfct1_K *((1 + (qvt1_kgkg(:,0,:) / (R_dry/R_vap))) / \
                    (1 + qvt1_kgkg(:,0,:)))

  tvst1@description = "surface virtual temperature"

  tvst2 = tsfct2_K
  tvst2 = tsfct2_K *((1 + (qvt2_kgkg(:,0,:) / (R_dry/R_vap))) / \
                    (1 + qvt2_kgkg(:,0,:)))

  tvst2@description = "surface virtual temperature"

; Height of level 0
  zctl_0 = (((R_dry * tvsctl) / g) * log(psfcctl_Pa / presctl_Pa(:,0,:)))
  zt1_0  = (((R_dry * tvst1) / g) * log(psfct1_Pa / prest1_Pa(:,0,:)))
  zt2_0  = (((R_dry * tvst2) / g) * log(psfct2_Pa / prest2_Pa(:,0,:)))

; Use hypsometric equation, with height at level 0, to find full model height
  hgtctl_m      = new((/ dimsizes(prest1_Pa(:,:,0)) /),float)
  hgtt1_m       = new((/ dimsizes(prest1_Pa(:,:,0)) /),float)
  hgtt2_m       = new((/ dimsizes(prest1_Pa(:,:,0)) /),float)
  hgtctl_m(:,0) = zctl_0(:,0)
  hgtt1_m(:,0)  = zt1_0(:,0)
  hgtt2_m(:,0)  = zt2_0(:,0)
  do i=1,dimsizes(hgtt1_m(0,:))-1,1
    hgtctl_m(:,i) = hgtctl_m(:,i-1) + (((R_dry * tvctl(:,i,0)) / g) * \
                    log(presctl_Pa(:,i-1,0) / presctl_Pa(:,i,0)))
    hgtt1_m(:,i)  = hgtt1_m(:,i-1) + (((R_dry * tvt1(:,i,0)) / g) * \
                    log(prest1_Pa(:,i-1,0) / prest1_Pa(:,i,0)))
    hgtt2_m(:,i)  = hgtt2_m(:,i-1) + (((R_dry * tvt2(:,i,0)) / g) * \
                    log(prest2_Pa(:,i-1,0) / prest2_Pa(:,i,0)))
  end do
  hgtctl_m!1           = "vert_dim_layer"
  hgtctl_m@description = "model height"
  hgtctl_m@units       = "m"
  hgtt1_m!1            = "vert_dim_layer"
  hgtt1_m@description  = "model height"
  hgtt1_m@units        = "m"
  hgtt2_m!1            = "vert_dim_layer"
  hgtt2_m@description  = "model height"
  hgtt2_m@units        = "m"

; Calculating SCM saturation vapor pressure and saturation specific
; humidity
  svpctl_hPa = 6.11 * (10^((7.5*(tmpctl_K - 273.15)) / (273.3 + \
                          (tmpctl_K-273.15))))
  svpctl_Pa  = svpctl_hPa * 100.

  qvsatctl_kgkg = 0.622 * (svpctl_Pa / presctl_Pa)

  svpt1_hPa = 6.11 * (10^((7.5*(tmpt1_K - 273.15)) / (273.3 + \
                          (tmpt1_K-273.15))))
  svpt1_Pa  = svpt1_hPa * 100.

  qvsatt1_kgkg = 0.622 * (svpt1_Pa / prest1_Pa)

  svpt2_hPa = 6.11 * (10^((7.5*(tmpt2_K - 273.15)) / (273.3 + \
                          (tmpt2_K-273.15))))
  svpt2_Pa  = svpt2_hPa * 100.

  qvsatt2_kgkg = 0.622 * (svpt2_Pa / prest2_Pa)

; Calculating SCM relative humidity
  rhctl_127 = (qvctl_kgkg / qvsatctl_kgkg) * 100.

  rhctl_127@description = "Relative Humidity"
  rhctl_127@units       = "%"
  rhctl_127!0           = "Time"
  rhctl_127!1           = "bottom_top"
  rhctl_127!2           = "horiz_dim_layer"

  rht1_127 = (qvt1_kgkg / qvsatt1_kgkg) * 100.

  rht1_127@description = "Relative Humidity"
  rht1_127@units       = "%"
  rht1_127!0           = "Time"
  rht1_127!1           = "bottom_top"
  rht1_127!2           = "horiz_dim_layer"

  rht2_127 = (qvt2_kgkg / qvsatt2_kgkg) * 100.

  rht2_127@description = "Relative Humidity"
  rht2_127@units       = "%"
  rht2_127!0           = "Time"
  rht2_127!1           = "bottom_top"
  rht2_127!2           = "horiz_dim_layer"

; SCM temperature to potential temperature
  thtctl_K = tmpctl_K
  thtctl_K = tmpctl_K * ((P00 / presctl_Pa)^(R_dry / c_p))

  thtt1_K = tmpt1_K
  thtt1_K = tmpt1_K * ((P00 / prest1_Pa)^(R_dry / c_p))

  thtt2_K = tmpt2_K
  thtt2_K = tmpt2_K * ((P00 / prest2_Pa)^(R_dry / c_p))

; Convert specific humidities to mixing ratios for theta-l calculation
  rvctl_kgkg = qvctl_kgkg
  rvctl_kgkg = qvctl_kgkg / (1 - qvctl_kgkg)

  rcctl_kgkg = qcctl_kgkg
  rcctl_kgkg = qcctl_kgkg / (1 - qcctl_kgkg)

  rcctl_conv = qcctl_conv
  rcctl_conv = qcctl_conv / (1 - qcctl_conv)

  rictl_kgkg = qictl_kgkg
  rictl_kgkg = qictl_kgkg / (1 - qictl_kgkg)

  rvt1_kgkg = qvt1_kgkg
  rvt1_kgkg = qvt1_kgkg / (1 - qvt1_kgkg)

  rct1_kgkg = qct1_kgkg
  rct1_kgkg = qct1_kgkg / (1 - qct1_kgkg)

  rct1_conv = qct1_conv
  rct1_conv = qct1_conv / (1 - qct1_conv)

  rit1_kgkg = qit1_kgkg
  rit1_kgkg = qit1_kgkg / (1 - qit1_kgkg)

  rvt2_kgkg = qvt2_kgkg
  rvt2_kgkg = qvt2_kgkg / (1 - qvt2_kgkg)

  rct2_kgkg = qct2_kgkg
  rct2_kgkg = qct2_kgkg / (1 - qct2_kgkg)

  rct2_conv = qct2_conv
  rct2_conv = qct2_conv / (1 - qct2_conv)

  rit2_kgkg = qit2_kgkg
  rit2_kgkg = qit2_kgkg / (1 - qit2_kgkg)

; Liquid water mixing ratio
  rlctl_kgkg = rcctl_kgkg
  rlctl_kgkg = rcctl_kgkg + rcctl_conv

  rlt1_kgkg = rct1_kgkg
  rlt1_kgkg = rct1_kgkg + rct1_conv

  rlt2_kgkg = rct2_kgkg
  rlt2_kgkg = rct2_kgkg + rct2_conv

; Total water mixing ratio (vapor+cloud+cloudconv+ice)
  rtotctl = rvctl_kgkg
  rtotctl = rvctl_kgkg + rlctl_kgkg + rictl_kgkg

  rtott1  = rvt1_kgkg
  rtott1  = rvt1_kgkg + rlt1_kgkg + rit1_kgkg

  rtott2  = rvt2_kgkg
  rtott2  = rvt2_kgkg + rlt2_kgkg + rit2_kgkg

; Calculate gamma and chi for liquid water potential temperature equation
  gammactl = (rtotctl * R_vap) / (c_p + (rtotctl * c_pv))
  gammat1  = (rtott1 * R_vap) / (c_p + (rtott1 * c_pv))
  gammat2  = (rtott2 * R_vap) / (c_p + (rtott2 * c_pv))

  chictl = new((/ dimsizes(qvctl_kgkg) /),float)
  chit1  = new((/ dimsizes(qvt1_kgkg) /),float)
  chit2  = new((/ dimsizes(qvt2_kgkg) /),float)
  do i=0,dimsizes(hscm)-1,1
    do j=0,dimsizes(presctl_Pa(0,:,0))-1,1
      if (rhctl_127(i,j,:) .lt. 100.) then
        chictl(i,j,:) = R_dry / c_p
      else
        chictl(i,j,:) = (R_dry / c_p) * ((1 + (rvctl_kgkg(i,j,:) / eps)) / \
                                        (1 + ((rvctl_kgkg(i,j,:)*c_pv) / \
                                        (c_p))))
      end if
      if (rht1_127(i,j,:) .lt. 100.) then
        chit1(i,j,:) = R_dry / c_p
      else
        chit1(i,j,:) = (R_dry / c_p) * ((1 + (rvt1_kgkg(i,j,:) / eps)) / \
                                        (1 + ((rvt1_kgkg(i,j,:)*c_pv) / \
                                        (c_p))))
      end if
      if (rht2_127(i,j,:) .lt. 100.) then
        chit2(i,j,:) = R_dry / c_p
      else
        chit2(i,j,:) = (R_dry / c_p) * ((1 + (rvt2_kgkg(i,j,:) / eps)) / \
                                        (1 + ((rvt2_kgkg(i,j,:)*c_pv) / \
                                        (c_p))))
      end if
    end do
  end do

; Calculate liquid water potential temperature for the SCM experiments
  thtlctl_K = thtctl_K *((((eps + rvctl_kgkg) / (eps + rtotctl))^chictl) * \
                        ((rvctl_kgkg / rtotctl)^(gammactl * -1)) * \
                        exp(((l_v * -1) * rlctl_kgkg) / ((c_p + (rtotctl * \
                        c_pv)) * tmpctl_K)))
  thtlt1_K  = thtt1_K *((((eps + rvt1_kgkg) / (eps + rtott1))^chit1) * \
                       ((rvt1_kgkg / rtott1)^(gammat1 * -1)) * \
                       exp(((l_v * -1) * rlt1_kgkg) / ((c_p + (rtott1 * \
                       c_pv)) * tmpt1_K)))
  thtlt2_K  = thtt2_K *((((eps + rvt2_kgkg) / (eps + rtott2))^chit2) * \
                       ((rvt2_kgkg / rtott2)^(gammat2 * -1)) * \
                       exp(((l_v * -1) * rlt2_kgkg) / ((c_p + (rtott2 * \
                       c_pv)) * tmpt2_K)))

; Averaging, starting from 0900 AEST (on day 33) until 0000 AEST the next day
; (day 34). 3 hour averages
  thtlctl_06 = thtlctl_K(closest_val(21600.,tscm):closest_val(32400.,tscm)-1, \
                         :,0)
  thtlt1_06  = thtlt1_K(closest_val(21600.,tscm):closest_val(32400.,tscm)-1, \
                        :,0)
  thtlt2_06  = thtlt2_K(closest_val(21600.,tscm):closest_val(32400.,tscm)-1, \
                        :,0)
  
  thtlctl_tavg06 = dim_avg_n_Wrap(thtlctl_06,0)
  thtlt1_tavg06  = dim_avg_n_Wrap(thtlt1_06,0)
  thtlt2_tavg06  = dim_avg_n_Wrap(thtlt2_06,0)

  thtlctl_09 = thtlctl_K(closest_val(32400.,tscm):closest_val(43200.,tscm)-1, \
                         :,0)
  thtlt1_09  = thtlt1_K(closest_val(32400.,tscm):closest_val(43200.,tscm)-1, \
                        :,0)
  thtlt2_09  = thtlt2_K(closest_val(32400.,tscm):closest_val(43200.,tscm)-1, \
                        :,0)

  thtlctl_tavg09 = dim_avg_n_Wrap(thtlctl_09,0)
  thtlt1_tavg09  = dim_avg_n_Wrap(thtlt1_09,0)
  thtlt2_tavg09  = dim_avg_n_Wrap(thtlt2_09,0)

  thtlctl_12 = thtlctl_K(closest_val(43200.,tscm):closest_val(54000.,tscm)-1, \
                         :,0)
  thtlt1_12  = thtlt1_K(closest_val(43200.,tscm):closest_val(54000.,tscm)-1, \
                        :,0)
  thtlt2_12  = thtlt2_K(closest_val(43200.,tscm):closest_val(54000.,tscm)-1, \
                        :,0)

  thtlctl_tavg12 = dim_avg_n_Wrap(thtlctl_12,0)
  thtlt1_tavg12  = dim_avg_n_Wrap(thtlt1_12,0)
  thtlt2_tavg12  = dim_avg_n_Wrap(thtlt2_12,0)

  thtlctl_15 = thtlctl_K(closest_val(54000.,tscm):closest_val(64800.,tscm)-1, \
                         :,0)
  thtlt1_15  = thtlt1_K(closest_val(54000.,tscm):closest_val(64800.,tscm)-1, \
                        :,0)
  thtlt2_15  = thtlt2_K(closest_val(54000.,tscm):closest_val(64800.,tscm)-1, \
                        :,0)
                        
  thtlctl_tavg15 = dim_avg_n_Wrap(thtlctl_15,0)
  thtlt1_tavg15  = dim_avg_n_Wrap(thtlt1_15,0)
  thtlt2_tavg15  = dim_avg_n_Wrap(thtlt2_15,0)

  thtlctl_18 = thtlctl_K(closest_val(64800.,tscm):closest_val(75600.,tscm)-1, \
                         :,0)
  thtlt1_18  = thtlt1_K(closest_val(64800.,tscm):closest_val(75600.,tscm)-1, \
                        :,0)
  thtlt2_18  = thtlt2_K(closest_val(64800.,tscm):closest_val(75600.,tscm)-1, \
                        :,0)
                        
  thtlctl_tavg18 = dim_avg_n_Wrap(thtlctl_18,0)
  thtlt1_tavg18  = dim_avg_n_Wrap(thtlt1_18,0)
  thtlt2_tavg18  = dim_avg_n_Wrap(thtlt2_18,0)

  thtlctl_21 = thtlctl_K(closest_val(75600.,tscm):closest_val(86400.,tscm)-1, \
                         :,0)
  thtlt1_21  = thtlt1_K(closest_val(75600.,tscm):closest_val(86400.,tscm)-1, \
                        :,0)
  thtlt2_21  = thtlt2_K(closest_val(75600.,tscm):closest_val(86400.,tscm)-1, \
                        :,0)
                        
  thtlctl_tavg21 = dim_avg_n_Wrap(thtlctl_21,0)
  thtlt1_tavg21  = dim_avg_n_Wrap(thtlt1_21,0)
  thtlt2_tavg21  = dim_avg_n_Wrap(thtlt2_21,0)

; End time averaging
; Height range for plotting
  h6k = closest_val(3000.,hgtctl_m(1,:))

; END VARIABLE input and calculations
; Open workstation for plotting
  wks = gsn_open_wks(gtype,"Wangara33FF_Theta-L")

; Resources
  res = True

  res@gsnMaximize      = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@xyLineThicknessF = 3.4
  res@gsnLeftString    = "(a) 0600-0900 AEST"
  res@gsnRightString   = "Wangara Day 33"
  res@tiXAxisString    = "Liquid Water Potential Temp. (K)"
  res@tiYAxisString    = "Height (m)"
  res@vpWidthF         = 0.99
  res@vpHeightF        = 0.5
  res@trXMinF          = 274.
  res@trXMaxF          = 312.
  res@trYMaxF          = 3000.

; Ranges for all 6 hour heights
  thtlctl_havg06 = thtlctl_tavg06(0:h6k(0))
  thtlt1_havg06  = thtlt1_tavg06(0:h6k(0))
  thtlt2_havg06  = thtlt2_tavg06(0:h6k(0))
; Define THE height range
  h_range    = hgtctl_m(1,0:h6k(0))

; Plot 06 hour EDMF
  pltctl = gsn_csm_xy(wks,thtlctl_havg06,h_range,res)

; Resources and plot for ACM
  res@xyLineColor      = "red"

  pltt1 = gsn_csm_xy(wks,thtlt1_havg06,h_range,res)

; Resources and plot for YSU
  res@xyLineColor      = "blue"

  pltt2 = gsn_csm_xy(wks,thtlt2_havg06,h_range,res)

; Overlay plots
  overlay(pltctl,pltt1)
  overlay(pltctl,pltt2)

; Resources for legend
  genres                    = True
  genres@XPosPercent        = 18.5
  genres@YPosPercent        = 71
  lineres                   = True
  lineres@lgLineThicknesses = 3.4
  lineres@lgLineColors      = (/ "black","red","blue" /)
  textres                   = True
  textres@lgLabels          = pbl_list

  draw(pltctl)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

; Resources for 09 hour plots
  res@gsnLeftString    = "(b) 0900-1200 AEST"
  res@xyLineColor      = "black"

; Ranges for all 9 hour heights
  thtlctl_havg09 = thtlctl_tavg09(0:h6k(0))
  thtlt1_havg09  = thtlt1_tavg09(0:h6k(0))
  thtlt2_havg09  = thtlt2_tavg09(0:h6k(0))

; Plot 09 hour EDMF
  pltctl = gsn_csm_xy(wks,thtlctl_havg09,h_range,res)

; Resources and plot for ACM
  res@xyLineColor      = "red"

  pltt1 = gsn_csm_xy(wks,thtlt1_havg09,h_range,res)

; Resources and plot for YSU
  res@xyLineColor      = "blue"

  pltt2 = gsn_csm_xy(wks,thtlt2_havg09,h_range,res)

; Overlay plots
  overlay(pltctl,pltt1)
  overlay(pltctl,pltt2)

  draw(pltctl)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

; Resources for 12 hour plots
  res@gsnLeftString    = "(c) 1200-1500 AEST"
  res@xyLineColor      = "black"

; Ranges for all 12 hour heights
  thtlctl_havg12 = thtlctl_tavg12(0:h6k(0))
  thtlt1_havg12  = thtlt1_tavg12(0:h6k(0))
  thtlt2_havg12  = thtlt2_tavg12(0:h6k(0))

; Plot 12 hour EDMF
  pltctl = gsn_csm_xy(wks,thtlctl_havg12,h_range,res)

; Resources and plot for ACM
  res@xyLineColor      = "red"

  pltt1 = gsn_csm_xy(wks,thtlt1_havg12,h_range,res)

; Resources and plot for YSU
  res@xyLineColor      = "blue"

  pltt2 = gsn_csm_xy(wks,thtlt2_havg12,h_range,res)

; Overlay plots
  overlay(pltctl,pltt1)
  overlay(pltctl,pltt2)

  draw(pltctl)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

; Resources for 15 hour plots
  res@gsnLeftString    = "(d) 1500-1800 AEST"
  res@xyLineColor      = "black"

; Ranges for all 15 hour heights
  thtlctl_havg15 = thtlctl_tavg15(0:h6k(0))
  thtlt1_havg15  = thtlt1_tavg15(0:h6k(0))
  thtlt2_havg15  = thtlt2_tavg15(0:h6k(0))

; Plot 15 hour EDMF
  pltctl = gsn_csm_xy(wks,thtlctl_havg15,h_range,res)

; Resources and plot for ACM
  res@xyLineColor      = "red"

  pltt1 = gsn_csm_xy(wks,thtlt1_havg15,h_range,res)

; Resources and plot for YSU
  res@xyLineColor      = "blue"
  
  pltt2 = gsn_csm_xy(wks,thtlt2_havg15,h_range,res)

; Overlay plots
  overlay(pltctl,pltt1)
  overlay(pltctl,pltt2)

  draw(pltctl)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

; Resources for 18 hour plots
  res@gsnLeftString    = "(e) 1800-2100 AEST"
  res@xyLineColor      = "black"

; Ranges for all 18 hour heights
  thtlctl_havg18 = thtlctl_tavg18(0:h6k(0))
  thtlt1_havg18  = thtlt1_tavg18(0:h6k(0))
  thtlt2_havg18  = thtlt2_tavg18(0:h6k(0))

; Plot 18 hour EDMF
  pltctl = gsn_csm_xy(wks,thtlctl_havg18,h_range,res)

; Resources and plot for ACM
  res@xyLineColor      = "red"

  pltt1 = gsn_csm_xy(wks,thtlt1_havg18,h_range,res)

; Resources and plot for YSU
  res@xyLineColor      = "blue"

  pltt2 = gsn_csm_xy(wks,thtlt2_havg18,h_range,res)

; Overlay plots
  overlay(pltctl,pltt1)
  overlay(pltctl,pltt2)

  draw(pltctl)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

; Resources for 21 hour plots
  res@gsnLeftString    = "(f) 2100-0000 AEST"
  res@gsnRightString   = "Wangara Day 33+34"
  res@xyLineColor      = "black"
  
; Ranges for all 21 hour heights
  thtlctl_havg21 = thtlctl_tavg21(0:h6k(0))
  thtlt1_havg21  = thtlt1_tavg21(0:h6k(0))
  thtlt2_havg21  = thtlt2_tavg21(0:h6k(0))

; Plot 18 hour EDMF
  pltctl = gsn_csm_xy(wks,thtlctl_havg21,h_range,res)

; Resources and plot for ACM
  res@xyLineColor      = "red"

  pltt1 = gsn_csm_xy(wks,thtlt1_havg21,h_range,res)

; Resources and plot for YSU
  res@xyLineColor      = "blue"

  pltt2 = gsn_csm_xy(wks,thtlt2_havg21,h_range,res)

; Overlay plots
  overlay(pltctl,pltt1)
  overlay(pltctl,pltt2)

  draw(pltctl)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

  print("Complete: "+systemfunc("date"))

end

