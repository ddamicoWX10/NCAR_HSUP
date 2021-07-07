load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; User modification area
  gtype = "x11"

  top_dir = "/glade/u/home/damico/gmtb-scm-GABLS3/scm/bin/output_gabls3_" + \
            "SCM_"

  gfs16 = "GFS_v16beta"
  gfs15 = "GFS_v15p2"

  date  = "2006-07-01 12:00:00"

  fileall = "output.nc"

; End user modification area
; Constants
  R_dry = 287.    ; ideal gas constant for dry air (J kg^-1 K^-1)
  R_vap = 461.5   ; gas constant for water vapor (J kg^-1 K^-1)
  c_p   = 1004.   ; specific heat at const pressure (J kg^-1 K^-1)
  l_v   = 2.5e6   ; latent heat of vaporization (J kg^-1)
  g     = 9.80665 ; acceleration due tp gravity (m s^-2)
  P00   = 100000. ; initial pressure (Pa)
  pi    = 4.*atan(1.)

; Add files
  g16 = addfile(top_dir+gfs16+"/"+fileall,"r")
  g15 = addfile(top_dir+gfs15+"/"+fileall,"r")

; Time
  time_s = g16->time

  time       = time_s
  time       = time_s / 60. / 60.
  time@units = "hours since "+date

; Add surface and soil temperature variables
  T_sfc_16   = g16->T_s
  T_sfc_15   = g15->T_s
  T_soil1_16 = g16->T_soil_1
  T_soil1_15 = g15->T_soil_1
  T_soil2_16 = g16->T_soil_2
  T_soil2_15 = g15->T_soil_2
  T_soil3_16 = g16->T_soil_3
  T_soil3_15 = g15->T_soil_3
  T_soil4_16 = g16->T_soil_4
  T_soil4_15 = g15->T_soil_4

; Add soil moisture variables
  smois1_16 = g16->soil_moisture_1
  smois1_15 = g15->soil_moisture_1
  smois2_16 = g16->soil_moisture_2
  smois2_15 = g15->soil_moisture_2
  smois3_16 = g16->soil_moisture_3
  smois3_15 = g15->soil_moisture_3
  smois4_16 = g16->soil_moisture_4
  smois4_15 = g15->soil_moisture_4

; Begin plotting
  wks = gsn_open_wks(gtype,"GABLS_SfcSoil")

; Resources for soil temperature
  res = True

  res@gsnMaximize      = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@xyLineThicknessF = 3.4
  res@gsnLeftString    = "Sfc and Soil Temp"
  res@gsnRightString   = "K"
  res@tiXAxisString    = time@units + " UTC"
  res@tiYAxisString    = "Temperature ("+T_sfc_16@units+")"
  res@trYMinF          = 282.
  res@trYMaxF          = 306.

  sfc16 = gsn_csm_xy(wks,time,T_sfc_16(:,0),res)

  res@xyDashPattern    = 1

  sfc15 = gsn_csm_xy(wks,time,T_sfc_15(:,0),res)

; Noah soil layer 1
  res@xyDashPattern    = 0
  res@xyLineColor      = "darkorange1"

  s116 = gsn_csm_xy(wks,time,T_soil1_16(:,0),res)

  res@xyDashPattern    = 1

  s115 = gsn_csm_xy(wks,time,T_soil1_15(:,0),res)

; Noah soil layer 2
  res@xyDashPattern    = 0
  res@xyLineColor      = "purple4"

  s216 = gsn_csm_xy(wks,time,T_soil2_16(:,0),res)

  res@xyDashPattern    = 1

  s215 = gsn_csm_xy(wks,time,T_soil2_15(:,0),res)

; Noah soil layer 3
  res@xyDashPattern    = 0
  res@xyLineColor      = "dodgerblue"
  
  s316 = gsn_csm_xy(wks,time,T_soil3_16(:,0),res)
  
  res@xyDashPattern    = 1
  
  s315 = gsn_csm_xy(wks,time,T_soil3_15(:,0),res)

; Noah soil layer 4
  res@xyDashPattern    = 0
  res@xyLineColor      = "deeppink"
  
  s416 = gsn_csm_xy(wks,time,T_soil4_16(:,0),res)
  
  res@xyDashPattern    = 1
  
  s415 = gsn_csm_xy(wks,time,T_soil4_15(:,0),res)

; Legend resources
  genres                    = True
  genres@XPosPercent        = 43
  genres@YPosPercent        = 90
  lineres                   = True
  lineres@lgLineColors      = (/ "black","black",             \
                                 "darkorange1","darkorange1", \
                                 "purple4","purple4",         \
                                 "dodgerblue","dodgerblue",   \
                                 "deeppink","deeppink" /)
  lineres@lgDashIndexes     = (/ 0,1,0,1,0,1,0,1,0,1 /)
  lineres@lgLineThicknesses = 3.4
  textres                   = True
  textres@lgLabels          = (/ gfs16+" Sfc",gfs15+" Sfc",   \
                                 gfs16+" ~F33~D~F21~Z 10cm",  \
                                 gfs15+" ~F33~D~F21~Z 10cm",  \
                                 gfs16+" ~F33~D~F21~Z 30cm",  \                 
                                 gfs15+" ~F33~D~F21~Z 30cm",  \
                                 gfs16+" ~F33~D~F21~Z 60cm",  \                 
                                 gfs15+" ~F33~D~F21~Z 60cm",  \
                                 gfs16+" ~F33~D~F21~Z 100cm", \                 
                                 gfs15+" ~F33~D~F21~Z 100cm" /)

; Overlay plots
  overlay(sfc16,sfc15)
  overlay(sfc16,s116)
  overlay(sfc16,s115)
  overlay(sfc16,s216)
  overlay(sfc16,s215)
  overlay(sfc16,s316)
  overlay(sfc16,s315)
  overlay(sfc16,s416)
  overlay(sfc16,s415)

; Draw plots and legend, frame workstation
  draw(sfc16)
;  drawNDCGrid(wks)
  simple_legend_ndc(wks,genres,lineres,textres)
  frame(wks)

; Resources for soil moisture
  res@gsnLeftString    = "Soil Moisture"
  res@gsnRightString   = "m~S~3~N~ m~S~-3~N~"
  res@tiYAxisString    = "Soil Moisture (m~S~3~N~ m~S~-3~N~)"
  res@trYMinF          = 0.186
  res@trYMaxF          = 0.196
  res@xyLineColor      = "darkorange1"
  res@xyDashPattern    = 0

; Soil level 1
  s116 = gsn_csm_xy(wks,time,smois1_16(:,0),res)

  res@xyDashPattern    = 1

  s115 = gsn_csm_xy(wks,time,smois1_15(:,0),res)

; Soil level 2
  res@xyLineColor      = "purple4"
  res@xyDashPattern    = 0

  s216 = gsn_csm_xy(wks,time,smois2_16(:,0),res)

  res@xyDashPattern    = 1

  s215 = gsn_csm_xy(wks,time,smois2_15(:,0),res)

; Soil level 3
  res@xyLineColor      = "dodgerblue"
  res@xyDashPattern    = 0
  
  s316 = gsn_csm_xy(wks,time,smois3_16(:,0),res)
  
  res@xyDashPattern    = 1

  s315 = gsn_csm_xy(wks,time,smois3_15(:,0),res)

; Soil level 4
  res@xyLineColor      = "deeppink"
  res@xyDashPattern    = 0
  
  s416 = gsn_csm_xy(wks,time,smois4_16(:,0),res)
  
  res@xyDashPattern    = 1

  s415 = gsn_csm_xy(wks,time,smois4_15(:,0),res)

; Overlay plots
  overlay(s116,s115)
  overlay(s116,s216)
  overlay(s116,s215)
  overlay(s116,s316)
  overlay(s116,s315)
  overlay(s116,s416)
  overlay(s116,s415)

;; Legend resources
;  genres@XPosPercent    = 50
;  genres@YPosPercent    = 50
;  lineres@lgLineColors  = (/ "darkorange1","darkorange1", \
;                             "purple4","purple4",         \
;                             "dodgerblue","dodgerblue",   \
;                             "deeppink","deeppink" /)
;  lineres@lgDashIndexes = (/ 0,1,0,1,0,1,0,1 /)
;  textres@lgLabels      = (/ gfs16+" ~F33~D~F21~Z 10cm",  \
;                             gfs15+" ~F33~D~F21~Z 10cm",  \
;                             gfs16+" ~F33~D~F21~Z 30cm",  \
;                             gfs15+" ~F33~D~F21~Z 30cm",  \
;                             gfs16+" ~F33~D~F21~Z 60cm",  \
;                             gfs15+" ~F33~D~F21~Z 60cm",  \
;                             gfs16+" ~F33~D~F21~Z 100cm", \
;                             gfs15+" ~F33~D~F21~Z 100cm" /)

; Draw plots and legend, frame workstation
;  draw(s116)
;  drawNDCGrid(wks)
;  simple_legend_ndc(wks,genres,lineres,textres)
;  frame(wks)

end

