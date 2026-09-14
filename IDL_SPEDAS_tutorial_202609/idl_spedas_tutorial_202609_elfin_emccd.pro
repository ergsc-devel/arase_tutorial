;----------------------------------------------------------------------------------
pro idl_spedas_tutorial_202609_elfin_emccd, no_load_data=no_load_data
;----------------------------------------------------------------------------------
;
;Description
;
; --- This procedure demonstrates a conjunction observation between the PsA/PWING
;     EMCCD all-sky camera and the ELFIN-A satellite on September 28, 2019.
;
;Example
;
; --- IDL> idl_spedas_tutorial_202609_elfin_emccd
;
;Input
;
; --- N/A
;
;Keyword
;
; --- no_load_data :[0 or 1]
;                   If set, the procedure of data loading is skipped.
;
;Output
;
; --- N/A
;
;History
;
; --- 2026/Sep/13 :[Drafted by KY]
;
;----------------------------------------------------------------------------------
;Setting

; Compile option
compile_opt idl2 

;----------------------------------------------------------------------------------
;Calculation


  ;
  ; --- (1) Load auroral image at Tjautjas
  ;
 
  ; Set time range ( A interval of <5 min is recommended.)
  timespan, '2019-09-28/00:43:00', 2, /min
  
  ; Option for data loading
  if keyword_set(no_load_data) then goto, skip_load_data

  ; Load 10-Hz EMCCD All sky image at Tjautjas with a table for mapping to geographic coordinates
  erg_load_camera_emccd_asf, site='tja', /mapping_table


  ;
  ; --- Mapping to geographic coordinates at a 120 km altitude
  ;

  ; Show message
  print, ' Convert pixel to geographic coordinates'

  ; Keywords of TASF2GMAP.PRO
  ; --- altitude          : mapping altitude [km (from 80 to 300)]
  ; --- grid_x            : longitude resolution [deg]
  ; --- grid_y            : latitude resolution [deg]
  ; --- fill_empty        : fill data gaps
  ; --- max_fill_distance : maximal number of the gap inteval to be filled
  tasf2gmap, 'emccd_asf_tja_image_raw' $
           , 'emccd_asf_tja_mapping_table' $
           , altitude=120, grid_x=0.05, grid_y=0.05, fill_empty=1, max_fill_distance=3

stop

  ;
  ; --- (2) Load ELFIN data
  ;

  ; load orbit data
  elf_load_state, probe='a'
  ; Calculate geophysical latitude, longitude, and altitude
  spd_cotrans, 'ela_pos_gei', 'ela_pos_geo', in_coord='gei', out_coord='geo'
  ttrace2iono, 'ela_pos_geo' $
             , in_coord='geo', out_coord='geo', newname='ela_pos_geo_alt120' $
             , internal_model='igrf', external_model='none' $
             , R0=120+6371.2, /standard_mapping, /km
  xyz_to_polar, 'ela_pos_geo_alt120', /ph_0_360 ; spherical coordinates

  ; load electron data
  elf_load_epd, probe='a', level='l1', cdf_version='v01', trange=timerange(/current)

stop 

  skip_load_data:

  ;
  ; --- (3) Plot multiple panels of all sky image
  ;

  get_data, 'emccd_asf_tja_image_raw', data=raw

  ; Load color table
  initct, 0, line_clr=8

  ; Plot setting
  idx       = nn(raw,'2019-09-28/00:43:35') ; plot data from 00:43:35 UT
  nstep     = 100 ; 10 Hz data --> 10 s interval
  xgrid     = 115 ; pixel of horizontal grid on the all sky images and Keogram
  zrange    = [2000,2600] ; range of raw count
  linestyle = 2

  ;==========================================================================
  ; Uncomment this part if you would like to create figures of all sky images
  ;==========================================================================
  ;window, xs=700, ys=600
  ;
  ;; for check
  ;for i = 0, 1800-1, 5 do begin
  ;  plotxyz, findgen(256), findgen(256), reform(raw.y[i,*,*]) $
  ;         , zrange=zrange $
  ;         , title=time_string(raw.x[i],tformat='YYYY MTH DD/hh:mm:ss.ff') $
  ;         , xtitle='X [pixel]', ytitle='Y [pixel]', ztitle='[Raw Count]'
  ;  makepng, 'erg_srvy_swg202609_aurora_ask_'+string(i,format='(i4.4)')
  ;endfor
  ;stop

  ; Set plot device
  window, 0, xs=1000, ys=600
  ; Plot all sky image in the first panel
  plotxyz, findgen(256), findgen(256), reform(raw.y[idx,*,*]) $
         , zrange=zrange $
         , mtitle='PsA EMCCD Camera (Tjautjas)' $
         , title=time_string(raw.x[idx],tformat='YYYY MTH DD/hh:mm:ss.ff') $
         , xtitle='X [pixel]', ytitle='Y [pixel]', ztitle='[Raw Count]', multi='3,2'
  ; Add grid
  oplot, [  0,256], [xgrid,xgrid], linestyle=linestyle
  oplot, [128,128], [  0,256], linestyle=linestyle

  ; Loop for time
  for i = 1, 5 do begin
 
    ; Plot all sky image in the other panels
    plotxyz, findgen(256), findgen(256), reform(raw.y[idx+i*nstep,*,*]) $
           , zrange=zrange $
           , title=time_string(raw.x[idx+i*nstep],tformat='YYYY MTH DD/hh:mm:ss.ff') $
           , xtitle='X [pixel]', ytitle='Y [pixel]', ztitle='[Raw Count]', addpanel=1
    ; Add grid
    oplot, [  0,256], [xgrid,xgrid], linestyle=linestyle  
    oplot, [128,128], [  0,256], linestyle=linestyle

  endfor

stop

  ;
  ; --- (4) Create Keogram in pixel
  ;

  x = raw.x
  y = reform(raw.y[*,128,*])
  v = findgen(256)

  store_data, 'tja_keogram_pix_x128' $
            , data={x:x,y:y,v:v}, lim={spec:1, no_interp:1, ystyle:1}

  ; Set a window for the keogram
  window, 1, xs=800, ys=1000

  options, 'tja_keogram_pix_x128', constant=xgrid

  tplot, 'tja_keogram_pix_x128', window=1
  timebar, raw.x[idx+indgen(6)*nstep], linestyle=linestyle

stop

  ;
  ; --- (5) Plot the footprint of ELFIN on auroral image in geographic coordinates
  ;

  get_data, 'emccd_asf_tja_image_raw_gmap_120', data=raw_gmap

  window, 2, xs=800, ys=600

  ; --- First panel

  idx1 =  nn(raw,'2019-09-28/00:44:19') ; index of all sky camera data for 00:44:19 UT

  plotxyz, raw_gmap.glon, raw_gmap.glat, reform(raw_gmap.y[idx1,*,*]) $
         , xrange=[15,25], yrange=[65,69], zrange=zrange $
         , noisotropic=1 $
         , title=time_string(raw.x[idx1],tformat='YYYY MTH DD/hh:mm:ss.ff') $
         , xtitle='Longitude [deg.]', ytitle='Latitude [deg.]', ztitle='[Raw Count]', multi='2,1'
 
  get_data, 'ela_pos_geo_alt120_phi' , data=ela_glon ; geographic longitude of ELFIN
  get_data, 'ela_pos_geo_alt120_th'  , data=ela_glat ; geographic latitude of ELFIN

  oplot, ela_glon.y, ela_glat.y
  idx_ela = nn(ela_glon,raw.x[idx1]) ; index of ELFIN data for 00:44:19 UT 
  oplot, [ela_glon.y[idx_ela]], [ela_glat.y[idx_ela]], psym=1, symsize=3, thick=3

  ; --- Second panel

  idx2 = nn(raw,'2019-09-28/00:44:28') ; index of all sky camera data for 00:44:28 UT

  plotxyz, raw_gmap.glon, raw_gmap.glat, reform(raw_gmap.y[idx2,*,*]) $
         , xrange=[15,25], yrange=[65,69], zrange=zrange $
         , noisotropic=1 $
         , title=time_string(raw.x[idx2],tformat='YYYY MTH DD/hh:mm:ss.ff') $
         , xtitle='Longitude [deg.]', ytitle='Latitude [deg.]', ztitle='[Raw Count]', addpanel=1
 
  oplot, ela_glon.y, ela_glat.y
  idx_ela = nn(ela_glon,raw.x[idx2]) ; index of ELFIN data for 00:44:28 UT
  oplot, [ela_glon.y[idx_ela]], [ela_glat.y[idx_ela]], psym=1, symsize=3, thick=3

stop

  ;
  ; --- (6) Create Keogram in the geophysical coordinates
  ;

  x   = raw_gmap.x
  tmp = min(abs(raw_gmap.glon-21.0),ilon)
  y   = reform(raw_gmap.y[*,ilon,*])
  v   = raw_gmap.glat

  store_data, 'tja_keogram_glon21' $
            , data={x:x,y:y,v:v}, lim={spec:1, no_interp:1, ystyle:1}

  options, 'tja_keogram_glon21', ytitle='TJA!CGLON = 21 deg.!CGLAT [deg.]'

  ; Create multi tplot variable
  store_data, 'tja_keogram_glon21_multi', data=['tja_keogram_glon21', 'ela_pos_geo_alt120_th']
  options, 'tja_keogram_pix_x128', yrange=[30,180]
  options, 'ela_pos_geo_alt120_th', thick=2
  options, 'tja_keogram_glon21_multi', yrange=[65,68], ystyle=1, constant=[67.0,67.55]
  options, 'ela_pef_nflux', ylog=1, zlog=1, color_table=1080, ystyle=1, ysubtitle='[keV]', ztitle='[/s/cm!u2!n/MeV/sr]'

  tplot, ['ela_pef_nflux','tja_keogram_glon21_multi'], add=1
  timebar, raw_gmap.x[[idx1,idx2]]

stop

  ;
  ; --- (7) Plot ERG data
  ;

  ; Load Arase data
  erg_load_pwe_ofa
  erg_load_mgf
  erg_load_orb_l3, model='ts04'
  split_vec, 'erg_orb_l3_pos_iono_north_TS04' 

  ; Set a window for ERG data
  window, 3, xs=800, ys=600

  ; Load a different color table
  loadct_sd, 48
  options, 'erg_orb_l3_pos_iono_north_TS04_0', yrange=[0,0], ystyle=2

  timespan, '2019-09-28/00:00:00', 6, /hour

  ; Plot chrous waves and L-shell
  tplot, window=3 $
       , ['erg_pwe_ofa_l2_spec_E_spectra_merged' $
         ,'erg_pwe_ofa_l2_spec_B_spectra_merged' $
         ,'erg_orb_l3_pos_iono_north_TS04_0']

  ; Add a vertical bar at the timing of the first panel (idx1)
  timebar, raw_gmap.x[idx1]

;----------------------------------------------------------------------------------
;Output

;----------------------------------------------------------------------------------

stop
end
