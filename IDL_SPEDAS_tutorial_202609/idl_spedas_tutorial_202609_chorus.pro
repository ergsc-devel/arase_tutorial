;----------------------------------------------------------------------------------
pro idl_spedas_tutorial_202609_chorus, no_load_data=no_load_data
;----------------------------------------------------------------------------------
;
;Description
;
; --- This procedure calculates the resonant as a function of pitch angle
;     for the chrous waves shown in Kurita et al. (2018).
;
;     Kurita, S., Miyoshi, Y., Kasahara, S.,Yokota, S., Kasahara, Y., Matsuda, S.,
;     et al.(2018). Deformation of electron pitchangle distributions caused by
;     upperband chorus observed by the Arasesatellite. Geophysical Research Letters,
;     45, 7996–8004. https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018GL079104
;
;Example
;
; --- IDL> idl_spedas_tutorial_202609_chorus
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
; --- 2026/Sep/12 :[Drafted by KY] 
;    
;----------------------------------------------------------------------------------
;Setting

  ; Compile option
  compile_opt idl2

  ; Observational parameters reported in Kurita et al. (2018)
  fpfc = 3.30 ; fpe/fce = 3.3
  xl   = 0.50 ; lower cutoff frequency
  xu   = 0.66 ; upper cutoff frequency

;----------------------------------------------------------------------------------
;Calculation

  ;
  ; --- (1) Calculate resonance energy
  ;

  ; Pitch angle 
  pa = findgen(start=1,179) ; in degree

  ; --- Lower cutoff frequency

  ; Refractive index
  n = sqrt(1 + fpfc^2/(xl*(1-xl)))

  ; Beta
  ; kpara < 0
  bl_1 = (    n*xl^2* cosd(pa) + sqrt( n^2*xl^2*cosd(pa)^2 - xl^2 + 1) ) $
       / ( n^2*xl^2*cosd(pa)^2 + 1 )  
  ; kpara > 0
  bl_2 = ( -1*n*xl^2* cosd(pa) + sqrt( n^2*xl^2*cosd(pa)^2 - xl^2 + 1) ) $
       / ( n^2*xl^2*cosd(pa)^2 + 1 )  

  ; remove beta < 0 or beta > 1 or 1-n*beta*cosPA > 0 (for kpara > 0)
  sub_1 = where(bl_1 lt 0 or bl_1 gt 1 or 1-n*bl_1*cosd(pa) lt 0, cnt_1)  
  ; remove beta < 0 or beta > 1 or 1+n*beta*cosPA > 0 (for kpara < 0)
  sub_2 = where(bl_2 lt 0 or bl_2 gt 1 or 1+n*bl_2*cosd(pa) lt 0, cnt_2)  
  if cnt_1 gt 0 then bl_1[sub_1] = !values.f_nan
  if cnt_2 gt 0 then bl_2[sub_2] = !values.f_nan

  ; Gamma
  gl_1 = 1/sqrt(1-bl_1^2)
  gl_2 = 1/sqrt(1-bl_2^2)

  ; Resonance energy for each pitch angle
  El_1 = 511*(gl_1 - 1) ; in KeV
  El_2 = 511*(gl_2 - 1) ; in KeV
  
  ; --- Upper cutoff frequency

  ; Refractive index
  n = sqrt(1 + fpfc^2/(xu*(1-xu)))

  ; Beta
  ; kpara < 0
  bu_1 = (    n*xu^2* cosd(pa) + sqrt( n^2*xu^2*cosd(pa)^2 - xu^2 + 1) ) $
       / ( n^2*xu^2*cosd(pa)^2 + 1 )  
  ; kpara > 0
  bu_2 = ( -1*n*xu^2* cosd(pa) + sqrt( n^2*xu^2*cosd(pa)^2 - xu^2 + 1) ) $
       / ( n^2*xu^2*cosd(pa)^2 + 1 )  
 
  ; remove beta < 0 or beta > 1 or 1-n*beta*cosPA > 0 (for kpara > 0)
  sub_1 = where(bu_1 lt 0 or bu_1 gt 1 or 1-n*bu_1*cosd(pa) lt 0, cnt_1)  
  ; remove beta < 0 or beta > 1 or 1+n*beta*cosPA > 0 (for kpara < 0)
  sub_2 = where(bu_2 lt 0 or bu_2 gt 1 or 1+n*bu_2*cosd(pa) lt 0, cnt_2)  
  if cnt_1 gt 0 then bu_1[sub_1] = !values.f_nan
  if cnt_2 gt 0 then bu_2[sub_2] = !values.f_nan

  ; remove beta < 0 or beta > 1
  sub_1 = where(bu_1 lt 0 or bu_1 gt 1, cnt_1)  
  sub_2 = where(bu_2 lt 0 or bu_2 gt 1, cnt_2)  
  if cnt_1 gt 0 then bu_1 = !values.f_nan
  if cnt_2 gt 0 then bu_2 = !values.f_nan

  ; Gamma
  gu_1 = 1/sqrt(1-bu_1^2)
  gu_2 = 1/sqrt(1-bu_2^2)

  ; Resonance energy for each pitch angue
  Eu_1 = 511*(gu_1 - 1) ; in KeV
  Eu_2 = 511*(gu_2 - 1) ; in KeV


  ;
  ; --- (2) Calculate pitch angle distribution
  ;

  timespan, '20170408'

  ; Load MEP-e data
  if keyword_set(no_load_data) then goto, skip_load_data
  erg_load_mepe, level='l3', datatype='3dflux'
  skip_load_data:

  ; Get IDL structure
  get_data, 'erg_mepe_l3_3dflux_FEDU'      , data=flux
  get_data, 'erg_mepe_l3_3dflux_FEDU_alpha', data=pa_3d

  ; parameters
  nt  = n_elements(flux.y[*,0,0,0])
  nsp = n_elements(flux.y[0,*,0,0])
  nen = n_elements(flux.y[0,0,*,0])
  nch = n_elements(flux.y[0,0,0,*])
  enbin = flux.v2
  ; time for panel c
  time = '2017-04-08/19:20:33' 
  ; index for time
  idx  = nn(flux.x,time)

  ; Setting for binning of pitch angle distribution
  dpa      = 5.0
  flux_pad = fltarr(180.0/dpa,nen)

  ; loop for energy
  for ien = 0, nen-1 do begin

    ; Create 1D data
    flux_1d = reform(flux.y[idx,*,ien,*],[nsp*nch])
    pa_1d   = reform(pa_3d.y[idx,*,ien,*],[nsp*nch])

    ; Binning with a 5 deg pitch angle bin
    bin1d, pa_1d, flux_1d, 0.0, 180.0, dpa, binnm, pabin, ave_flux
    flux_pad[*,ien] = ave_flux

  endfor


  ;
  ; --- (3) Plot pitch angle distribution and resonance energy
  ;

  ; Set color table
  loadct_sd, 48

  ; Plot color amp
  plotxyz, pabin, enbin, flux_pad $
         , title=time $
         , xrange=[0,180], xstyle=1, xtitle='Local Pitch Angle [deg.]' $
         , yrange=[6.5,35], ystyle=1, ylog=1, ytitle='Electron Energy [keV/q]' $
         , zlog=1, zrange=[5.e4,2.e7], ztitle='[s!u-1!ncm!u-2!nkeV!u-1!nsr!u-1!n]' $
         , xtickinterval=30, /noisotropic, zticklen=-0.4, xmargin=[0.17,0.20], ymargin=[0.17,0.10]

  ; Plot resonance energy for parallel propagation
  oplot, pa, Eu_1, thick=3
  oplot, pa, El_1, lines=2, thick=3

  ; Plot resonance energy for anti-parallel propagation
  oplot, pa, Eu_2, thick=3
  oplot, pa, El_2, lines=2, thick=3

;----------------------------------------------------------------------------------
;Output
;----------------------------------------------------------------------------------

stop
end

