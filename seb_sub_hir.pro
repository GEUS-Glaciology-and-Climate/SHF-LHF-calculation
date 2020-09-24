; Surface energy and mass budget model for ice sheets, by Dirk van As.
; The model can be run for a single location with a time series of p, T, RH, WS, SR, and LRin,
; or for a transect for which the input variables are height-dependent.
; The current version lacks:
; - retention by capillary forces
; - dry snow densification
; - reduction of sub-surface density by sub-surface melt
; -
; This is model version 2009/12, with plotting bits taken out and minute column introduced and year-2000 and fixed T_ice and snowthick ini.
; -
; Update November 2015, by Robert S. Fausto (RSF)
; - HIRHAM5 subsurface scheme has replaced the former one. The subsurface scheme now includes retention by capilary forces and dry snow densification. 
; - The subsurface subroutine is called "subsurface_hirham". See subroutine for description of parameters and variables. 
; - Layers mass (or water equivalent thickness) is constant in time. Each layer has a part which is snow (snowc), ice (snic) and water (slwc).
; - Total water equivalent thickness of layer n is thickness_weq(n) = snowc(n)+snic(n)+slwc(n). This number is constant in time.
; - Partitioning into snow, water and ice varies from time to time, and so does the density of the snow fraction (the variable rhofirn). 
; - This means that the actual thickness (as opposed to water eqiv), let’s call that ”act” (as opposed to “weq”), is:
;   thickness_act(n) = snowc(n)*(rho_w/rhofirn(n)) + snic*(rho_w/rho_ice) + slwc


pro SEB_sub_hir

  print,'Running...'

  ; AWS INPUT FILE INFO -----------------------------------------------------------------------
  stationname = 'QAS_L_hour_trimmed1'  ; filename without '.txt'
  directory = 'M:\Geus\work\e-drev\journal_papers\extreme_melt_south Greenland\SEB\'
  runID   = '_pr=10_SRcor_test_hirham'  ; addition to output file names
  dt_obs    = 3600.   ; observational time step in seconds
  rows    = 3816;15469;11100   ; number of rows of data in AWS input file QAS_L_hour_trimmed.txt
  ;rows    = 21325   ; number of rows of data in AWS input file KAN_L_hour_trimmed.txt
  ;rows   = 5000
  columns   = 44    ; number of columns in AWS input file

  col_year  = 1     ; Column designation: year
  col_day   = 5     ; day of year
  col_hour  = 4     ; time of day
  ;col_minutes  = 4     ; minute of hour
  col_pres  = 7     ; air pressure (hPa)
  col_T   = 8     ; air temperature (C)
  col_RH    = 11      ; relative humidity (%)
  col_WS    = 12      ; wind speed (m/s)
  col_SRin  = 15    ; incoming shortwave radiation (W/m2)
  col_SRout = 17    ; outgoing shortwave radiation (W/m2)
  col_LRin  = 19    ; incoming longwave radiation (W/m2)
  col_LRout = 20    ; outgoing longwave radiation (W/m2) - only for validation purposes
  col_Haws  = 23    ; height of sonic ranger on AWS (m) - only for validation purposes
  col_Hst   = 24    ; height of sonic ranger on stake assembly (m) - only for validation purposes
  col_Hpt   = 26    ; depth of pressure transducer in the ice (m) - only for validation purposes

  ; STATION&TRANSECT-SPECIFIC CONSTANTS FOR INITIALIZATION -----------------------------------------------------------------------
  ELA         = 1000.     ; equilibrium line altitude to approximate initial snow layer thickness in accumulation zone (m)
  elev_AWS      =  300      ; AWS elevation above sea level
  elev_bins     =   1     ; number of grid points along transect
  elev_binsize    =  100      ; vertical grid spacing (m)
  elev_start      =    0      ; lowest calculation level (m above sea level) -> if elev_bins = 1 then this will be set to elev_AWS
  gradT       =   -6.699e-3 ; -6.814e-3 ; along-slope vertical temperature gradient (C/m)
  gradT_Tdep      =    4.865e-5 ; temperature (in C) dependence of along-slope vertical temperature gradient (C/m/C)
  gradRH        =   17.504e-3 ; along-slope vertical relative humidity gradient (%/m)
  gradWS        =    3.145e-3 ; along-slope vertical wind speed gradient (m/s/m)
  gradLRin      =  -18.911e-3 ; -17.815e-3  ; along-slope vertical incoming longwave radiation gradient (W/m2/m)
  gradLRin_Tdep   =   -4.626e-4 ; temperature (in C) dependence of along-slope vertical incoming longwave radiation gradient (W/m2/m/C)
  gradsnowthick   =    0.1    ; along-slope vertical gradient in initial snow&firn layer thickness above ELA (m/m)
  gradTice      =   -0.001    ; along-slope vertical gradient in initial ice temperatures (C/m)
  H_T         =    2.6    ; initial height of temperature and humidity sensor over ice horizon (m)
  H_WS        =    3.1    ; initial height of anemometer over ice horizon (m)
  snowthick_ini   =    0.0      ; initial snow layer thickness below ELA (m)
  T_ice_AWS     =   -0.1    ; initial snow & ice temperature (C)

  ; PROGRAM CONSTANTS -----------------------------------------------------------------------
  dev         =   1   ; time step interpolation factor: 1,2,...,n (to keep sub-surface calculations stable)
  dTsurf_ini      =  10   ; initial surface temperature step in search for EB=0 (C)
  dz_ice        =   0.1    ; distance between vertical sub-surface levels (m)
  z_max         =   40.    ; maximum depth for sub-surface levels (m)
  EB_max        =   0.1   ; surface temperature iteration will stop when the imbalance in the energy budget is less then EB_max (W/m2)
  iter_max_EB     =  60;30    ; max number of iterations in reaching EB=0
  iter_max_flux   =  60;30    ; max number of iterations in calculating flux profiles
  L_dif       =   0.01  ; limit of relative change in Obukhov length during flux iterations
  prec_cutoff     = 0   ; elevation below which there is half precipitation
  prec_rate     =   1.e-3;13e-3 ; artificial precipitation rate (m of water equivalent per hour)
  RH_min        =  20.    ; lowest possible relative humidity over snow&ice after extrapolating over the DEM (%)
  WS_lim        =   0.8   ; wind speed below which turbulent fluxes are set to zero to prevent model instability (m/s)
  z_ice_max     = z_max/dz_ice   ; number of sub-surface levels

  ; CONSTANTS -----------------------------------------------------------------------
  alb_ice     =    0.45   ; albedo of ice -> if elev_bins = 1 then albedo will be determined from measurements
  alb_snow    =    0.8    ; albedo of snow -> if elev_bins = 1 then albedo will be determined from measurements
  eps       =    0.622    ;
  es_0      =    6.1071   ; saturation vapour pressure at the melting point (hPa)
  es_100      = 1013.246    ; saturation vapour pressure at steam point temperature (hPa)
  em        =    1.     ; longwave surface emissivity
  ext_air     =    1.5e-4   ; extinction coefficient of shortwave radiation in air for clear-sky conditions (m^-1)
  ext_air_cl    =    1.5e-4   ; cloud-dependent add-on to extinction coefficient of shortwave radiation in air (m^-1)
  ext_ice     =    1e10 ;1.4    ; extinction coefficient of shortwave radiation in ice (m^-1) (no penetration if set to 1e10)
  ext_snow    =    1e10 ;20.      ; extinction coefficient of shortwave radiation in snow (m^-1) (no penetration if set to 1e10)
  g       =    9.82   ; gravitational constant (at sea level). 80 N -> 9.83, 60 N -> 9.82
  gamma     =   16.     ; flux profile correction (Paulson & Dyer)
  kappa     =    0.4    ; Von Karman constant (0.35-0.42)
  L_sub     =    2.83e6   ; latent heat of sublimation (J/kg)
  L_fus     =    3.34e5   ; latent heat of fusion/melting (J/kg)
  L_vap     =    2.50e6   ; latent heat of vaporization (J/kg)
  R_d       =  287.05   ; gas constant of dry air
  R_v       =  461.51   ; gas constant of water vapour
  rho_ice     =  900.     ; density of ice (kg/m3)
  ;rho_snow   =  500.     ; density of snow (kg/m3)
  rho_water   =  999.8395   ; density of water at the melting point (kg/m3)
  sigma     =    5.6704e-8  ; Stefan-Boltzmann's constant
  T_solidprecip =    0.     ; Near-surface temperature below which precipitation is solid (C)
  T_0       =  273.15   ; melting point temperature (K)
  T_100     =  373.15   ; steam point temperature (K)
  z0_ice      =    0.005    ; surface roughness length for ice (m)
  z0_snow     =    0.0001    ; surface roughness length for snow (m)
  aa        =    0.7    ; flux profile correction constants (Holtslag & De Bruin '88)
  bb        =    0.75
  cc        =    5.
  dd        =    0.35
  beta      = 2317.
  c_pd      = 1005.     ; specific heat of dry air (J/kg/K)
  c_w       = 4210.     ; specific heat of water at 0 C (J/kg/K) (Source www.engineeringtoolbox.com)
  ch1       = fltarr(4)
  ch2       = fltarr(4)
  ch3       = fltarr(4)
  cq1       = fltarr(4)
  cq2       = fltarr(4)
  cq3       = fltarr(4)
  ch1       = ([1.25,0.149,0.317])  ; used in calculating roughness length for heat z_h over smooth surfaces (Andreas 1987)
  ch2       = ([0,-0.55,-0.565])
  ch3       = ([0,0,-0.183])
  cq1       = ([1.61,0.351,0.396])  ; used in calculating roughness length for moisture z_q over smooth surfaces (Andreas 1987)
  cq2       = ([0,-0.628,-0.512])
  cq3       = ([0,0,-0.180])

  starttime = SYSTIME(1)

  ; DECLARING ARRAYS AND READING AWS DATA -----------------------------------------------------------------------
  line  = fltarr(columns)
  year  = fltarr(rows)
  day   = fltarr(rows)
  hour  = fltarr(rows)
  minutes = fltarr(rows)
  pres  = fltarr(elev_bins,rows)
  T   = fltarr(elev_bins,rows)
  RH    = fltarr(elev_bins,rows)
  q   = fltarr(elev_bins,rows)
  WS    = fltarr(elev_bins,rows)
  SRin  = fltarr(elev_bins,rows)
  SRout = fltarr(elev_bins,rows)
  LRin  = fltarr(elev_bins,rows)
  LRout = fltarr(rows)
  H_aws = fltarr(rows)
  H_st  = fltarr(rows)
  H_pt  = fltarr(rows)

  filename = directory+stationname+'.txt'
  header=''
  openr,1,filename
  readf,1,header
  for i=0L,rows-1 do begin
    if eof(1) then break
    readf,1,line
    year[i]   = line[col_year-1]-2000
    day[i]    = line[col_day-1]
    hour[i]   = line[col_hour-1]
    ;  minutes[i] = line[col_minutes-1]
    pres[0,i]   = line[col_pres-1]
    T[0,i]    = line[col_T-1]
    RH[0,i]   = line[col_RH-1]
    WS[0,i]   = line[col_WS-1]
    SRin[0,i]   = line[col_SRin-1]
    SRout[0,i]  = line[col_SRout-1]
    LRin[0,i]   = line[col_LRin-1]
    LRout[i]    = line[col_LRout-1]
    H_aws[i]    = line[col_Haws-1]
    H_st[i]   = line[col_Hst-1]
    H_pt[i]   = line[col_Hpt-1]
  endfor
  close,1

  ;minutes = hour - 100*fix(hour/100)
  ;hour = minutes/60. + hour
  time = year + (day-1 + hour/24.)/365.
  leapyear = where(year/4 eq fix(year/4))
  if total(leapyear) gt -1 then time[leapyear] = year[leapyear]+(day[leapyear]-1+hour[leapyear]/24.)/366
  ; correctly calculated leapyear
  ;leapyear = (((year MOD 4) EQ 0) AND ((year MOD 100) NE 0)) OR ((year MOD 400) EQ 0)
  ;if total(leapyear) gt 0 then time[leapyear] = year[leapyear]+(day[leapyear]-1+hour[leapyear]/24.)/366
  ;if total(leapyear) gt 0 then year_sec=3600.*24.*366. else year_sec= 3600.*24.*365.

  ; REPLACING MISSING AWS VALUES THROUGH INTERPOLATION -----------------------------------------------------------------------
  missingvalues = intarr(8)
  missing = where(pres[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(pres[0,*] ne -999)
    pres[0,missing] = interpol(pres[0,notmissing],notmissing,missing)
    missingvalues[0] = size(missing,/n_elements)
  endif
  missing = where(T[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(T[0,*] ne -999)
    T[0,missing] = interpol(T[0,notmissing],notmissing,missing)
    missingvalues[1] = size(missing,/n_elements)
  endif
  missing = where(RH[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(RH[0,*] ne -999)
    RH[0,missing] = interpol(RH[0,notmissing],notmissing,missing)
    missingvalues[2] = size(missing,/n_elements)
  endif
  missing = where(WS[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(WS[0,*] ne -999)
    WS[0,missing] = interpol(WS[0,notmissing],notmissing,missing)
    missingvalues[3] = size(missing,/n_elements)
  endif
  missing = where(SRin[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(SRin[0,*] ne -999)
    SRin[0,missing] = interpol(SRin[0,notmissing],notmissing,missing)
    missingvalues[4] = size(missing,/n_elements)
  endif
  missing = where(SRout[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(SRout[0,*] ne -999)
    SRout[0,missing] = interpol(SRout[0,notmissing],notmissing,missing)
    missingvalues[5] = size(missing,/n_elements)
  endif
  missing = where(LRin[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(LRin[0,*] ne -999)
    LRin[0,missing] = interpol(LRin[0,notmissing],notmissing,missing)
    missingvalues[6] = size(missing,/n_elements)
  endif
  missing = where(LRout[0,*] eq -999)
  if total(missing) gt -1 then begin
    notmissing = where(LRout[0,*] ne -999)
    LRout[0,missing] = interpol(LRout[0,notmissing],notmissing,missing)
    missingvalues[6] = size(missing,/n_elements)
  endif
  print,'Missing values for p,T,RH,WS,SRin,SRout,LRin,LRout: ',missingvalues

  T = T+T_0       ; Temperature in degrees Kelvin
  SRout_mdl = SRout   ; modelled instead of measured outgoing shortwave radiation, unless elev_bins = 1

  ; INTERPOLATION TO OBTAIN HIGHER TIME RESOLUTION -----------------------------------------------------------------------
  M = long(rows)*dev
  if dev ne 1 then begin
    ipfac1 = M+1
    time        = interpol(time,ipfac1)
    pres[0,*]     = interpol(pres[0,*],ipfac1)
    T[0,*]      = interpol(T[0,*],ipfac1)
    RH[0,*]     = interpol(RH[0,*],ipfac1)
    WS[0,*]     = interpol(WS[0,*],ipfac1)
    SRin[0,*]     = interpol(SRin[0,*],ipfac1)
    SRout_mdl[0,*]  = interpol(SRout_mdl[0,*],ipfac1)
    LRin[0,*]     = interpol(LRin[0,*],ipfac1)
  endif
  T_AWS = T           ; needed for precipitation estimate
  LRin_AWS = LRin
  for i=0L,M-1 do begin
    T_AWS[*,i] = T_AWS[0,i]
    LRin_AWS[*,i] = LRin_AWS[0,i]
  endfor

  ; EXTRAPOLATING AWS DATA OVER THE TRANSECT -----------------------------------------------------------------------
  if elev_bins eq 1 then elev_start = elev_AWS
  elev = findgen(elev_bins)*elev_binsize + elev_start
  for i=0l,M-1 do begin   ; if a for-loop over all time steps is too slow, perhaps use the rebin function
    pres[*,i] = pres[0,i]*(T[0,i]/(T[0,i]+(gradT+gradT_Tdep*(T[0,i]-T_0))*(elev-elev_AWS)))^(g/R_d/(gradT+gradT_Tdep*(T[0,i]-T_0)))
    T[*,i]  = T[0,i] + (elev-elev_AWS)*(gradT + gradT_Tdep*(T[0,i]-T_0))
    RH[*,i] = RH[0,i] + (elev-elev_AWS)*gradRH
    WS[*,i] = WS[0,i] + (elev-elev_AWS)*gradWS
    SRin[*,i] = SRin[0,i]*exp(-ext_air*(elev_AWS-elev))   ; this needs inclusion of cloud effect and air mass dependence
    LRin[*,i] = LRin[0,i] + (elev-elev_AWS)*(gradLRin + gradLRin_Tdep*(T[0,i]-T_0))
  endfor
  RHtoolow = where(RH lt RH_min)
  if total(RHtoolow) gt -1 then RH[RHtoolow] = RH_min   ; setting a lower limit for relative humidity (for supersaturation: see below)
  WStoolow = where(WS lt 0)
  if total(WStoolow) gt -1 then WS[WStoolow] = 0      ; setting a lower limit for wind speed

  rho_atm = 100*pres/R_d/T                ; atmospheric density
  mu = 18.27e-6*(291.15+120)/(T+120)*(T/291.15)^1.5   ; dynamic viscosity of air (Pa s) (Sutherlands' formula using C = 120 K)
  nu = mu/rho_atm                     ; kinematic viscosity of air (m^2/s)
  sum_dH_surf=0
  ; DECLARATION OF ARRAYS OF MODELLED VARIABLES -----------------------------------------------------------------------
  ci =fltarr(z_ice_max+1)
  ki = fltarr(z_ice_max+1)
  D = fltarr(z_ice_max+1)
  l = fltarr(z_ice_max+1)
  alpha_t = fltarr(z_ice_max+1)
  beta_t = fltarr(z_ice_max+1)
  gamma_t =fltarr(z_ice_max+1)
  br = fltarr(z_ice_max+1)
  ar = fltarr(z_ice_max+1)
  cr = fltarr(z_ice_max+1)
  b=fltarr(z_ice_max+1)
  betadot=fltarr(z_ice_max+1)
  bdot=fltarr(z_ice_max+1)
  sol_stable=fltarr(z_ice_max+1)
  snowthickAWS    = fltarr(M)
  c_i         = fltarr(z_ice_max+1)
  dH_internal     = fltarr(elev_bins)
  dH_melt_internal  = fltarr(z_ice_max+1,M)
  mass  = fltarr(z_ice_max+1,M)
  mass_i  = fltarr(z_ice_max+1,M)
  mass_w  = fltarr(z_ice_max+1,M)
  mass_rho  = fltarr(z_ice_max+1,M)
  es_ice_surf     = fltarr(elev_bins,M)
  GF          = fltarr(z_ice_max+1)
  GFsurf        = fltarr(elev_bins,M)
  H_melt        = fltarr(elev_bins,M)
  H_rain        = fltarr(elev_bins,M)
  H_subl        = fltarr(elev_bins,M)
  H_surf        = fltarr(elev_bins,M)
  H_surfp        = fltarr(elev_bins,M)
  H_snow        = fltarr(elev_bins,M)
  k_eff       = fltarr(z_ice_max+1)
  L         = fltarr(elev_bins,M)
  zeta         = fltarr(elev_bins,M)
  LHF         = fltarr(elev_bins,M)
  meltflux      = fltarr(elev_bins,M)
  meltflux_internal = fltarr(elev_bins,M)
  precip        = fltarr(elev_bins,M)
  q_surf        = fltarr(elev_bins,M)
  rainfall      = fltarr(elev_bins,M)
  refreezing      = fltarr(z_ice_max+1)
  retention     = fltarr(z_ice_max+1)
  liq_max_V     = fltarr(z_ice_max+1)
  liq_V         = fltarr(z_ice_max+1)
  refreezing_tot    = fltarr(elev_bins,M)
  retention_tot   = fltarr(elev_bins,M)
  rho         = fltarr(z_ice_max+1,M)
  rho_d         = fltarr(z_ice_max+1)
  rho_snow      = fltarr(elev_bins,M)
  runoff        = fltarr(elev_bins,M)
  SHF         = fltarr(elev_bins,M)
  SRnet       = fltarr(z_ice_max+1)
  ;subl_ice     = fltarr(elev_bins,M)
  ;subl_sn      = fltarr(elev_bins,M)
  snowfall      = fltarr(elev_bins,M)
  snowthick     = fltarr(elev_bins,M)
  T_ice       = fltarr(elev_bins,z_ice_max+1,M)   ; sub-surface temperatures (snow and ice)
  T_iced       = fltarr(z_ice_max+1)   ; dummy sub-surface temperatures (snow and ice)
  T_ice_d       = fltarr(z_ice_max+1)   ; dummy sub-surface temperatures (snow and ice)
  T_rain        = fltarr(elev_bins,M)
  Tsurf       = fltarr(elev_bins,M)
  jpgrnd=z_ice_max+1
  tsoil      =fltarr(jpgrnd)
  grndc      =fltarr(jpgrnd)
  grndd      =fltarr(jpgrnd)
  slwc       =fltarr(jpgrnd)
  snic       =fltarr(jpgrnd)
  zrfrz      =fltarr(jpgrnd)
  snowc      =fltarr(jpgrnd)
  rhofirn    =fltarr(jpgrnd)
  zd1V       =fltarr(jpgrnd)
  cdelV      =fltarr(jpgrnd)
  zdz1       =fltarr(jpgrnd)
  zdz2       =fltarr(jpgrnd)
  zkappa     =fltarr(jpgrnd)
  zcapa      =fltarr(jpgrnd)

  tsoil[*]=T_ice_AWS+T_0
  ;pT_ice=tsoil
  rhofirn[*]=rho_ice
  grndc[*]=T_ice_AWS+T_0;0
  grndd[*]=0.00001;1;0
  grndcapc=0
  slwc[*]=0.;1
  snic[*]=dz_ice;0.1;rho_ice*dz_ice;0
  snowc[*]=0.;1
  zrogl=0.
  zrfrz[*]=0
  sn=0
  Slush=0.
  supimp=0.
  H_surf[*,*] = 0
  H_surfp[*,*] = 0
  SHF[*,*] = 0    ; no turbulent heat fluxes for wind speeds below the threshold WS_lim
  LHF[*,*] = 0
  L[*,*] = -999
  zeta[*,*] = 0

  mean_T        = fltarr(elev_bins)
  mean_RH       = fltarr(elev_bins)
  mean_WS       = fltarr(elev_bins)
  mean_SRin     = fltarr(elev_bins)
  mean_LRin     = fltarr(elev_bins)
  mean_SHF      = fltarr(elev_bins)
  mean_LHF      = fltarr(elev_bins)
  mean_GFsurf     = fltarr(elev_bins)
  mean_rainHF     = fltarr(elev_bins)
  mean_Tsurf      = fltarr(elev_bins)
  mean_Tice10m    = fltarr(elev_bins)
  mean_snowthick    = fltarr(elev_bins)
  runoff_tot      = fltarr(elev_bins)

  for j=0,elev_bins-1 do begin
    mean_T[j] = mean(T[j,*])
    rho_snow[j,*] = 350.;625. + 18.7*(mean_T[j]-T_0) + 0.293*(mean_T[j]-T_0)^2    ; surface snow density by Reeh et al. 2005
    ;  rho_snow[j,*] = 481. + 4.834*(mean_T[j]-T_0)    ; surface snow density by Kuipers-Munneke et al. 2015
  endfor
  ;T_ice_AWS =mean_T[0]
  ; SPECIFIC HUMIDITY & SATURATION -----------------------------------------------------------------------
  es_wtr = 10.^(-7.90298*(T_100/T-1.) + 5.02808 * ALOG10(T_100/T) $   ; saturation vapour pressure above 0 C (hPa)
    - 1.3816E-7 * (10.^(11.344*(1.-T/T_100))-1.) $
    + 8.1328E-3*(10.^(-3.49149*(T_100/T-1)) -1.) + ALOG10(es_100))
  es_ice = 10.^(-9.09718 * (T_0 / T - 1.) - 3.56654 * ALOG10(T_0 / T) + 0.876793 * (1. - T / T_0) + ALOG10(es_0))   ; saturation vapour pressure below 0 C (hPa)
  q_sat = eps * es_wtr/(pres-(1-eps)*es_wtr)    ; specific humidity at saturation (incorrect below melting point)
  freezing = where(T lt T_0)            ; replacing saturation specific humidity values below melting point
  if total(freezing) gt -1 then q_sat[freezing] = eps * es_ice[freezing]/(pres[freezing]-(1-eps)*es_ice[freezing])
  supersaturated = where(RH gt 100)       ; replacing values of supersaturation by saturation
  if total(supersaturated) gt -1 then RH[supersaturated] = 100.
  q = RH*q_sat/100                ; specific humidity in kg/kg

  ; PRECIPITATION -----------------------------------------------------------------------
  snowfall[*,*] = 0
  rainfall[*,*] = 0
  T_rain[*,*] = T_0
  dpd=1.0
  snowing = where(LRin_AWS ge dpd*sigma*T_AWS^4 and T-T_0 le T_solidprecip)   ; assuming equal precipitation rate over entire transect
  raining = where(LRin_AWS ge dpd*sigma*T_AWS^4 and T-T_0 gt T_solidprecip)
  ;snowing = where(LRin ge sigma*T^4 and T-T_0 le T_solidprecip)        ; assuming extrapolated values can be used for precipitation estimates
  ;raining = where(LRin ge sigma*T^4 and T-T_0 gt T_solidprecip)
  if total(snowing) gt -1 then snowfall[snowing] = prec_rate*rho_water/rho_snow[snowing]/3600.*dt_obs/dev ; in m of snow
  if total(raining) gt -1 then rainfall[raining] = prec_rate/3600.*dt_obs/dev           ; in m of water
  if total(raining) gt -1 then T_rain[raining] = T[raining]
  ; NB: The temperature of the newly accumulated snow is set by the surface energy balance.
  ; The temperature of the rain water is assumed equal to the air temperature (or 0 C when air temperature is below freezing).

  ; START OF SPATIAL LOOP -----------------------------------------------------------------------
  prevgress = 0
  for j=0,elev_bins-1 do begin

    ; INITIAL SUB-SURFACE PROFILES OF TEMPERATURE AND DENSITY -----------------------------------------------------------------------
    if elev_bins eq 1 then snowthick[j,0] = snowthick_ini else if elev[j] gt ELA then snowthick[j,0] = snowthick_ini + gradsnowthick*(elev[j]-ELA)
    if elev_bins eq 1 then T_ice[j,*,0] = T_ice_AWS + T_0 else T_ice[j,*,0] = T_ice_AWS + T_0 + gradTice*(elev[j]-elev_AWS)
    subsurfmelt = where(T_ice[j,*,0] gt T_0)
    if total(subsurfmelt gt -1) then T_ice[j,subsurfmelt,0] = T_0 ; removing non-freezing temperatures
    rho[*,0] = rho_ice
    z_icehorizon = floor(snowthick[j,0]/dz_ice)
    ;depth = findgen(z_ice_max)*dz_ice
    if snowthick[j,0] gt dz_ice*z_ice_max then z_icehorizon = z_ice_max
    if snowthick[j,0] gt 0 then begin
      for z=0,z_ice_max do begin
        rho[z,0] = rho_snow[j,0] ;+ (rho_ice-rho_snow[j,0])*z*dz_ice/snowthick[j,0]   ; including snow in the density array
        if z+1 ge z_icehorizon then break
      endfor
    endif
    ;toodense = where(rho gt rho_ice)
    ;if total(toodense) gt -1 then rho[toodense] = rho_ice
    c_i = 152.456 + 7.122*T_ice[j,*,0]                ; specific heat of ice (perhaps a slight overestimation for near-melt ice temperatures (max 48 J/kg/K))
    Tsurf[*,*] = T_0                        ; initial surface temperature in iteration for balanced energy components

    ; START OF TIME LOOP -----------------------------------------------------------------------
    z_0=z0_ice
    for k=0l,M-1 do begin
      EB_prev = 1.0                 ; reset value after previous time step
      dTsurf = dTsurf_ini             ; reset value after previous time step
      if k gt 0 then begin
        snowthick[j,k] = snowthick[j,k-1]       ; will be updated at end of time loop
        if elev_bins gt 1 then snowthick_AWS = interpol(snowthick[*,k],elev,elev_AWS) - snowthick[*,0] else snowthick_AWS = snowthick[j,k] - snowthick[j,0]
        z_WS = H_WS - snowthick_AWS
        z_T  = H_T - snowthick_AWS
        if z_T le 0.1 then z_T=0.1
        if z_WS le 0.1 then z_WS=0.1
      endif else begin
        z_WS = H_WS ; - snowthick_ini
        z_T  = H_T  ; - snowthick_ini
        if z_T le 0 then z_T=0
        if z_WS le 0 then z_WS=0
      endelse

      if k gt 0 then rho[*,k] = rho[*,k-1]
      K_eff = 0.021 + 2.5e-6*rho[*,k]^2         ; effective conductivity by Anderson 1976, is ok at limits

      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization
      ; Here comes the albedo parameterization

      ; shortwave radiation snow&ice penetration
      if elev_bins ne 1 and k gt 0 then begin         ; only used measured albedo at the AWS site
        if snowthick[j,k] gt 0 then SRout_mdl[j,k] = alb_snow*SRin[j,k] else SRout_mdl[j,k] = alb_ice*SRin[j,k]
      endif
      SRnet[0:z_icehorizon] = (SRin[j,k]-SRout_mdl[j,k])*exp(-ext_snow*findgen(z_icehorizon+1)*dz_ice)
      if z_icehorizon lt z_ice_max then begin
        SRnet[z_icehorizon+1:z_ice_max] = (SRin[j,k]-SRout_mdl[j,k])*exp(-ext_snow*snowthick[j,k])*exp(-ext_ice*((findgen(z_ice_max-z_icehorizon)+z_icehorizon+1)*dz_ice-snowthick[j,k]))
      endif ; else SRnet = (SRin[j,k]-SRout_mdl[j,k])*exp(-ext_ice*findgen(z_ice_max+1)*dz_ice)
      ; snow&ice temperature rise due to shortwave radiation absorption
      if k gt 0 then begin
        c_i = 152.456 + 7.122*T_ice[j,*,k-1]        ; Specific heat of ice (a slight overestimation for near-melt T (max 48 J kg-1 K-1))
        T_ice[j,0:z_ice_max-1,k] = T_ice[j,0:z_ice_max-1,k-1] + dt_obs/float(dev)/rho[0:z_ice_max-1,k]/c_i[0:z_ice_max-1]*(SRnet[0:z_ice_max-1]-SRnet[1:z_ice_max])/dz_ice
        T_ice[j,z_ice_max,k] = T_ice[j,z_ice_max,0]
      endif
      subsurfmelt = where(T_ice[j,*,k] gt T_0)
      nosubsurfmelt = where(T_ice[j,*,k] le T_0)
      dT_ice = T_ice[j,*,k] - T_0
      dT_ice[nosubsurfmelt] = 0
      meltflux_internal_temp = rho[*,k]*c_i*dT_ice/dt_obs*float(dev)*dz_ice
      meltflux_internal[j,k] = total(meltflux_internal_temp[1:z_ice_max])
      dH_melt_internal[*,k] = -c_i*dT_ice/L_fus*dz_ice
      dH_melt_internal[0,k] = 0
      if total(subsurfmelt gt -1) then T_ice[j,subsurfmelt,k] = T_0   ; removing non-freezing temperatures
      ; Reduce sub-surface density due to melt? Will most likely cause model instability...

      theta = T[j,k] + z_T*g/c_pd

      for findbalance=0,iter_max_EB-1 do begin              ; start iteration of surface temperature

        ; SENSIBLE AND LATENT HEAT FLUX -----------------------------------------------------------------------
; z_0 new parameterization start
;        k1=0.005;001
;        k2=0;.015*sin(!pi*day(k)/367);0.0075;0.01
; z_0 new parameterization end
        ; Roughness length scales for snow or ice - initial guess
        if snowthick[j,k] gt 0 then begin
;          if SRin[j,k] eq 0 then SRin[j,k]=0.1
;          if SRout[j,k] eq 0 then SRout[j,k]=0.1
              z_0 = z0_snow
;          if hour[k] eq 12 then z_0 = k1+k2*tanh(1-(SRout[j,k]/SRin[j,k]))*tanh(1*(1-(SRout[j,k]/SRin[j,k])));
          u_star = kappa*WS[j,k]/alog(z_WS/z_0)
          Re = u_star * z_0 / nu[j,k]
          ;    print,u_star
          if Re le 0.135 then range = 0
          if Re gt 0.135 and Re lt 2.5 then range = 1
          if Re ge 2.5 then range = 2
          z_h = z_0 * exp(ch1[range] + ch2[range]*alog(Re) + ch3[range]*(alog(Re))^2)   ; smooth surfaces: Andreas 1987
          z_q = z_0 * exp(cq1[range] + cq2[range]*alog(Re) + cq3[range]*(alog(Re))^2)
        endif else begin
;          if SRin[j,k] eq 0 then SRin[j,k]=0.1
;          if SRout[j,k] eq 0 then SRout[j,k]=0.1
              z_0 = z0_ice
          ;if hour[k] eq 12 then z_0 = k1+k2*tanh(1-(SRout[j,k]/SRin[j,k]))*tanh(1*(1-(SRout[j,k]/SRin[j,k])));
          u_star = kappa*WS[j,k]/alog(z_WS/z_0)
          Re = u_star * z_0 / nu[j,k]
          z_h = z_0 * exp(1.5 - 0.2*alog(Re) - 0.11*(alog(Re))^2)   ; rough surfaces: Smeets & Van den Broeke 2008
          z_q = z_h
        endelse
        if WS[j,k] eq 0 then begin
          z_h = 1e-10
          z_q = 1e-10
        endif

        es_ice_surf[j,k] = 10.^(-9.09718 * (T_0 / Tsurf[j,k] - 1.) - 3.56654 * ALOG10(T_0 / Tsurf[j,k]) + 0.876793 * (1. - Tsurf[j,k] / T_0) + ALOG10(es_0))
        q_surf[j,k]  = eps * es_ice_surf[j,k]/(pres[j,k]-(1-eps)*es_ice_surf[j,k])

        L[j,k] = 10e4
        if theta ge Tsurf[j,k] and WS[j,k] ge WS_lim then begin   ; stable stratification
          for i=0,iter_max_flux-1 do begin
            psi_m1 = -(aa*z_0/L[j,k]  +  bb*(z_0/L[j,k]-cc/dd)*exp(-dd*z_0/L[j,k])  + bb*cc/dd)
            psi_m2 = -(aa*z_WS/L[j,k] + bb*(z_WS/L[j,k]-cc/dd)*exp(-dd*z_WS/L[j,k]) + bb*cc/dd)
            psi_h1 = -(aa*z_h/L[j,k]  +  bb*(z_h/L[j,k]-cc/dd)*exp(-dd*z_h/L[j,k])  + bb*cc/dd)
            psi_h2 = -(aa*z_T/L[j,k]  +  bb*(z_T/L[j,k]-cc/dd)*exp(-dd*z_T/L[j,k])  + bb*cc/dd)
            psi_q1 = -(aa*z_q/L[j,k]  +  bb*(z_q/L[j,k]-cc/dd)*exp(-dd*z_q/L[j,k])  + bb*cc/dd)
            psi_q2 = -(aa*z_T/L[j,k]  +  bb*(z_T/L[j,k]-cc/dd)*exp(-dd*z_T/L[j,k])  + bb*cc/dd)
            if snowthick[j,k] gt 0 then begin
                        z_0 = z0_snow
              ;if hour[k] eq 12 then z_0 = k1+k2*tanh(1-(SRout[j,k]/SRin[j,k]))*tanh(1*(1-(SRout[j,k]/SRin[j,k])));
              u_star  = kappa * WS[j,k] / (alog(z_WS / z_0) - psi_m2 + psi_m1)
              Re = u_star * z_0 / nu[j,k]
              if Re le 0.135 then range = 0
              if Re gt 0.135 and Re lt 2.5 then range = 1
              if Re ge 2.5 then range = 2
              z_h = z_0 * exp(ch1[range] + ch2[range]*alog(Re) + ch3[range]*(alog(Re))^2)   ; smooth surfaces: Andreas 1987
              z_q = z_0 * exp(cq1[range] + cq2[range]*alog(Re) + cq3[range]*(alog(Re))^2)
              if (z_h lt 1e-10) then z_h = 1e-10
              if (z_q lt 1e-10) then z_q = 1e-10
            endif else begin
                        z_0 = z0_ice
              ;if hour[k] eq 12 then z_0 = k1+k2*tanh(1-(SRout[j,k]/SRin[j,k]))*tanh(1*(1-(SRout[j,k]/SRin[j,k])));
              u_star  = kappa * WS[j,k] / (alog(z_WS / z_0) - psi_m2 + psi_m1)
              Re = u_star * z_0 / nu[j,k]
              z_h = z_0 * exp(1.5 - 0.2*alog(Re) - 0.11*(alog(Re))^2)   ; rough surfaces: Smeets & Van den Broeke 2008
              if (z_h lt 1e-10) then z_h = 1e-10
              z_q = z_h
            endelse
            th_star = kappa * (theta - Tsurf[j,k]) / (alog(z_T / z_h) - psi_h2 + psi_h1)
            q_star  = kappa * (q[j,k] - q_surf[j,k]) / (alog(z_T / z_q) - psi_q2 + psi_q1)
            SHF[j,k]  = rho_atm[j,k] * c_pd  * u_star * th_star
            LHF[j,k]  = rho_atm[j,k] * L_sub * u_star * q_star
            L_prev  = L[j,k]
            L[j,k]    = u_star^2 * theta*(1 + ((1 - eps)/eps)*q[j,k]) / (g * kappa * th_star*(1 + ((1 - eps)/eps)*q_star))
            if abs((L_prev - L[j,k]) / L_prev) lt L_dif then break
          endfor
        endif
        if theta lt Tsurf[j,k] and WS[j,k] ge WS_lim then begin   ; unstable stratification
          for i=0,iter_max_flux-1 do begin
            x1      = (1 - gamma * z_0  / L[j,k])^0.25
            x2      = (1 - gamma * z_WS / L[j,k])^0.25
            y1      = (1 - gamma * z_h  / L[j,k])^0.5
            y2      = (1 - gamma * z_T  / L[j,k])^0.5
            yq1     = (1 - gamma * z_q  / L[j,k])^0.5
            yq2     = (1 - gamma * z_T  / L[j,k])^0.5
            psi_m1  = alog( ((1 + x1)/2)^2 * (1 + x1^2)/2) - 2*atan(x1) + !pi/2
            psi_m2  = alog( ((1 + x2)/2)^2 * (1 + x2^2)/2) - 2*atan(x2) + !pi/2
            psi_h1  = alog( ((1 + y1)/2)^2 )
            psi_h2  = alog( ((1 + y2)/2)^2 )
            psi_q1  = alog( ((1 + yq1)/2)^2 )
            psi_q2  = alog( ((1 + yq2)/2)^2 )
            if snowthick[j,k] gt 0 then begin
                        z_0 = z0_snow
              ;if hour[k] eq 12 then z_0 = k1+k2*tanh(1-(SRout[j,k]/SRin[j,k]))*tanh(1*(1-(SRout[j,k]/SRin[j,k])));
              u_star  = kappa * WS[j,k] / (alog(z_WS / z_0) - psi_m2 + psi_m1)
              Re = u_star * z_0 / nu[j,k]
              if Re le 0.135 then range = 0
              if Re gt 0.135 and Re lt 2.5 then range = 1
              if Re ge 2.5 then range = 2
              z_h = z_0 * exp(ch1[range] + ch2[range]*alog(Re) + ch3[range]*(alog(Re))^2)   ; smooth surfaces: Andreas 1987
              z_q = z_0 * exp(cq1[range] + cq2[range]*alog(Re) + cq3[range]*(alog(Re))^2)
              if (z_h lt 1e-10) then z_h = 1e-10
              if (z_q lt 1e-10) then z_q = 1e-10
            endif else begin
                        z_0 = z0_ice
              ;if hour[k] eq 12 then z_0 = k1+k2*tanh(1-(SRout[j,k]/SRin[j,k]))*tanh(1*(1-(SRout[j,k]/SRin[j,k])));
              u_star  = kappa * WS[j,k] / (alog(z_WS / z_0) - psi_m2 + psi_m1)
              Re = u_star * z_0 / nu[j,k]
              z_h = z_0 * exp(1.5 - 0.2*alog(Re) - 0.11*(alog(Re))^2)   ; rough surfaces: Smeets & Van den Broeke 2008
              if (z_h lt 1e-10) then z_h = 1e-10
              z_q = z_h
            endelse
            th_star = kappa * (theta - Tsurf[j,k]) / (alog(z_T / z_h) - psi_h2 + psi_h1)
            q_star  = kappa * (q[j,k]  - q_surf[j,k]) / (alog(z_T / z_q) - psi_q2 + psi_q1)
            SHF[j,k]  = rho_atm[j,k] * c_pd  * u_star * th_star
            LHF[j,k]  = rho_atm[j,k] * L_sub * u_star * q_star
            L_prev  = L[j,k]
            L[j,k]    = u_star^2 * theta*(1 + ((1 - eps)/eps)*q[j,k]) / (g * kappa * th_star*(1 + ((1 - eps)/eps)*q_star))
            if abs((L_prev - L[j,k]) / L_prev) lt L_dif then break
          endfor
        endif
        if WS[j,k] lt WS_lim then begin
          u_star = -999
          th_star = -999
          q_star = -999
          L[j,k] = -999
          z_h = 1e-10
          z_q = 1e-10
          psi_m1 = 0
          psi_m2 = -999
          psi_h1 = 0
          psi_h2 = -999
          psi_q1 = 0
          psi_q2 = -999
        endif

        ; SURFACE ENERGY BUDGET -----------------------------------------------------------------------

        meltflux[j,k] = SRnet[0] - SRnet[1] + LRin[j,k] - em*sigma*Tsurf[j,k]^4-(1-em)*LRin[j,k] + SHF[j,k] + LHF[j,k]  $
          -(k_eff[1])*(Tsurf[j,k]-T_ice[j,1,k])/dz_ice + rho_water*c_w*rainfall[j,k]*dev/dt_obs*(T_rain[j,k]-Tsurf[j,k])
        if meltflux[j,k] ge 0 and Tsurf[j,k] eq T_0 then break    ; stop iteration for melting surface
        if abs(meltflux[j,k]) lt EB_max then begin
          meltflux[j,k] = 0
          break                           ; stop iteration when energy components in balance
        endif
        if meltflux[j,k]/EB_prev lt 0 then dTsurf = 0.5*dTsurf    ; make surface temperature step smaller when it overshoots EB=0
        EB_prev = meltflux[j,k]
        if meltflux[j,k] lt 0 then Tsurf[j,k] = Tsurf[j,k] - dTsurf else Tsurf[j,k] = Tsurf[j,k] + dTsurf
      endfor
      if findbalance eq iter_max_EB and abs(meltflux[j,k]) ge 10*EB_max then print,'Problem closing energy budget:',j,k,meltflux[j,k],Tsurf[j,k]-T_0,T[j,k]-T_0,RH[j,k],WS[j,k],SRnet[0],LRin[j,k],SHF[j,k],LHF[j,k],-(k_eff[1])*(Tsurf[j,k]-T_ice[j,1,k])/dz_ice
      if findbalance eq iter_max_EB and abs(meltflux[j,k]) ge 10*EB_max then print,u_star,th_star,q_star,z_0,z_h,z_q,Re,L[j,k],snowthick[j,k]
      zeta[j,k]=z_WS/L[j,k]
      ; MASS BUDGET -----------------------------------------------------------------------
      ; Melt in m of snow and ice, surface and internal
      dH_melt = -meltflux[j,k]*dt_obs/dev/L_fus;/rhofirn[0]
      ;  dH_melt = -meltflux[j,k]*dt_obs/dev/L_fus/rho[0,k]
      if k gt 0 then H_melt[j,k] = H_melt[j,k-1] + dH_melt/rho_water ;+ total(dH_melt_internal[*,k])
      ;if (dH_melt lt -dz_ice) then stop
      ; Sublimation / deposition in m of snow and ice         ; no evaporation / condensation
      ;  dH_subl = LHF[j,k]*dt_obs/dev/L_sub/rho[0,k]          ; positive LHF -> deposition -> dH_subl positive
      dH_subl = LHF[j,k]*dt_obs/dev/L_sub;/rhofirn[0]          ; positive LHF -> deposition -> dH_subl positive
      if k gt 0 then H_subl[j,k] = H_subl[j,k-1] + dH_subl/rho_water

      ; Total precipitation in m of water
      if elev[j] lt prec_cutoff then begin
        snowfall[j,k] = snowfall[j,k]/2
        rainfall[j,k] = rainfall[j,k]/2
      endif

      if k gt 0 then H_snow[j,k] = H_snow[j,k-1] + snowfall[j,k]*rho_snow[j,k]/rho_water
      if k gt 0 then H_rain[j,k] = H_rain[j,k-1] + rainfall[j,k]

      ; Surface height change in m of snow and ice
      dH_surf = (dH_melt + dH_subl + snowfall[j,k]*rho_snow[j,k])/rho_water   ; no surface height change due to internal melt
      ;  dH_surf = dH_melt + dH_subl + snowfall[j,k] ;+ (total(rho[*,k-1])-total(rho[*,k]))*dz_ice       ; no surface height change due to internal melt
      if k gt 0 then H_surf[j,k] = H_surf[j,k-1] + (dH_melt + dH_subl + snowfall[j,k]*rho_snow[j,k])/rho_water ;+ total(dH_melt_internal[*,k])    ; but internal melt does contribute (eventually)
      ;  if k gt 0 then H_surf[j,k] = H_surf[j,k-1] + dH_melt + dH_subl + snowfall[j,k] ;+ total(dH_melt_internal[*,k])    ; but internal melt does contribute (eventually)
      if k gt 0 then begin
        if snowthick[j,k-1] gt 0 then snowthick[j,k] = snowthick[j,k-1] + (snowfall[j,k]*rho_snow[j,k] + dH_melt + dH_subl)/rho_water
        ;    if snowthick[j,k-1] gt 0 then snowthick[j,k] = snowthick[j,k-1] + snowfall[j,k] + dH_melt + dH_subl
        if snowthick[j,k-1] eq 0 then snowthick[j,k] = snowfall[j,k]
        if snowthick[j,k-1] lt 0 then print,'Alarm! Negative snow depth!'
      endif
      if snowthick[j,k] lt 0 then snowthick[j,k]=0
      ;  T_ice[j,*,k]=pT_ice

      ; SUB-SURFACE TEMPERATURES -----------------------------------------------------------------------
      GF[1:z_ice_max] = -k_eff[1:z_ice_max]*(T_ice[j,0:z_ice_max-1,k]-T_ice[j,1:z_ice_max,k])/dz_ice
      GFsurf[j,k] = GF[1]
      pTsurf=Tsurf[j,k]
      pT_ice=T_ice[j,*,k]
      grndhflx=GFsurf[j,k]
      zsn= (dH_subl + snowfall[j,k]*rho_snow[j,k])/rho_water;dH_subl+(snowfall[j,k])*rho_snow[j,k]/rho_water;snowfall[j,k]
      snmel=(-dH_melt)/rho_water;*rhofirn[0]/rho_water;*dz_ice
      raind=rainfall[j,k]
      Tdeep=T_ice_AWS+T_0
      ElevGrad=0.01

      pzrogl=0.
      subsurface_hirham, pTsurf, pT_ice, grndc, grndd, grndcapc,$
        grndhflx, slwc, snic, pzrogl, snowc, zsn, snmel, raind,   $
        zrfrz, sn, rhofirn, Tdeep, Slush, ElevGrad, supimp, jpgrnd, dz_ice

      if k eq 0 then H_surfp[j,k]=0.
      if k gt 0 then H_surfp[j,k]=H_surfp[j,k-1] - snmel+zsn;*dz_ice;/dz_ice
      T_ice[j,*,k]=pT_ice
      rho[*,k]=rhofirn
      refreezing_tot[j,k] = total(zrfrz[*])
      retention_tot[j,k] = total(slwc[*])
      runoff[j,k] = pzrogl ;- refreezing_tot[j,k] - retention_tot[j,k]
      z_icehorizon = floor(snowthick[j,k]/dz_ice)


      ; MODEL RUN PROGRESS -----------------------------------------------------------------------
      ;  progress = fix(10.*(float(j)+float(k+1)/M)/elev_bins)
      ;  runtimehrs = (SYSTIME(1)-starttime)*10/progress/3600
      ;  runtimemins = (runtimehrs - fix(runtimehrs))*60
      ;  runtimesecs = (runtimemins - fix(runtimemins))*60
      ;  if progress ne prevgress then print,10*progress," %. Estimated run time (h:m:s):",fix(runtimehrs),fix(runtimemins),fix(runtimesecs)
      ;  prevgress = progress
      ; if (day[k] gt 152 and year[k] eq 12) then stop
      if (k mod 24 eq 0) then print,time(k)+2000,"  day of year:", day(k) ; print daily (24) time progress for k being hourly
    endfor  ; END OF TIME LOOP -----------------------------------------------------------------------

    rainHF = rho_water*c_w*rainfall*dev/dt_obs*(T_rain-Tsurf)

    ;----------------------------------------------------------------------------------------------------
    ; INTERPOLATION BACK TO ORIGINAL TIME RESOLUTION
    M = M/dev
    if dev ne 1 then begin
      ipfac2 = M+1
      time        = interpol(time,ipfac2)
      pres[j,*]     = interpol(pres[j,*],ipfac2)
      T[j,*]      = interpol(T[j,*],ipfac2)
      RH[j,*]     = interpol(RH[j,*],ipfac2)
      WS[j,*]     = interpol(WS[j,*],ipfac2)
      SRin[j,*]     = interpol(SRin[j,*],ipfac2)
      SRout_mdl[j,*]  = interpol(SRout_mdl[j,*],ipfac2)
      LRin[j,*]     = interpol(LRin[j,*],ipfac2)
      Tsurf[j,*]    = interpol(Tsurf[j,*],ipfac2)
      L[j,*]      = interpol(L[j,*],ipfac2)
      GFsurf[j,*]   = interpol(GFsurf[j,*],ipfac2)
      SHF[j,*]      = interpol(SHF[j,*],ipfac2)
      LHF[j,*]      = interpol(LHF[j,*],ipfac2)
      rainHF[j,*]   = interpol(rainHF[j,*],ipfac2)
      meltflux[j,*]   = interpol(meltflux[j,*],ipfac2)
      meltflux_internal[j,*] =  interpol(meltflux_internal[j,*],ipfac2)
      H_aws[j,*]    = interpol(H_aws[j,*],ipfac2)
      H_surf[j,*]   = interpol(H_surf[j,*],ipfac2)
      H_melt[j,*]   = interpol(H_melt[j,*],ipfac2)
      H_snow[j,*]   = interpol(H_snow[j,*],ipfac2)
      H_subl[j,*]   = interpol(H_subl[j,*],ipfac2)
      runoff[j,*]   = interpol(runoff[j,*],ipfac2)*dev
      rainfall[j,*]   = interpol(rainfall[j,*],ipfac2)*dev
      snowfall[j,*]   = interpol(snowfall[j,*],ipfac2)*dev
      snowthick[j,*]  = interpol(snowthick[j,*],ipfac2)
    endif

    ;----------------------------------------------------------------------------------------------------
    ; WRITING TO FILES -----------------------------------------------------------------------

    openw,10,directory+'output\'+stationname+runID+'.txt'
    printf,10,' Time Year Day Hour P_hPa T_C RH_% WS_ms-1 SRin_Wm-2 SRout_Wm-2 SRout_mdl_Wm-2 LRin_Wm-2 LRout_Wm-2 SHF_Wm-2 LHF_Wm-2 GF_Wm-2 rainHF_Wm-2 MEsurf_Wm-2 MEint_Wm-2 Hsurf_m Hmelt_m Hsubl_m Runoff_m SnowThck_m SnowAcc_m Rain_m Albedo Tsurf_C L_Ob'
    for k=0l,rows-1 do begin
      printf,10, format='(100F15.5)', $
        time[k],year[k],day[k],hour[k],pres[j,k],T[j,k]-T_0,RH[j,k],WS[j,k],SRin[j,k],SRout[j,k],SRout_mdl[j,k], $
        LRin[j,k],em*sigma*Tsurf[j,k]^4+(1-em)*LRin[j,k],SHF[j,k],LHF[j,k],GFsurf[j,k],rainHF[j,k],meltflux[j,k],meltflux_internal[j,k], $
        H_surf[j,k],H_melt[j,k],H_subl[j,k],runoff[j,k],snowthick[j,k],snowfall[j,k],rainfall[j,k],SRout[j,k]/SRin[j,k],Tsurf[j,k]-T_0,L[j,k]
    endfor
    close,10

    ; PLOTTING -----------------------------------------------------------------------

    set_plot,'ps'
    !p.background = 255
    !p.color = 0
    !p.thick = 2
    ;!p.charsize = 2
    !p.charthick = 2
    ;!x.range = [160,280]
    ;!y.range = [0,200]
    ;!x.title='Time (years)'
    ;!y.title=''
    ;!x.minor = 4
    ;!x.tickinterval = 20
    ;!y.minor = 4
    ;!y.tickinterval = 10
    !y.margin=[2,1]

    if elev_bins eq 1 then elev_start = elev_AWS
    elev_str = string(j*elev_binsize+elev_start)

    !p.multi = [0,1,8,0,0]
    device,filename=directory+'output\'+stationname+runID+'_meteo.ps',xsize=18,ysize=27,yoffset=0.5
    plot,time,pres[j,*],ytitle='Air Pressure'
    plot,time,T[j,*]-T_0,ytitle='Temperature'
    plot,time,RH[j,*],ytitle='Relative Humidity'
    plot,time,WS[j,*],ytitle='Wind Speed'
    plot,time,SRin[j,*],ytitle='Shortwave Radiation'
    oplot,time,SRout[j,*]
    oplot,time,SRout_mdl[j,*]
    plot,time,LRin[j,*],ytitle='Longwave Radiation'
    plot,time,H_aws[*],ytitle='Surface height',min_value=-998
    plot,time,SRout[j,*]/SRin[j,*],ytitle='Albedo'
    device,/close

    !p.multi = [0,1,9,0,0]
    device,filename=directory+'output\'+stationname+runID+'_SEB.ps',xsize=18,ysize=27,yoffset=0.5
    plot,time,SRin[j,*]-SRout[j,*],ytitle='Shortwave Radiation'
    plot,time,LRin[j,*]-em*sigma*Tsurf[j,*]^4+(1-em)*LRin[j,*],ytitle='Longwave Radiation'
    plot,time,SHF[j,*],ytitle='Sensible Heat Flux'
    plot,time,LHF[j,*],ytitle='Latent Heat Flux'
    plot,time,GFsurf[j,*],ytitle='Sub-surface Heat Flux'
    plot,time,rainHF[j,*],ytitle='Rain Heat Flux'
    plot,time,meltflux[j,*],ytitle='Surface melt energy'
    plot,time,meltflux_internal[j,*],ytitle='Sub-surface melt energy'
    plot,time,meltflux[j,*]+meltflux_internal[j,*],ytitle='Total melt energy'
    device,/close

    time_AWS = time
    H_AWS = H_AWS[0] - H_AWS

    ; Start validation at 0 m 
    H_st = H_st[0] - H_st
    H_pt = H_pt - H_pt[0]

    !p.multi = [0,1,6,0,0]
    device,filename=directory+'output\'+stationname+runID+'_SMB.ps',xsize=18,ysize=27,yoffset=0.5
    plot,time,H_surf[j,*],ytitle='Surface height (m of snow and ice)',yrange=[-20,2]
    ;oplot,time_AWS,H_st[*],linestyle=1
    ;oplot,time_AWS,H_aws[*],linestyle=3
    oplot,time_AWS,H_pt[*],linestyle=2
    plot,time,H_melt[j,*],ytitle='Melt (m of snow and ice)'
    plot,time,H_subl[j,*],ytitle='Sublimation / deposition (m of snow and ice)'
    plot,time,snowthick[j,*],ytitle='Snow layer thickness (m of snow)'
    plot,time,runoff[j,*],ytitle='Runoff (m of water)'
    plot,time,H_snow[j,*],ytitle='Precipitation (m of water equivalent)'
    oplot,time,H_rain[j,*],linestyle=1
    device,/close

    !p.multi = [0,1,2,0,0]
    device,filename=directory+'output\'+stationname+runID+'_validation.ps',xsize=18,ysize=27,yoffset=0.5
    plot,time,H_surf[j,*]*rho_water/rho_ice,ytitle='Surface height (m of snow and ice)',yrange=[-10,2]
    ;oplot,time_AWS,H_st,psym=3,linestyle=3
    ;oplot,time_AWS,H_pt,psym=3,linestyle=4
    oplot,time_AWS,H_pt,psym=3,linestyle=4
    ;oplot,time_AWS,H_st,linestyle=3
    ;oplot,time_AWS,H_pt,linestyle=4
    plot,((LRout-(1-em)*LRin[j,*])/em/sigma)^0.25-T_0,Tsurf[j,*]-T_0,xtitle='Measured surface temperature (C)',ytitle='Modelled surface temperature (C)',psym=3 ;,xrange=[150,350],yrange=[150,350]
    oplot,[-100,100],[-100,100]
    device,/close

    M = long(rows)*dev
    !p.multi = [0,1,2,0,0]
    device,filename=directory+'output\'+stationname+runID+'_subsurf.ps',xsize=18,ysize=27,yoffset=0.5
    plot,rho[*,0],-findgen(z_ice_max+1)*dz_ice,xtitle='Snow and ice density (kg/m3)',xrange=[0,1000]
    oplot,rho[*,M/4],-findgen(z_ice_max+1)*dz_ice,linestyle=1
    oplot,rho[*,M/2],-findgen(z_ice_max+1)*dz_ice,linestyle=2
    oplot,rho[*,M*3/4],-findgen(z_ice_max+1)*dz_ice,linestyle=3
    oplot,rho[*,M-1],-findgen(z_ice_max+1)*dz_ice,linestyle=4
    if M gt 365*24*dev then oplot,rho[*,long(dev)*365*24],-findgen(z_ice_max+1)*dz_ice,linestyle=0
    plot,T_ice[j,*,0]-T_0,-findgen(z_ice_max+1)*dz_ice,xtitle='Temperature snow and ice (C)',xrange=[-20,0]
    oplot,T_ice[j,*,M/4]-T_0,-findgen(z_ice_max+1)*dz_ice,linestyle=1
    oplot,T_ice[j,*,M/2]-T_0,-findgen(z_ice_max+1)*dz_ice,linestyle=2
    oplot,T_ice[j,*,M*3/4]-T_0,-findgen(z_ice_max+1)*dz_ice,linestyle=3
    oplot,T_ice[j,*,M-1]-T_0,-findgen(z_ice_max+1)*dz_ice,linestyle=4
    if M gt 365*24*dev then oplot,T_ice[j,*,long(dev)*365*24]-T_0,-findgen(z_ice_max+1)*dz_ice,linestyle=0
    device,/close
;    !p.multi = [0,1,1,0,0]
;    device,filename=directory+'output\'+stationname+runID+'_obu_L.ps',xsize=18,ysize=27,yoffset=0.5
;    plot,(time-12)*366,zeta[j,*],ytitle='L',yrange=[0,1],xrange=[153,244]
;    device,/close

    mean_RH[j] = mean(RH[j,*])
    mean_WS[j] = mean(WS[j,*])
    mean_SRin[j] = mean(SRin[j,*])
    mean_LRin[j] = mean(LRin[j,*])
    mean_SHF[j] = mean(SHF[j,*])
    mean_LHF[j] = mean(LHF[j,*])
    mean_GFsurf[j] = mean(GFsurf[j,*])
    mean_rainHF[j] = mean(rainHF[j,*])
    mean_Tsurf[j] = mean(Tsurf[j,*])
    mean_Tice10m[j] = mean(T_ice[j,10/dz_ice,*])
    mean_snowthick[j] = mean(snowthick[j,*])
    runoff_tot[j] = total(runoff[j,*])

  endfor  ; END OF SPATIAL LOOP -----------------------------------------------------------------------

  runtimehrs = fix((SYSTIME(1)-starttime)/3600)
  runtimemins = fix((SYSTIME(1)-starttime)/60 - 60*runtimehrs)
  runtimesecs = fix((SYSTIME(1)-starttime) - 3600*runtimehrs - 60*runtimemins)
  print,'Run time (h:m:s):', runtimehrs,runtimemins,runtimesecs
  print,"Done..."
  beep
END
pro subsurface_hirham, pts, ptsoil    $
  ,pgrndc,     pgrndd,    pgrndcapc   $
  ,pgrndhflx,  pslwc ,    psnic       $
  ,zrogl,      psnowc,    zsn         $
  ,zsnmel,     zraind,    zrfrz       $
  ,psn,        prhofirn,  pTdeep      $
  ,pSlush,     pElevGrad, zsupimp     $
  ,jpgrnd,     pcdel

  ;kproma    =horizontal loop/index parameter
  ;kbdim     =horizontal dimension
  ;pts       =surface temperature [K]
  ;ptsoil    =temperature in ice and snow [K]
  ;pgrndc    =Coefficient used in the ptsoil calculation
  ;pgrndd    =Coefficient used in the ptsoil calculation
  ;pgrndcapc =Heat capacity of the uppermost layer [j/m^2/K]
  ;pgrndhflx =Heat flux at the surface [W/m^2]
  ;ldglac    =ice mask(1=ice point)
  ;pslwc     =Snow liquid water content (kg/m^2)
  ;psnic     =Snow refrozen ice content (kg/m^2)
  ;zrogl     =Runoff at glacier points (rain and melt, but no calving)
  ;psnowc    =Snow content (kg/m^2)
  ;zsn       =Snow budget (snowfall-subl-melt)
  ;zsnmel    =Snow/ice melt ()
  ;zraind    =Total rain
  ;zrfrz     =Refrozen meltwater (kg/m^2)
  ;psn       =Snow depth [m water equivalent]
  ;prhofirn  =density of snow
  ;pTdeep    =fixed temperature below model layers
  ;pSlush    =Liquid water available for runoff
  ;pElevGrad =Elevation gradient at point (m/m) used in runoff parameterization
  ;zsupimp   =Superimposed ice content (kg/m^2)
  ;pjrow     =Grid point index rows
  ;  use declare, only : jpgrnd
  ;   jpgrnd=25

  ;   ! Physical parameters
  stbo  = 5.67e-8   ;! Stephan-Boltzmann constant in W/m2/K4
  alf = 333700.  ;! latent heat for fusion in J/kg
  tmelt = 273.15  ;! melting temperature of ice/snow
  rhoice = 900.  ;! standard density of glacier ice (fixed at 917 kg m-3)
  rhosnow = 330. ;! density of freshly fallen snow (currently fixed at 330 kg m-3) ! PLA densification (Feb 2015)
  rhoh2o = 1000. ;! density of water (1000 kg m-3)
  rh2oice  = rhoh2o / rhoice
  ;! PLA densification (Feb 2015)
  rho_pco = 830.     ;! Pore close-off density (kg/m3) (determines when layer becomes impermeable)
  ;!pemi = 0.98    ! standard value of surface emissivity over glaciers (Oerlemans reference)
  pemi = 0.996    ;! HIRHAM value of surface emissivity
  zdifiz = 12.e-7        ;! temperature diffusivity of ice  [m**2/s]
  g_grav = 9.82          ;! acceleration of gravity (m/s2)
  ;! PLA densification (Feb 2015)
  liqmax = 0.02       ;! Maximum liquid water fraction in the snow (pr pore space vol, SOMARS/RACMO uses 2%)
  ;! PLA densification (Feb 2015)
  icemax = 0.85       ;! Ice fraction that determines when to stop counting snow depth
  ;! ## Albedo parameterizations ##
  ;! albedo_hhcr:
  calbmns_cr = 0.65      ;! minimum (glacier, snow on ice) albedo Van de Wal and Oerlemans 1994
  ;!(they get there through ageing, we do by warming as in Roeckner 2003)
  calbmxs_cr = 0.85      ;! maximum (glacier, snow on ice) albedo Van de Wal and Oelremans 1994
  calbbare_cr = 0.40     ;! bare ice albedo (Cuffey and Paterson say 0.35 for clean ice,
  ;! van de Wal and Oerlemans say 0.55.
  ;! Oerlemans and Knap:
  calb_ok_ds =  0.032*rhosnow/rhoh2o   ;! Oerlemans and Knap suggest 3.2 cm snow depth
  ;! PLA superimposed IF (Feb 2015):
  cro_1 = 0.33 *24.*3600. ;! (units seconds) Parameter c1 from Zuo & oerlemans and Reijmer et al
  ;! controlling runoff timescale. Z&O say 1.5 but Reijmer (MAR) say 0.33
  cro_2 = 25.  *24.*3600. ;! (units seconds) Parameter c2 from Zuo & oerlemans
  cro_3 = 140.                ;! (unitless)      Parameter c3 from Zuo & oerlemans
  ;! PLA RFO fresh snow density (May 2015) rhosnow = a+b*elev+c*lat+d*lon
  a_rho = 328.35327      ;! a parameter in RFO's surface snow density parameterization
  b_rho =  -0.049376258  ;! b parameter in RFO's surface snow density parameterization
  c_rho =   1.0426993    ;! c parameter in RFO's surface snow density parameterization
  d_rho =  -0.11186342   ;! d parameter in RFO's surface snow density parameterization

  delta_time = 3600. ;!1200.   ;! Model time step (s)
  delta_time_flux = 21600.   ;! Flux average time step (s)
  delta_time_total= 86400.    ;! Integration time
  nstep    = delta_time_total / delta_time
  nstepflx = delta_time_flux  / delta_time
  ;! Tuning parameters
  eps = 0.1
  ;! Missing value
  missing_value = -9999.
  ;! Small number for layer calculations
  smallno = 1.e-12
  cdel=fltarr(jpgrnd) ;! Soil layer thicknesses
  cmid=fltarr(jpgrnd)  ;! Mid point of soil layers
  rcdel=fltarr(jpgrnd) ;! Reciprocal cdel
  RHOS=prhofirn
  ;;  use tools,   only : ice_heats, calc_snowdepth1D
  ;cdel (0) = 0.0650d0
  ;cdel (1) = 0.2540d0
  ;cdel (2) = 0.2984d0
  ;cdel (3) = 0.3507d0
  ;cdel (4) = 0.4120d0
  ;cdel (5) = 0.4842d0
  ;cdel (6) = 0.5689d0
  ;cdel (7) = 0.6684d0
  ;cdel (8) = 0.7854d0
  ;cdel (9) = 0.9229d0
  ;cdel (10) = 1.0844d0
  ;cdel (11) = 1.2741d0
  ;cdel (12) = 1.4971d0
  ;cdel (13) = 1.7591d0
  ;cdel (14) = 2.0669d0
  ;cdel (15) = 2.4286d0
  ;cdel (16) = 2.8537d0
  ;cdel (17) = 3.3530d0
  ;cdel (18) = 3.9398d0
  ;cdel (19) = 4.6293d0
  ;cdel (20) = 5.4394d0
  ;cdel (21) = 6.3913d0
  ;cdel (22) = 7.5098d0
  ;cdel (23) = 8.8240d0
  ;cdel (24) = 10.3682d0
  cdel(*)=pcdel
  cdelsum = 0.
  for jk = 0,jpgrnd-1 do begin
    cdelsum = cdelsum + cdel(jk)
    cmid(jk) = cdelsum - ( cdel(jk) / 2. )
    rcdel(jk) = 1./cdel(jk)
  endfor

  ; Arguments

  ;INTEGER ::
  ;  kproma, kbdim, pjrow

  ;;  REAL(8) ::                                                             &
  ;;  pts=fltarr(kbdim)
  ;  ptsoil=fltarr(jpgrnd)
  ;  pgrndc=fltarr(jpgrnd)
  ;  pgrndd=fltarr(jpgrnd)
  ;;  pgrndcapc=fltarr(kbdim)
  ;;  pgrndhflx=fltarr(kbdim)
  ;  pslwc=fltarr(jpgrnd)
  ;  psnic=fltarr(jpgrnd)
  ;;  zrogl=fltarr(kbdim)
  ;  zrfrz=fltarr(jpgrnd)
  ;;  zraind=fltarr(kbdim)
  ;  psnowc=fltarr(jpgrnd)
  ;;  zsn=fltarr(kbdim)
  ;;  zsnmel=fltarr(kbdim)
  ;;  psn=fltarr(kbdim)
  ;  prhofirn=fltarr(jpgrnd)
  ;;  pTdeep=fltarr(kbdim)
  ;;  pSlush=fltarr(kbdim)
  ;;  pElevGrad=fltarr(kbdim)
  ;;  zsupimp=fltarr(kbdim)
  ;
  ;  ;LOGICAL ::  ldglac(kbdim)
  ;;  ldglac=intarr(kbdim)
  ;
  ;  ; Local variables
  ;;  zso_cond=fltarr(kbdim)
  ;;  zso_capa=fltarr(kbdim)
  ;;  Nslush=intarr(kbdim)     ; Slush layer index
  ;;  z1=fltarr(kbdim)
  zd1V=fltarr(jpgrnd)
  cdelV=fltarr(jpgrnd)
  zdz1=fltarr(jpgrnd)
  zdz2=fltarr(jpgrnd)
  zkappa=fltarr(jpgrnd)
  zcapa=fltarr(jpgrnd)
  zcapaF=fltarr(jpgrnd)
  zsn_capaF=fltarr(jpgrnd)
  zsn_condF=fltarr(jpgrnd)
  ptsoild=fltarr(jpgrnd)
  ;  ptsoild(*)=0;fltarr(jpgrnd)
  ;  liqexcess=0.;fltarr(jpgrnd)
  ;  DpslwcB = 0. ;
  ;  DpsnowB = 0. ;
  ;  DpsnicB = 0. ;
  ;  liqout = 0.
  ;  liqro  = 0.
  ;  zdel  = 0.
  ;rho_frozen =0.
  ;massadj=0.
  ;cfrozenp1=0.
  ;potret=0.
  ;DPSLWCMAB=0.
  ;DPSNICMAB=0.
  ;DPSNOWMAB=0.
  ;LIQIN=0.
  ;LOCRHOSNOW=0.
  ;MATOP=0.
  ;PSNOWC_OLD=0.
  ;TSNOWIN=0.
  ;ZDZ2SUB=0.
  ;DPSLWCT=0.
  ;DPSNICT=0.
  ;DPSNOWT=0.
  ;MABOT=0.
  ;TRAININ=0.
  ;Z1=0.
  ;NSLUSH=0.
  ;POTSIFORM=0.
  ;SIFORM=0.
  ;SNOWV=0.
  ;SNOWV1=0.
  ;SNOWV2=0.
  ;T_RUNOFF=0.
  ;ZSNIN=0.
  ;ZX1=0.
  ;ZX11=0.
  ;ZX12=0.
  ;ZX2=0.
  ;TOTALV=0.
  ;SLUSH_RUNOFF=0.
  ;ICEV=0.
  ;ICELAYER=0.
  ;DTDZ=0.
  ;ki=0.
  ;zsnout=0.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  call ice_heats, kproma, kbdim, zso_cond, zso_capa
  ;subroutine ice_heats (kproma, kbdim, zso_cond, zso_capa)
  ;use declare, only : rhoice, zdifiz
  ;
  ;IMPLICIT NONE

  ;  ARGUMENTS
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: zso_cond(kbdim), zso_capa(kbdim)

  ;  local Variables
  ;INTEGER :: jl
  ;REAL(8) :: zrici
  cpice = 2.108e3       ; Specific heat capacity of ice J/kg/K

  zrici = cpice(0.)*rhoice

  ;*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
  ;*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
  ;for jl = 0,kproma-1 do begin
  zso_capa = zrici
  zso_cond = zso_capa*zdifiz
  ;endfor
  ;end subroutine ice_heats

  ; ###########################################
  ; ###########################################
  ;function cpice(TEMP)
  ;  IMPLICIT NONE
  ;  real(8), intent(in)  :: TEMP
  ;  real(8) :: cpice

  ;  cpice = 2.108d+03       ; Specific heat capacity of ice J/kg/K

  ;end function cpice

  ; ###########################################
  ; PLA densification (Feb 2015)
  ;function zsn_capaF(RHOS)
  ;  IMPLICIT NONE
  ;  real(8), intent(in)  :: RHOS
  ;  real(8) :: zsn_capaF
  ; The F in the name zsn_capaF stands for "function"
  ;  zsn_capaF = cpice(0.0d0)*RHOS  ; snow vol. heat capacity   [J/m**3/K]
  zsn_capaF = cpice(0.)*prhofirn  ; snow vol. heat capacity   [J/m**3/K]

  ;end function zsn_capaF

  ; ###########################################
  ; PLA densification (Feb 2015
  ;function zsn_condF(RHOS)
  ;  use declare, only : rhoice, zdifiz, rhoh2o
  ;  IMPLICIT NONE
  ;  real(8), intent(in)  :: RHOS
  ;  real(8)  :: zsn_condF
  ; The F in the name zsn_condF stands for "function"
  ;  zsn_condF = cpice(0.0d0)*rhoice*zdifiz*((RHOS/rhoh2o)^1.88)
  zsn_condF = cpice(0.)*rhoice*zdifiz*((prhofirn/rhoh2o)^1.88)
  ; snow thermal conductivity [J/s/m/K]
  ; Yen, Y. (1981), Review of thermal properties of snow, ice and sea ice,
  ; Rep. 81-10, U. S. Army Cold Reg. Res. and Eng. Lab. (CRREL), Hanover, N. H

  ;end function zsn_condF

  ; Update temperatures based on previous time steps coefficients for heat conduction
  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  ; Upper layer
  ptsoil(0)= pts
  ; Lower layers
  ;for jk = 0,jpgrnd-2 do begin
  ;ptsoil(jk+1) = pgrndc(jk)+pgrndd(jk)*ptsoil(jk)
  ;endfor
  for jk = 0,jpgrnd-2 do begin
    ptsoil(jk+1) = pgrndc(jk)+pgrndd(jk)*ptsoil(jk)
  endfor
  ;stop;endif
  ;endfor

  ;END SUBROUTINE tsoil_diffusion
  ;densification, kproma, kbdim           $
  ;  , ldglac, pslwc, psnowc , prhofirn, ptsoil

  zdtime = delta_time

  ; PLA densification (Feb 2015)
  ; Densification of layers (RSF+PLA Jan-Feb 2015):
  ; Vionnet et al. (GMD, 2012) densification eqs 5-9
  ;for jl = 0,kproma-1 do begin
  ;  if (ldglac(jl) eq 1) then begin
  for jk = 0,jpgrnd-1 do begin
    ; Pressure in Pa (evaluated at layer mid-point).
    ; We have neglected surface slope in eqn 6.
    sigma_pres = cmid(jk)*rhoh2o*g_grav
    f1f2 = 4./( 1.+60.*(pslwc(jk)/(psnowc(jk)+smallno)) $
      * (prhofirn(jk) / rhoh2o) ) ; f2 = 4
    eta_firn   = f1f2 * 7.62237e6 * ( prhofirn(jk) / 250. ) $
      * exp( 0.1*(tmelt - ptsoil(jk)) + 0.023*prhofirn(jk) )
    prhofirn(jk) = prhofirn(jk) + zdtime*prhofirn(jk) * sigma_pres /eta_firn
    if (prhofirn(jk) gt rhoice) then prhofirn(jk) = rhoice
    if (psnowc(jk) lt 1.e-3) then prhofirn(jk) = rhoice
  endfor
  ;  endif
  ;endfor

  ;  call snowfall_shiftmass, kproma , kbdim     $
  ;  , ldglac, zsn, psnowc, psnic, pslwc    $
  ;  , prhofirn, ptsoil, pts, zrogl, pjrow
  ;SUBROUTINE snowfall_shiftmass (kproma , kbdim &
  ;, ldglac, zsn, psnowc, psnic, pslwc &
  ;, prhofirn, ptsoil, pts, zrogl, pjrow )
  ;
  ;use declare, only : jpgrnd, smallno, cdel, rcdel, tmelt, a_rho, b_rho, c_rho, d_rho &
  ;, lat, lon, elev
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim, pjrow
  ;REAL(8) :: zsn(kbdim) &
  ;, psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd), pslwc(kbdim,jpgrnd) &
  ;, prhofirn(kbdim,jpgrnd), ptsoil(kbdim,jpgrnd), pts(kbdim), zrogl(kbdim)
  ;
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: zsnin, zdel
  ;REAL(8) :: DpsnowB, DpsnowT, DpsnicB, DpsnicT, DpslwcB, DpslwcT
  ;REAL(8) :: Tsnowin
  ;REAL(8) :: locrhosnow

  ;!********************************************************************************
  ;! Move the layer interfaces upwards due to snowfall
  ;!-----------------------------------------------------------------------------
  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  if (zsn gt 0.) then begin ;! zsn means snowfall in timestep
    zsnin = zsn
    ;! Special Treatment for the first layer
    while (zsnin gt smallno ) do begin
      ptsoild=ptsoil
      ;! 1st layer
      if (zsnin lt cdel(0)) then zdel=zsnin else zdel=cdel(0);
      ;zdel = min(zsnin, minval(cdel)) ;! Change in surface height (layer mass) in loop

      jk = 0
      DpsnowB = -(zdel*rcdel(jk) * psnowc(jk)) ;! Snow:  Change through bottom:+gain, -lost
      DpsnicB = -(zdel*rcdel(jk) * psnic(jk))  ;! Ice content
      DpslwcB = -(zdel*rcdel(jk) * pslwc(jk))  ;! Liquid water content
      ;! Mass update
      psnowc(jk) = psnowc(jk) + zdel + DpsnowB
      psnic(jk)  = psnic(jk)         + DpsnicB
      pslwc(jk)  = pslwc(jk)         + DpslwcB
      ;! Temperature update
      if (pts lt tmelt) then Tsnowin=pts else Tsnowin=tmelt
      ;Tsnowin = min(pts(jl),tmelt) ; PLA: May be updated later to use tas
      ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(Tsnowin-ptsoil(jk))
      ;! PLA densification (Feb 2015)
      ;! Density update. New surface snow has rhosnow dependent on place and altitude
      if (zdel gt smallno and psnowc(jk) gt smallno) then begin
        locrhosnow = 350.;a_rho + b_rho*elev(jl,pjrow) + c_rho*lat(jl,pjrow) + d_rho*lon(jl,pjrow)
        prhofirn(jk) = psnowc(jk) / ( (psnowc(jk)-zdel )/prhofirn(jk) + zdel/locrhosnow)
      endif

      ;! 2nd -  layer
      for jk = 1,jpgrnd-1 do begin
        ;! Deeper layers need to have snic and slwc shuffled too
        DpsnowT = - DpsnowB ;! Coming from the layer above
        DpsnicT = - DpsnicB ;! Coming from the layer above
        DpslwcT = - DpslwcB ;! Coming from the layer above
        DpsnowB = -(zdel*rcdel(jk) * psnowc(jk)) ;! Change:-gain, +lost
        DpsnicB = -(zdel*rcdel(jk) * psnic(jk))
        DpslwcB = -(zdel*rcdel(jk) * pslwc(jk))
        ;! Mass update
        psnowc(jk) = psnowc(jk) + DpsnowT + DpsnowB
        psnic(jk)  = psnic(jk)  + DpsnicT + DpsnicB
        pslwc(jk)  = pslwc(jk)  + DpslwcT + DpslwcB
        ;! Temperature update (zdel in at top, zdel out at bottom)
        ;ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(ptsoil(jk-1)-ptsoil(jk))
        ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(ptsoild(jk-1)-ptsoil(jk))
        ;! PLA densification (Feb 2015)
        ;! Density update
        if (abs(DpsnowT) gt smallno and psnowc(jk) gt smallno) then begin
          prhofirn(jk) = psnowc(jk) / ( (psnowc(jk)-DpsnowT )/prhofirn(jk) + DpsnowT/prhofirn(jk-1) )
        endif
      endfor
      zsnin = zsnin - zdel
      ;! The liquid out of the bottom layer is added to runoff, DpslwcB is negative number:
      ;! PLA superimposed IF (Feb 2015): This mass is shifted out of bottom of model
      ;! due to addition of mass at top. Does not contribute to slush bucket.
      ;! PLA (Feb 27 2015) mass-shifted water out of bottom no longer counts as runoff
      ;! zrogl(jl) = zrogl(jl) - DpslwcB
      ;stop
    endwhile ;END DO ; END OF WHILE LOOP
  endif
  ;endif ; IF LGLAC
  ;endfor ; HORIZONTAL JL (KPROMA) LOOP

  ;  call rainfall_shiftmass, kproma , kbdim     $
  ;  , ldglac, zraind, psnowc, psnic, pslwc $
  ;  , prhofirn, ptsoil, pts, zrogl
  ;SUBROUTINE rainfall_shiftmass (kproma , kbdim &
  ;, ldglac, zraind, psnowc, psnic, pslwc &
  ;, prhofirn, ptsoil, pts, zrogl )
  ;
  ;use declare, only : jpgrnd, smallno, cdel, rcdel, tmelt
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: zraind(kbdim) &
  ;, psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd), pslwc(kbdim,jpgrnd) &
  ;, prhofirn(kbdim,jpgrnd), ptsoil(kbdim,jpgrnd), pts(kbdim), zrogl(kbdim)
  ;
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: zsnin, zdel
  ;REAL(8) :: DpsnowB, DpsnowT, DpsnicB, DpsnicT, DpslwcB, DpslwcT
  ;REAL(8) :: Trainin
  ;
  ;!********************************************************************************
  ;! Move the layer interfaces upwards due to rain
  ;! Supersaturation is allowed. Handled in melt_perc
  ;!-----------------------------------------------------------------------------
  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  if (zraind gt 0.) then begin ;! zsn means rainfall in timestep
    zsnin = zraind
    ;! Special Treatment for the first layer
    while (zsnin gt smallno ) do begin
      ptsoild=ptsoil
      ;! 1st layer
      if (zsnin lt cdel(0)) then zdel=zsnin else zdel=cdel(0);  zdel = min(klud, dz_ice)
      ;if (zsnin lt minval(cdel)) then zdel=zsnin else zdel=minval(cdel)
      ;zdel = min(zsnin, minval(cdel))

      jk = 0
      DpsnowB = -(zdel*rcdel(jk) * psnowc(jk)) ;! Change through bottom:+gain, -lost
      DpsnicB = -(zdel*rcdel(jk) * psnic(jk))
      DpslwcB = -(zdel*rcdel(jk) * pslwc(jk))
      ;! Mass update
      psnowc(jk) = psnowc(jk)       + DpsnowB
      psnic(jk)  = psnic(jk)        + DpsnicB
      pslwc(jk)  = pslwc(jk) + zdel + DpslwcB
      ;! Temperature update
      if (pts gt tmelt) then Trainin=pts else Trainin=tmelt
      ;Trainin = max(pts(jl),tmelt)  ;! PLA: May be updated later to use tas
      ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(Trainin-ptsoil(jk))
      ;! PLA densification (Feb 2015)
      ;! Density update: No change in top layer due to rain.

      ;! 2nd -  layer
      for jk = 1,jpgrnd-1 do begin
        ;! Deeper layers need to have snic and slwc shuffled too
        DpsnowT = - DpsnowB ;! Coming from the layer above
        DpsnicT = - DpsnicB ;! Coming from the layer above
        DpslwcT = - DpslwcB ;! Coming from the layer above
        DpsnowB = -(zdel*rcdel(jk) * psnowc(jk)) ;! Change:-gain, +lost
        DpsnicB = -(zdel*rcdel(jk) * psnic(jk))
        DpslwcB = -(zdel*rcdel(jk) * pslwc(jk))
        ;! Mass update
        psnowc(jk) = psnowc(jk) + DpsnowT + DpsnowB
        psnic(jk)  = psnic(jk)  + DpsnicT + DpsnicB
        pslwc(jk)  = pslwc(jk)  + DpslwcT + DpslwcB
        ;! Temperature update (zdel in at top, zdel out at bottom)
        ;ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(ptsoil(jk-1)-ptsoil(jk))
        ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(ptsoild(jk-1)-ptsoil(jk))
        ;! PLA densification (Feb 2015)
        ;! Density update
        if (abs(DpsnowT) gt smallno and psnowc(jk) gt smallno) then begin
          prhofirn(jk) = psnowc(jk) / ( (psnowc(jk)-DpsnowT )/prhofirn(jk) + DpsnowT/prhofirn(jk-1) )
        endif
      endfor;ENDDO
      zsnin = zsnin - zdel
      ;! The liquid out of the bottom layer is added to runoff, DpslwcB is negative number:
      ;! PLA superimposed IF (Feb 2015): This mass is shifted out of bottom of model
      ;! due to addition of mass at top. Does not contribute to slush bucket.
      ;! PLA (Feb 27 2015) mass-shifted water out of bottom no longer counts as runoff
      ;!zrogl(jl) = zrogl(jl) - DpslwcB
    endwhile;END DO ;! END OF WHILE
  endif
  ;endif ;! IF LGLAC
  ;endfor ;! HORIZONTAL JL (KPROMA) LOOP

  ;  call melt_perc, kproma , kbdim      $
  ;  , ldglac, prhofirn             $
  ;  , psnowc, psnic, pslwc         $
  ;  , ptsoil, zsnmel               $
  ;  , pSlush, pTdeep
  ;SUBROUTINE melt_perc (kproma , kbdim &
  ;, ldglac, prhofirn              &
  ;, psnowc, psnic, pslwc          &
  ;, ptsoil, zsnmel                &
  ;, pSlush, pTdeep )
  ;
  ;use declare, only : jpgrnd, smallno, liqmax, rhoh2o, rhoice, rcdel, tmelt, rho_pco
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: prhofirn(kbdim,jpgrnd)
  ;REAL(8) :: psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd), pslwc(kbdim,jpgrnd)
  ;REAL(8) :: ptsoil(kbdim,jpgrnd), zsnmel(kbdim)
  ;REAL(8) :: pSlush(kbdim),  pTdeep(kbdim)
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: zdel
  ;REAL(8) :: liqexcess, zsnout, cfrozen, cfrozenp1
  ;REAL(8) :: DpsnowB, DpsnowT, DpsnicB, DpsnicT, DpslwcB, DpslwcT
  ;REAL(8) :: potret, massadj
  ;REAL(8) :: DpsnowMAB, DpsnicMAB, DpslwcMAB
  ;REAL(8) :: liqro, liqin, liqout, MAbot, MAtop
  ;REAL(8) :: psnowc_old, rho_frozen, liqmaxM
  ;
  ;!======================================================
  ;!Here we do surface melt and liquid water percolation
  ;!======================================================
  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  jk = 0
  ;! PLA densification (Feb 2015): Both of following lines
  liqmaxM = liqmax*rhoh2o/rhoice*(rhoice/prhofirn(jk) - 1.)
  liqexcess = pslwc(jk) - liqmaxM*psnowc(jk)

  if ( (zsnmel gt 0.) or (liqexcess gt 0.0) ) then begin
    ;! Here is surface melt OR supersaturation due to rain (see above)
    ;! HIRHAM: In HIRHAM remember to write out the different phases to
    ;! restart and output files
    zsnout = zsnmel
    ;if zsnmel gt 0.001 then stop
    ;! Special Treatment for the first layer
    ;while (.true.) do begin;! Checking below to allow for entry for supersaturation
    while (zsnout ge smallno) do begin;! Checking below to allow for entry for supersaturation
      ptsoild=ptsoil
      ;!under none melting conditions
      ;! 1st layer
      jk=0
      ;! How much frozen mass do we have in top layer available for melting?
      cfrozen = psnowc(jk) + psnic(jk)
      if (zsnout lt cfrozen) then zdel=zsnout else zdel=cfrozen;
      ;if zsnmel gt 0.001 then stop
      ;zdel = min(zsnout, cfrozen)
      ;! We are moving up the layer by subtracting ice and snow
      ;! The melting is subtracted by the existing fraction from the ice and snow
      ;! DpsnowT indicates mass moves through top of the layer
      ;! DpsnowB indicates mass through the bottom of the layer
      ;! Top layer: here is where the melting takes place. Take mass from frozen phases and
      ;! add to liquid phase. Too much liquid after this? Let it percolate down.

      ;! If less melt than there is snow, take it all from the snow
      if ( zdel le psnowc(jk)) then begin
        DpsnowT = -zdel
        DpsnicT = 0.
        psnowc(jk) = psnowc(jk) + DpsnowT
      endif else begin;! If more melt than snow, take all snow and some ice
        DpsnowT = -psnowc(jk)
        DpsnicT = -1.*(zdel - psnowc(jk) )
        psnowc(jk) = 0.  ;! Include this line here to make it identically 0.
      endelse

      DpslwcT = (DpsnowT + DpsnicT)*(-1.0)
      ;! Snow update was handled inside the if loop
      psnic(jk) = psnic(jk) + DpsnicT
      pslwc(jk) = pslwc(jk) + DpslwcT


      ;! PLA densification (Feb 2015): Both of following lines
      liqmaxM = liqmax*rhoh2o/rhoice*(rhoice/prhofirn(jk) - 1.)
      potret    = liqmaxM* psnowc(jk)
      ;! liqexcess is the liquid water in a layer above the retention threshold.
      ;! Moves to layer below
      liqexcess = pslwc(jk) - potret

      DpslwcB = 0. ;
      DpsnowB = 0. ;
      DpsnicB = 0. ;! Preset the bottom fluxes
      liqout = 0.
      if ( liqexcess gt 0. ) then begin
        ;! Yes: too much liquid
        DpslwcB = -liqexcess
        cfrozenp1 =  psnowc(jk+1) + psnic(jk+1) ;!cfrozen for the next layer
        ;! here the frozen parts of the layers are shifted upwards
        ;! for downward liquid water percolation
        DpsnowB = liqexcess/cfrozenp1 * psnowc(jk+1)
        DpsnicB = liqexcess/cfrozenp1 * psnic(jk+1)
        ;! liquid down to next layer (used in temperature calc)
        liqout = liqexcess
      endif

      ;! Mass update
      psnowc(jk) = psnowc(jk) + DpsnowB
      psnic(jk)  = psnic(jk)  + DpsnicB
      pslwc(jk)  = pslwc(jk)  + DpslwcB

      ;! Temperature update, assume liquid run-down is at tmelt
      ptsoil(jk) = ptsoil(jk) - liqout*rcdel(jk)*( tmelt - ptsoil(jk+1) )

      ;! PLA densification (Feb 2015)
      ;! Density update:
      if (abs(DpsnowB) gt smallno and psnowc(jk) gt smallno) then begin
        prhofirn(jk) = psnowc(jk) / $
          ( (psnowc(jk) - DpsnowB ) /prhofirn(jk) + DpsnowB /prhofirn(jk+1))
      endif

      ;! Intervening layers: Receive liquid from above, shift ice and snow up.
      ;! Check for saturation etc.
      ;! If a perched ice layer is encountered in the column, liqexcess goes to runoff.
      ;! In that case, mass is removed from the column and subsequent interfaces must be shifted upward.
      massadj = 0.
      ;! Preset the mass adjustment bottom fluxes
      DpslwcMAB = 0. ;
      DpsnowMAB = 0. ;
      DpsnicMAB = 0.
      for jk = 1,jpgrnd-2 do begin
        liqin = liqout ;! Liquid out of bottom of previous layer
        MAtop = massadj ;! Mass adjustment upwards out of top (for temperature calc)
        ;! PLA densification (Feb 2015):
        psnowc_old = psnowc(jk) ;! saved for density update at bottom of loop

        DpsnowT = -DpsnowB - DpsnowMAB ;! Change through top: +gain, -lost
        DpsnicT = -DpsnicB - DpsnicMAB
        DpslwcT = -DpslwcB - DpslwcMAB

        psnowc(jk) = psnowc(jk) + DpsnowT
        psnic(jk)  = psnic(jk)  + DpsnicT
        pslwc(jk)  = pslwc(jk)  + DpslwcT
        ;! Vertical coordinates shifted if runoff has occurred due to perched ice layers
        ;! Mass is taken from underlying layer to balance the loss to runoff.
        DpslwcMAB = massadj * pslwc(jk+1)  * rcdel(jk+1)
        DpsnowMAB = massadj * psnowc(jk+1) * rcdel(jk+1)
        DpsnicMAB = massadj * psnic(jk+1)  * rcdel(jk+1)
        psnowc(jk) = psnowc(jk) + DpsnowMAB
        psnic(jk)  = psnic(jk)  + DpsnicMAB
        pslwc(jk)  = pslwc(jk)  + DpslwcMAB

        ;! PLA densification (Feb 2015): Both of following lines
        liqmaxM = liqmax*rhoh2o/rhoice*(rhoice/prhofirn(jk) - 1.)
        potret    = liqmaxM* psnowc(jk)
        ;! liqexcess is the liquid water in a layer above the retention threshold.
        ;! Moves to layer below
        liqexcess = pslwc(jk) - potret
        DpslwcB = 0. ;
        DpsnowB = 0. ;
        DpsnicB = 0. ;! Preset the bottom fluxes
        liqout = 0.
        liqro = 0.
        if ( liqexcess gt 0. ) then begin
          ;! PLA densification (Feb 2015)
          ;! New clause for impermeability: compare bulk rho to rho_pco
          ;!Special case: If more ice in layer than some threshold in next layer,
          rho_frozen = ( psnic(jk+1)*rhoice + psnowc(jk+1)*prhofirn(jk+1) ) / $
            (psnic(jk+1)+psnowc(jk+1))
          if (rho_frozen ge rho_pco and massadj le 0.) then begin
            ;! The next layer has too much ice to receive the liquid water. Give to runoff
            ;! Update layer interface due to removal of mass from column
            ;! PLA superimposed IF (Feb 2015):
            ;!zrogl(jl)    = zrogl(jl)    + liqexcess
            ;! GIve to slush bucket instead:
            pSlush = pSlush + liqexcess
            pslwc(jk) = pslwc(jk) - liqexcess
            DpslwcB = liqexcess * pslwc(jk+1)  * rcdel(jk+1)
            DpsnowB = liqexcess * psnowc(jk+1) * rcdel(jk+1)
            DpsnicB = liqexcess * psnic(jk+1)  * rcdel(jk+1)
            massadj = liqexcess
            liqout = 0.
            liqro = liqexcess
          endif else begin
            ;! Standard case
            ;! Yes: too much liquid, which may go to next layer
            DpslwcB = -liqexcess
            cfrozenp1 =  psnowc(jk+1) + psnic(jk+1) ;!cfrozen for the next layer
            ;! Here the frozen parts of the layers are shifted upwards for
            ;!  downward liquid water percolation
            DpsnowB = liqexcess/cfrozenp1 * psnowc(jk+1)
            DpsnicB = liqexcess/cfrozenp1 * psnic(jk+1)
            liqout = liqexcess
            liqro = 0.
          endelse
        endif
        psnowc(jk) = psnowc(jk) + DpsnowB
        psnic(jk)  = psnic(jk)  + DpsnicB
        pslwc(jk)  = pslwc(jk)  + DpslwcB

        ;! Temperature update
        MAbot = massadj ;! Mass adjustment bottom, in from below
        ptsoil(jk) = ptsoil(jk) + (   $
          liqin    * ( tmelt - ptsoil(jk)   ) $
          - liqout * ( tmelt - ptsoil(jk+1) ) $
          - liqro  *   tmelt $
          - MAtop  *  ptsoil(jk) $
          + MAbot  *  ptsoil(jk+1) $
          ) * rcdel(jk)

        ;! PLA densification (Feb 2015)
        ;! Density update:
        if (( abs(DpsnowB) gt smallno or abs(DpsnowT) gt smallno or abs(DpsnowMAB) gt smallno ) $
          and psnowc(jk) gt smallno) then begin
          prhofirn(jk) = psnowc(jk) / ( $
            psnowc_old / prhofirn(jk) + DpsnowT / prhofirn(jk) $
            + (DpsnowMAB+DpsnowB)/prhofirn(jk+1) $
            )
        endif
      endfor;END DO ;!jk = 2,jpgrnd-1

      ;! Bottom layer: Receive from above. Give excess liquid to runoff and replace by pure ice.
      jk = jpgrnd-1
      liqin = liqout ;! Liquid out of bottom of previous layer
      MAtop = massadj ;! Mass adjustment upwards out of top
      DpsnowT = -DpsnowB - DpsnowMAB ;! Change through top: +gain, -lost
      DpsnicT = -DpsnicB - DpsnicMAB
      DpslwcT = -DpslwcB - DpslwcMAB
      psnowc(jk) = psnowc(jk) + DpsnowT
      psnic(jk)  = psnic(jk)  + DpsnicT
      pslwc(jk)  = pslwc(jk)  + DpslwcT

      ;! If mass adjustment has occurred further up we need to include here:
      DpsnicMAB = massadj
      psnic(jk) = psnic(jk) + DpsnicMAB

      ;! PLA densification (Feb 2015): Both of following lines
      liqmaxM = liqmax*rhoh2o/rhoice*(rhoice/prhofirn(jk) - 1.)
      potret    = liqmaxM* psnowc(jk)
      ;! liqexcess is the liquid water in a layer above the retention threshold.
      ;! Moves to layer below
      liqexcess = pslwc(jk) - potret

      DpslwcB = 0. ;
      DpsnowB = 0. ;
      DpsnicB = 0. ;! Preset the bottom fluxes
      liqout = 0.
      liqro  = 0.
      if ( liqexcess gt 0. ) then begin
        ;! Yes: too much liquid. Replace runoff in layer with ice from infinite reservoir below.
        DpslwcB = -liqexcess
        DpsnicB =  liqexcess
        ;! PLA superimposed IF (Feb 2015):
        ;!zrogl(jl) = zrogl(jl) + liqexcess
        ;! Give to slush bucket instead:
        pSlush = pSlush + liqexcess
        liqout = liqexcess
      endif
      psnic(jk) = psnic(jk) + DpsnicB
      pslwc(jk) = pslwc(jk) + DpslwcB

      ;! Temperature update
      MAbot = massadj ;! Mass adjustment bottom, in from below
      ptsoil(jk) = ptsoil(jk) + (   $
        liqin    * ( tmelt - ptsoil(jk)   ) $
        - liqout * ( tmelt - pTdeep ) $
        - MAtop  *  ptsoil(jk) $
        + MAbot  *  pTdeep $
        ) * rcdel(jk)
      ;! PLA Tdeep (feb 2015) pTdeep(jl) changed two places in the above
      ;! PLA densification (Feb 2015)
      ;! Density update: No change in snow density, since snow can only leave the layer.
      ;if zsnmel gt 0.001 then stop

      zsnout = zsnout - zdel
      ;if (zsnout le smallno) then EXIT  ;! No more melting required

    endwhile;END DO ;! WHILE
  endif

  ;endif ;! IF LGLAC
  ;endfor ;! HORIZONTAL JL (KPROMA) LOOP

  ;END SUBROUTINE melt_perc

  ;  call calc_snowdepth1D, kproma , kbdim $
  ;  , ldglac, psn, psnowc, psnic
  ;! ###########################################
  ;subroutine calc_snowdepth1D (kproma , kbdim &
  ;, ldglac, psn, psnowc, psnic )
  ;
  ;use declare, only : jpgrnd, smallno, icemax, cdel
  ;
  ;implicit none
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd),  psn(kbdim)
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: notice

  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin

  ;! Diagnose snow depth (for use in albedo parameterizations and other places?)
  ;! Include only the snow above first perched ice layer (and include the snow of this layer)

  notice = 1. ;! Should be read as "this layer is not ice"
  ;! First the snow in the top layer
  if ( psnowc(0) lt smallno ) then begin
    ;! If "no snow" in top layer, set psn = 0 and
    ;! include no snow from subsequent layers
    psn = 0.
    notice = 0.
  endif else begin
    psn = psnowc(0)
  endelse

  ;! Then the next layers
  for jk=1,jpgrnd-1 do begin
    psn = psn + notice*psnowc(jk)
    if (psnic(jk) gt icemax*cdel(jk) ) then begin
      ;! Once this has been set to 0, it will not change again
      ;! and no more snow is counted
      notice = 0.
    endif
  endfor
  ;endif ;! IF LGLAC
  ;endfor ;! HORIZONTAL JL (KPROMA) LOOP
  ;end subroutine calc_snowdepth1D

  ;  call refreeze, kproma , kbdim       $
  ;  , ldglac                       $
  ;  , psnowc, psnic, pslwc         $
  ;  , ptsoil, zrfrz, pts
  ;SUBROUTINE refreeze (kproma , kbdim &
  ;, ldglac                       &
  ;, psnowc, psnic, pslwc         &
  ;, ptsoil, zrfrz, pts  )
  ;
  ;use declare, only : jpgrnd, tmelt, alf
  ;use tools,   only : cpice
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd), pslwc(kbdim,jpgrnd)
  ;REAL(8) :: ptsoil(kbdim,jpgrnd), zrfrz(kbdim,jpgrnd), pts(kbdim)
  ;
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: cfrozen, zpotref, coldcontent
  ;
  ;!======================================================
  ;!  Here we do the refreezing based on the cold content
  ;! of each layer converting mass from liquid to ice.
  ;!======================================================

  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  for jk=0,jpgrnd-1 do begin
    ;! Determine layer cold content and convert mass from liq to ice
    cfrozen = psnowc(jk) + psnic(jk)
    if ((tmelt-ptsoil(jk)) le 0.) then coldcontent=0. else coldcontent=(tmelt-ptsoil(jk));
    ;coldcontent = max(0., (tmelt-ptsoil(jl,jk)))
    zpotref = coldcontent * cpice(0.) * cfrozen / alf
    if (zpotref lt pslwc(jk)) then zrfrz(jk)=zpotref else zrfrz(jk)=pslwc(jk)
    ;zrfrz(jl,jk)= min(zpotref , pslwc(jl,jk))
    if ((pslwc(jk)-zrfrz(jk)) le 0.) then pslwc(jk)=0. else pslwc(jk)=(pslwc(jk)-zrfrz(jk))
    ;pslwc(jl,jk) = max(0., pslwc(jl,jk) - zrfrz(jl,jk)) ;! To avoid tiny negative numbers
    psnic(jk) = psnic(jk) + zrfrz(jk)
    ;! Update temperature with latent heat of freezing
    ptsoil(jk) = ptsoil(jk) + $
      zrfrz(jk)* alf/( cpice(0.)*(psnic(jk)-zrfrz(jk) + psnowc(jk)+smallno) )
  endfor

  ;! Update surface temperature
  pts = ptsoil(0)

  ;endif
  ;endfor
  ;END SUBROUTINE refreeze

  ;  call find_slushlevel, kproma , kbdim $
  ;  , ldglac,  psnic, Nslush
  ;SUBROUTINE find_slushlevel (kproma , kbdim &
  ;, ldglac,  psnic, Nslush)
  ;
  ;use declare, only : jpgrnd, icemax, cdel
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: psnic(kbdim,jpgrnd)
  ;LOGICAL :: ldglac(kbdim)
  ;INTEGER :: Nslush(kbdim)     ! Slush layer index. That is, layer-number of
  ;! the slush layer depth, just above first icemax-layer
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;LOGICAL :: icelayer

  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  ;! Diagnose layer on top of which the  slush bucket sits
  icelayer = 0;.false.
  for jk=0,jpgrnd-1 do begin
    if (psnic(jk) gt icemax*cdel(jk) and icelayer eq 0 ) then begin
      ;! First time icelayer is encountered, the level is recorded and
      ;! we do not enter into this if loop again:
      icelayer = 1;.true.
      Nslush = jk
    endif
  endfor
  ;! If we haven't seen ice yet, icelayer is still .false. and Nslush
  ;! is set to 0, which actually means that
  ;! it is sitting on top of the infinite subsurface layer
  if (icelayer eq 0) then Nslush = 0
  ;endif
  ;endfor
  ;END SUBROUTINE find_slushlevel

  ;  call runoff, kproma , kbdim $
  ;  , Nslush, pSLush, pElevGrad, zrogl
  ;SUBROUTINE runoff ( kproma , kbdim &
  ;, Nslush, pSLush, pElevGrad, zrogl )
  ;
  ;use declare, only : jpgrnd, smallno, cro_1, cro_2, cro_3, delta_time
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;INTEGER :: Nslush(kbdim)     ! Slush layer index. That is, layer-number of
  ;! the slush layer depth, just above first icemax-layer
  ;REAL(8) :: pSlush(kbdim), pElevGrad(kbdim), zrogl(kbdim)
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl
  ;REAL(8) :: slush_runoff, t_runoff, zdtime

  zdtime = delta_time
  ;!=======================================================================
  ;! First subtract water due to time-scale runoff (Zuo and Oerlemans 1996)
  ;! Parameters  are set as in Lefebre et al (JGR, 2003) = MAR value (Fettweis pers comm)
  ;! If top-layer is ice, we allow no slush layer and all runs off
  ;!=======================================================================

  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  if (Nslush eq 1) then begin
    slush_runoff = pSLush
  endif else begin
    t_runoff = cro_1 + cro_2 * exp(- cro_3 * pElevGrad)
    slush_runoff = pSLush/t_runoff * zdtime
  endelse
  pSlush = pSlush - slush_runoff
  if (pSlush le 0.) then pSlush=0. else pSlush=pSlush
  ;pSlush(jl) = max(0.d0, pSlush(jl)) ;! To avoid tiny negative numbers
  zrogl  = zrogl  + slush_runoff
  ;endif
  ;endfor
  ;END SUBROUTINE runoff

  ;  call superimposedice, kproma , kbdim       $
  ;  , ldglac, prhofirn, ptsoil            $
  ;  , psnowc, psnic, pslwc                $
  ;  , pSlush, Nslush, zso_cond            $
  ;  , zsupimp, zrogl
  ;SUBROUTINE superimposedice (kproma , kbdim &
  ;, ldglac, prhofirn, ptsoil            &
  ;, psnowc, psnic, pslwc                &
  ;, pSlush, Nslush, zso_cond            &
  ;, zsupimp, zrogl )
  ;
  ;use declare, only : jpgrnd, smallno, rhoh2o, rh2oice, tmelt, alf, delta_time, cdel, rcdel
  ;use tools,   only : cpice, zsn_condF
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;REAL(8) :: prhofirn(kbdim,jpgrnd)
  ;REAL(8) :: psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd), pslwc(kbdim,jpgrnd)
  ;REAL(8) :: ptsoil(kbdim,jpgrnd)
  ;REAL(8) :: pSlush(kbdim),  zso_cond(kbdim)
  ;REAL(8) :: zsupimp(kbdim), zrogl(kbdim)
  ;
  ;INTEGER :: Nslush(kbdim)     ! Slush layer index. That is, layer-number of the
  ;! slush layer depth, just above first icemax-layer
  ;LOGICAL :: ldglac(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: snowV, iceV, totalV, zx1, zx2
  ;REAL(8) :: ki, dTdz, zdtime
  ;REAL(8) :: SIform, potSIform
  ;REAL(8) :: zsnin, zdel
  ;REAL(8) :: DpsnowB, DpsnowT, DpsnicB, DpsnicT, DpslwcB, DpslwcT

  zdtime = delta_time
  ;! --------------------------------------
  ;! ----- SUPERIMPOSED ICE FORMATION -----
  ;! --------------------------------------
  ;! PLA superimposed IF (Feb 2015):

  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  ;! No SIF allowed onto inifinte sublayer, so only enter
  ;! loop if we have encountered an ice layer in the column:
  if ( (Nslush ge 2) ) then begin
    ;! Calculate SI-formation
    ;! #####################
    jk = Nslush
    snowV = 0.5*psnowc(jk)  *rhoh2o/prhofirn(jk)
    iceV  = 0.5*psnic(jk)   *rh2oice
    totalV = snowV + iceV
    zx1 = snowV / totalV
    zx2  = iceV   / totalV
    ;ki= 1.d0/( zx1/zsn_condF(prhofirn(jk)) + zx2/zso_cond )
    ki= 1./( zx1/zsn_condF(jk) + zx2/zso_cond )
    dTdz =  (tmelt-ptsoil(jk))  / totalV

    ;! The potential SIformation
    if ((ki*dTdz / (rhoh2o*alf) * zdtime) le 0.) then potSIform=0. else potSIform=(ki*dTdz / (rhoh2o*alf) * zdtime)
    ;potSIform = max(0.d0, ki*dTdz / (rhoh2o*alf) * zdtime) ;! Only interested in positive SIF

    ;! Make as much SI as there is slush available
    if (pSlush lt potSIform) then SIform=pSlush else SIform=potSIform
    ;SIform = min(pSlush(jl),potSIform)
    zsupimp = SIform

    ;! Remove from slush bucket
    ;! ########################
    pSlush = pSlush - SIform
    if (pSlush le 0.) then pSlush=0. else pSlush=pSlush
    ;pSlush(jl) = max(0.d0, pSlush(jl)) ;! To avoid tiny negative numbers
    ;! If it becomes very small, set pSlush to zero:
    IF (pSlush lt smallno ) then pSlush = 0.

    ;! Update temperature of ice layer
    ;! ###############################
    ptsoil(jk) = ptsoil(jk) + SIform * alf / (cdel(jk)*cpice(0.))

    ;! Handle addition of ice to layer and downward mass, temperature and density advection
    ;! #####################################################################################
    zsnin = SIform
    ;! Special Treatment for the first layer
    while (zsnin gt smallno ) do begin
      ptsoild=ptsoil
      ;! 1st layer
      if (zsnin lt cdel(0)) then zdel=zsnin else zdel=cdel(0)
      ;zdel = min(zsnin, minval(cdel)) ;! Mass addition split up in bits no larger than thinnest layer
      jk = Nslush
      DpsnowB = -(zdel*rcdel(jk) * psnowc(jk)) ;! Snow:  Change through bottom:+gain, -lost
      DpsnicB = -(zdel*rcdel(jk) * psnic(jk))  ;! Ice content
      DpslwcB = -(zdel*rcdel(jk) * pslwc(jk))  ;! Liquid water content
      ;! Mass update
      psnowc(jk) = psnowc(jk)        + DpsnowB
      psnic(jk)  = psnic(jk)  + zdel + DpsnicB
      pslwc(jk)  = pslwc(jk)         + DpslwcB
      ;! Temperature update
      ;! New ice added has temperature tmelt:
      ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(tmelt-ptsoil(jk))
      ;! PLA densification (Feb 2015)
      ;! Density update. No effect on snow density from added ice (snow only leaves)

      ;! Subsequent layers
      if (Nslush lt jpgrnd-1) then begin
        for jk = Nslush,jpgrnd-1 do begin ;for jk = Nslush(jl)+1,jpgrnd do begin
          ;! Deeper layers need to have snic and slwc shuffled too
          DpsnowT = - DpsnowB ;! Coming from the layer above
          DpsnicT = - DpsnicB ;! Coming from the layer above
          DpslwcT = - DpslwcB ;! Coming from the layer above
          DpsnowB = -(zdel*rcdel(jk) * psnowc(jk)) ;! Change:-gain, +lost
          DpsnicB = -(zdel*rcdel(jk) * psnic(jk))
          DpslwcB = -(zdel*rcdel(jk) * pslwc(jk))
          ;! Mass update
          psnowc(jk) = psnowc(jk) + DpsnowT + DpsnowB
          psnic(jk)  = psnic(jk)  + DpsnicT + DpsnicB
          pslwc(jk)  = pslwc(jk)  + DpslwcT + DpslwcB
          ;! Temperature update (zdel in at top, zdel out at bottom)
          ;ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(ptsoil(jk-1)-ptsoil(jk))
          ptsoil(jk) = ptsoil(jk) + (zdel*rcdel(jk))*(ptsoild(jk-1)-ptsoil(jk))
          ;! PLA densification (Feb 2015)
          ;! Density update
          if (abs(DpsnowT) gt smallno and psnowc(jk) gt smallno) then begin
            prhofirn(jk) = psnowc(jk) / $
              ( (psnowc(jk)-DpsnowT )/prhofirn(jk) + DpsnowT/prhofirn(jk-1) )
          endif
        endfor
      endif
      zsnin = zsnin - zdel
      ;! The liquid out of the bottom layer is added to runoff, DpslwcB is negative number:
      ;! PLA superimposed IF (Feb 2015): This mass is shifted out of bottom of model
      ;! due to addition of mass at top. Does not contribute to slush bucket.
      ;! PLA (Feb 27 2015) mass-shifted water out of bottom no longer counts as runoff
      ;! zrogl(jl) = zrogl(jl) - DpslwcB
    endwhile ;! END OF WHILE LOOP
  endif

  ;endif
  ;endfor
  ;END SUBROUTINE superimposedice

  ;  call update_tempdiff_params,   kproma , kbdim &
  ;  , ldglac, prhofirn, pTdeep                    &
  ;  , psnowc, psnic                               &
  ;  , ptsoil, zso_cond, zso_capa                  &
  ;  , pgrndc, pgrndd,   pgrndcapc, pgrndhflx
  ;SUBROUTINE update_tempdiff_params ( kproma , kbdim &
  ;, ldglac, prhofirn, pTdeep                    &
  ;, psnowc, psnic                               &
  ;, ptsoil, zso_cond, zso_capa                  &
  ;, pgrndc, pgrndd,   pgrndcapc, pgrndhflx )
  ;
  ;use declare, only : jpgrnd, smallno, delta_time, rhoh2o, rh2oice, cdel
  ;use tools,   only : zsn_capaF, zsn_condF
  ;
  ;IMPLICIT NONE
  ;
  ;! Arguments
  ;INTEGER :: kproma, kbdim
  ;LOGICAL :: ldglac(kbdim)
  ;REAL(8) :: prhofirn(kbdim,jpgrnd),pTdeep(kbdim)
  ;REAL(8) :: psnowc(kbdim,jpgrnd), psnic(kbdim,jpgrnd)
  ;REAL(8) :: ptsoil(kbdim,jpgrnd), zso_cond(kbdim), zso_capa(kbdim)
  ;REAL(8) :: pgrndc(kbdim,jpgrnd), pgrndd(kbdim,jpgrnd)
  ;REAL(8) :: pgrndcapc(kbdim)    ,pgrndhflx(kbdim)
  ;
  ;! Local variables
  ;INTEGER :: jl, jk
  ;REAL(8) :: zdtime
  ;REAL(8) :: z1(kbdim)
  ;REAL(8) :: zd1V(kbdim,jpgrnd), cdelV(kbdim,jpgrnd)
  ;REAL(8) :: zdz1(kbdim,jpgrnd),   zdz2(kbdim,jpgrnd)
  ;REAL(8) :: zkappa(kbdim,jpgrnd), zcapa(kbdim,jpgrnd)
  ;REAL(8) :: snowV, iceV, totalV, zx1, zx2
  ;REAL(8) :: zdz2sub
  ;REAL(8) :: snowV1, snowV2, zx11, zx12

  zdtime = delta_time

  ;!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
  for jk = 0,jpgrnd-1 do begin
    ;for jl = 0,kproma-1 do begin
    zkappa(jk) = zso_cond
    zcapa(jk)  = zso_capa
    ;endfor
  endfor

  ;!   ---------------------------------------------------------------
  ;!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
  ;!   ---------------------------------------------------------------
  ;!
  ;! PETER AND RUTH's VERSION (Note: we ignore the liquid content in all of the following)
  ;! We apparently need physical layer thicknesses (i.e., in snow and ice meters, not liquid meters):
  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  for jk = 0, jpgrnd-1 do begin
    ;! PLA densification (Feb 2015)
    cdelV(jk) = psnowc(jk)*rhoh2o/prhofirn(jk) + psnic(jk)*rh2oice + smallno
  endfor
  ;! Layer distances
  for jk = 0,jpgrnd-2 do begin
    zd1V(jk) = 1./( 0.5*(cdelV(jk+1)+cdelV(jk)) )
  endfor
  ;! Special new code for the sublayer (which is all-ice):
  zd1V(jpgrnd-1) = 1./( 0.5*(cdel(jk)*rh2oice+cdelV(jk)) )
  ;endif
  ;endfor
  ;! Calculate layerwise volume-weighted versions of capa and kappa
  ;for jl = 0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  ;! Special treatment of first layer
  ;! *Here*, the first layer is made up of layer 1 and half of layer 2:
  ;! PLA densification (Feb 2015)
  snowV1 =       psnowc(0)*rhoh2o/prhofirn(0)
  snowV2 = 0.5*psnowc(1)*rhoh2o/prhofirn(1)
  iceV  = (psnic(0) + 0.5*psnic(1)) *rh2oice
  totalV = snowV1 + snowV2 + iceV
  zx11 = snowV1 / totalV
  zx12 = snowV2 / totalV
  zx2 = iceV  / totalV
  ;zcapa(0) = zx11 * zsn_capaF(prhofirn(0)) + zx12 * zsn_capaF(prhofirn(1)) $
  ;+ zx2 * zso_capa
  zcapa(0) = zx11 * zsn_capaF(0) + zx12 * zsn_capaF(1) $
    + zx2 * zso_capa
  ;zkappa(0) = 1.d0/ $
  ;( zx11/zsn_condF(prhofirn(0)) + zx12/zsn_condF(prhofirn(1)) + zx2/zso_cond )
  zkappa(0) = 1./ $
    ( zx11/zsn_condF(0) + zx12/zsn_condF(1) + zx2/zso_cond )

  ;! Next layers:
  ;! Here, it is done midpoint-to-midpoint
  for jk = 1, jpgrnd-2 do begin;! original code only went to jpgrnd - 3
    ;! PLA densification (Feb 2015)
    snowV1 = 0.5*psnowc(jk)  *rhoh2o/prhofirn(jk)
    snowV2 = 0.5*psnowc(jk+1)*rhoh2o/prhofirn(jk+1)
    iceV  = 0.5*(psnic(jk) + psnic(jk+1)) *rh2oice
    totalV = snowV1 + snowV2 + iceV
    zx11 = snowV1 / totalV
    zx12 = snowV2 / totalV
    zx2 = iceV  / totalV
    ;zcapa(jk) = zx11 * zsn_capaF(prhofirn(jk)) + zx12 * zsn_capaF(prhofirn(jk+1)) $
    ;+ zx2 * zso_capa
    ;zkappa(jk) = 1.d0/ ( zx11/zsn_condF(prhofirn(jk)) $
    ;+ zx12/zsn_condF(prhofirn(jk+1)) + zx2/zso_cond )
    zcapa(jk) = zx11 * zsn_capaF(jk) + zx12 * zsn_capaF(jk+1) $
      + zx2 * zso_capa
    zkappa(jk) = 1./ ( zx11/zsn_condF(jk) $
      + zx12/zsn_condF(jk+1) + zx2/zso_cond )
  endfor
  ;! New code for bottom layer zcapa assuming below is ice to same thickness (at least)
  jk = jpgrnd-1
  ;! PLA densification (Feb 2015)
  snowV = 0.5 * psnowc(jk) *rhoh2o/prhofirn(jk)
  iceV = 0.5 * (psnic(jk) + cdel(jk)) *rh2oice
  totalV = snowV + iceV
  zx1 = snowV / totalV
  zx2 = iceV  / totalV
  ;zcapa(jk) = zx1*zsn_capaF(prhofirn(jk)) + zx2*zso_capa
  ;zkappa(jk) = 1.d0/( zx1/zsn_condF(prhofirn(jk)) + zx2/zso_cond )
  zcapa(jk) = zx1*zsn_capaF(jk) + zx2*zso_capa
  zkappa(jk) = 1./( zx1/zsn_condF(jk) + zx2/zso_cond )
  ;endif
  ;endfor

  ;for jl=0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  for jk=0,jpgrnd-1 do begin
    ;! Use volume version of cdel:
    zdz2(jk)=zcapa(jk)*cdelV(jk)/zdtime
  endfor
  ;endif
  ;endfor

  ;for jl=0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  ;! The following used to go jk=1,jpgrnd-1
  ;! Now we include jpgrnd, because it is needed below:
  for jk=0,jpgrnd-1 do begin
    ;! Use volume version of xd1
    zdz1(jk)=zd1V(jk)*zkappa(jk)
    ;if  ptsoil(jk) gt 273.15 then ptsoil(jk)=273.15
  endfor
  ;endif
  ;endfor

  ;! We have now calculated our new versions of zdz1, zdz2, zcapa and zkappa.
  ;! Before introducing diffusion with sub-model layer, the old code was used as it were.
  ;! Now, we do some more:
  ;
  ;! In the original version, this loop calculated c and d of jpgrnd-1 (which are in turn used
  ;! above (in next time step) to prognose tsoil(jpgrnd). Now we calculate c and d of jpgrnd.
  ;! These are not used to prognose t of jpgrnd+1 (because this isn't relevant) but
  ;! to initialize the upwards calculation of c and d:
  ;for jl=0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  ;! Sublayer is all ice (zcapa = zso_capa) and has same mass as jpgrnd-layer (cdel(jpgrnd)),
  ;! giving physical thickness cdel(jpgrnd)*rh2oice:
  zdz2sub = zso_capa*cdel(jpgrnd-1)*rh2oice/zdtime ;! corresponds to zdz2(jl,jpgrnd+1)
  z1=zdz2sub+zdz1(jpgrnd-1)
  pgrndc(jpgrnd-1)=zdz2sub*pTdeep/z1
  pgrndd(jpgrnd-1)=zdz1(jpgrnd-1)/z1

  ;! PLA Tdeep (feb 2015) pTdeep(jl) changed in the above
  ;endif
  ;endfor

  ;! This loop went jk=jpgrnd-1,2,-1 (ie, calculating c and d for 3,2,1)
  ;! Now, it goes   jk=jpgrnd,2,-1 (ie, calculating c and d for 4,3,2,1)
  ;! It thus needs zdz1(jl,jpgrnd) which is now also calculated above
  ;for jl=0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  for jk=jpgrnd-1,1,-1 do begin
    z1=1./(zdz2(jk)+zdz1(jk-1) +                    $
      zdz1(jk)*(1.-pgrndd(jk)))
    pgrndc(jk-1)=(ptsoil(jk)*zdz2(jk) +               $
      zdz1(jk)*pgrndc(jk))*z1
    pgrndd(jk-1)=zdz1(jk-1)*z1
  endfor
  ;endif
  ;endfor
  ;stop
  ;!   ---------------------------------------------------------
  ;!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
  ;!   CALORIFIC CAPACITY OF THE GROUND:
  ;!   ---------------------------------------------------------
  ;for jl=0,kproma-1 do begin
  ;if (ldglac(jl) eq 1) then begin
  pgrndhflx=zdz1(0)*(pgrndc(0)          $
    +(pgrndd(0)-1.)*ptsoil(0))
  pgrndcapc=(zdz2(0)*zdtime+                $
    zdtime * (1.-pgrndd(0)) * zdz1(0))
  ;endif
  ;endfor
  ;END SUBROUTINE update_tempdiff_params

  ;END SUBROUTINE subsurface
  ;stop
end