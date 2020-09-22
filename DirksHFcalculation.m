%----------------------------------------------------------------------------------------------------
% Calculating turbulent heat fluxes -----------------------------------------------------------------
% NB: Requires (hourly) averages. Only variables feeding into && out of this section have subscript "h", but all are hourly values.
addpath(genpath('C:\Users\bava\ownCloud\Phd_owncloud\Data\AWS\Input'))
[data] = ImportPROMICEData('KAN_U');

%%
emissivity = 0.97;

% Tsurf_obs = min(c.T_0, ((data.LongwaveRadiationUpWm2- (1 - c.em) * data.LongwaveRadiationDownWm2)/c.em/c.sigma).^(1/4));
% Tsurf_obs = (data.LongwaveRadiationUpWm2/c.em/c.sigma).^(1/4);

% % Hi B,
% %  
% % Excellent. Next step then is to compare how we calculate Tsurf. This is my calculation:
% %  
% % emissivity = 0.97
% % Tsurf = ((LRout-(1-emissivity)*LRin)/emissivity/5.67e-8)^0.25 - T_0
% % too_warm = where(Tsurf gt 0)
% % if total(too_warm) ne -1 then Tsurf[too_warm] = 0
% %  
% % Cheers, 
% % Dirk
param.station = 'KAN-U';
[c] = ImportConst(param);

Tsurf_obs = ((data.LongwaveRadiationUpWm2-(1-emissivity)*data.LongwaveRadiationDownWm2)/emissivity/5.67e-8).^0.25;
figure
hold on
plot(data.SurfaceTemperatureC(1:100)+273.15)
plot(Tsurf_obs(1:100))
% data.SurfaceTemperatureC = Tsurf_obs-273.15;


data.time =data.time+1-0*0.0417;
T_h = data.AirTemperature1C+273.15;
RH_h = data.RelativeHumidity1Perc;
WS_h = data.WindSpeed1ms;
H_instr_wind_h = data.WindSensorHeight1m;
H_instr_temp_h = data.TemperatureSensorHeight1m;
pres_h = data.AirPressurehPa;


% Constant declariation
z_0    =    0.001;  % aerodynamic surface roughness len>h for momention (assumed constant for all ice./snow surfaces)
eps    =    0.622;
es_0   =    6.1071; % saturation vapour pressure at the me<ing point (hPa)
es_100 = 1013.246  ;% saturation vapour pressure at steam point temperature (hPa)
g      =    9.82   ;% gravitational acceleration (m./s2)
gamma  =   16;     % flux profile correction (Paulson & Dyer)
kappa  =    0.4;    % Von Karman constant (0.35-0.42)
L_sub  =    2.83e6; % latent heat of sublimation (J./kg)
R_d    =  287.05   ;% gas constant of dry air
aa     =    0.7    ;% flux profile correction constants (Holtslag & De Bruin '88)
bb     =    0.75;
cc     =    5.;
dd     =    0.35;
c_pd   = 1005.   ;  % specific heat of dry air (J./kg./K)
WS_lim =    1.;
L_dif_max = 0.01;
T_0 = 273.15;
T_100 = 373.15;
% Array declaration && initial guesses
% z_WS = HaWS + 0.4;
% z_T = HaWS - 0.1
z_WS = H_instr_wind_h;
z_T = H_instr_temp_h;
Tsurf_h = data.SurfaceTemperatureC(1:length(T_h))+273.15;

rho_atm_h = 100.*pres_h./R_d./T_h ;% atmospheric density
mu_h = 18.27e-6*(291.15+120)./(T_h+120).*(T_h./291.15).^1.5; % dynamic viscosity of air (Pa s) (Sutherl&&s' equation using C = 120 K)
nu_h = mu_h./rho_atm_h ;% kinematic viscosity of air (m^2/s)
u_star = kappa.*WS_h./log(z_WS./z_0);
Re = u_star.*z_0./nu_h;
z_0h = z_0.*exp(1.5-0.2.*log(Re)-0.11.*(log(Re)).^2); % rough surfaces: Smeets & Van den Broeke 2008

z_0h( (WS_h <= 0)) = 1e-10;
es_ice_surf = 10.^(-9.09718.*(T_0./(Tsurf_h)-1.) ...
    - 3.56654.*log10(T_0./(Tsurf_h))...
    +0.876793.*(1. - Tsurf_h/T_0) + log10(es_0));
q_surf = eps.*es_ice_surf./(pres_h-(1-eps).*es_ice_surf);
es_wtr = 10.^(-7.90298.*(T_100./T_h-1.) ...
    + 5.02808 .* log10(T_100./T_h) ... % saturation vapour pressure above 0 C (hPa)
          - 1.3816E-7 .* (10.^(11.344.*(1.-T_h./T_100))-1.) ...
        + 8.1328E-3.*(10.^(-3.49149.*(T_100./T_h-1)) -1.) + log10(es_100));
es_ice = 10.^(-9.09718 .* (T_0 ./ T_h - 1.) ...
    - 3.56654 .* log10(T_0 ./ T_h) ...
    + 0.876793 .* (1. - T_h ./ T_0) + log10(es_0));   % saturation vapour pressure below 0 C (hPa)
q_sat = eps .* es_wtr./(pres_h-(1-eps).*es_wtr); % specific humidity at saturation (incorrect below me<ing point)
freezing =  (T_h < T_0); % replacing saturation specific humidity values below me<ing point
q_sat(freezing) = eps .* es_ice(freezing)./(pres_h(freezing)-(1-eps).*es_ice(freezing));
q_h = RH_h.*q_sat./100; % specific humidity in kg./kg
theta = T_h + z_T.*g./c_pd;
SHF_h = T_h ;
SHF_h(:) = 0; 
LHF_h = SHF_h;
L = SHF_h+1e5;
[~, q3] = SpecHumSat(RH_h,T_h,pres_h,c);

stable   =  and(theta >= Tsurf_h , WS_h > WS_lim); % && T ne -999 && Tsurf_h ne -999 && RH ne -999 && pres ne -999 && HaWS ne -999)
unstable =  and(theta < Tsurf_h , WS_h > WS_lim); % && T ne -999 && Tsurf_h ne -999 && RH ne -999 && pres ne -999 && HaWS ne -999)

% stable   =  (theta >= Tsurf_h && WS > WS_lim && T ne -999 && Tsurf_h ne -999 && RH ne -999 && pres ne -999 && HaWS ne -999)
% unstable =  (theta < Tsurf_h && WS > WS_lim && T ne -999 && Tsurf_h ne -999 && RH ne -999 && pres ne -999 && HaWS ne -999)
%no_wind  =  ( WS ne -999    && WS <= WS_lim && T ne -999 && Tsurf_h ne -999 && RH ne -999 && pres ne -999 && HaWS ne -999)

for i=1:30   % stable stratification
  psi_m1 = -(aa.*         z_0./L(stable) + bb.*(         z_0./L(stable)-cc./dd).*exp(-dd.*         z_0./L(stable)) + bb.*cc./dd);
  psi_m2 = -(aa.*z_WS(stable)./L(stable) + bb.*(z_WS(stable)./L(stable)-cc./dd).*exp(-dd.*z_WS(stable)./L(stable)) + bb.*cc./dd);
  psi_h1 = -(aa.*z_0h(stable)./L(stable) + bb.*(z_0h(stable)./L(stable)-cc./dd).*exp(-dd.*z_0h(stable)./L(stable)) + bb.*cc./dd);
  psi_h2 = -(aa.* z_T(stable)./L(stable) + bb.*( z_T(stable)./L(stable)-cc./dd).*exp(-dd.* z_T(stable)./L(stable)) + bb.*cc./dd);
  u_star(stable) = kappa.*WS_h(stable)./(log(z_WS(stable)./z_0)-psi_m2+psi_m1);
  Re(stable) = u_star(stable).*z_0./nu_h(stable);
  z_0h(stable) = z_0.*exp(1.5-0.2.*log(Re(stable))-0.11.*(log(Re(stable))).^2);
  if sum( (z_0h(stable) < 1e-6)) > 1 
      z_0h(stable( (z_0h(stable) < 1e-6))) = 1e-6;
  end
  th_star = kappa.*(theta(stable)-Tsurf_h(stable))./(log(z_T(stable)./z_0h(stable))-psi_h2+psi_h1);
  q_star  = kappa.*(  q_h(stable)- q_surf(stable))./(log(z_T(stable)./z_0h(stable))-psi_h2+psi_h1);
  SHF_h(stable) = rho_atm_h(stable).*c_pd .*u_star(stable).*th_star;
  LHF_h(stable) = rho_atm_h(stable).*L_sub.*u_star(stable).* q_star;
  L_prev = L(stable);
  L(stable) = u_star(stable).^2.*(theta(stable))...
      .*(1+((1-eps)./eps).*q_h(stable))./(g.*kappa.*th_star.*(1+((1-eps)./eps).*q_star));
  L_dif = abs((L_prev-L(stable))./L_prev);
%  print,"HF iterations stable stratification: ",i+1,sum( (L_dif > L_dif_max)),100..*sum( (L_dif > L_dif_max))/sum( (L_dif))
  if sum( (L_dif > L_dif_max)) == 1 
      break
  end
end

if sum(unstable) > 1  
  for i=1:30    % unstable stratification
    x1  = (1-gamma.*z_0           ./L(unstable)).^0.25;
    x2  = (1-gamma.*z_WS(unstable)./L(unstable)).^0.25;
    y1  = (1-gamma.*z_0h(unstable)./L(unstable)).^0.5;
    y2  = (1-gamma.*z_T(unstable) ./L(unstable)).^0.5;
    psi_m1 = log(((1+x1)./2).^2.*(1+x1.^2)./2)-2.*atan(x1)+pi./2;
    psi_m2 = log(((1+x2)./2).^2.*(1+x2.^2)./2)-2.*atan(x2)+pi./2;
    psi_h1 = log(((1+y1)./2).^2);
    psi_h2 = log(((1+y2)./2).^2);
    u_star(unstable) = kappa.*WS_h(unstable)./(log(z_WS(unstable)./z_0)-psi_m2+psi_m1);
    Re(unstable) = u_star(unstable).*z_0./nu_h(unstable);
    z_0h(unstable) = z_0.*exp(1.5-0.2.*log(Re(unstable))-0.11.*(log(Re(unstable))).^2);
    if sum( (z_0h(unstable) < 1e-6)) > 1 
        z_0h(unstable( (z_0h(unstable) < 1e-6))) = 1e-6;
    end
    th_star = kappa.*(theta(unstable)-Tsurf_h(unstable))./(log(z_T(unstable)./z_0h(unstable))-psi_h2+psi_h1);
    q_star  = kappa.*(  q_h(unstable)- q_surf(unstable))./(log(z_T(unstable)./z_0h(unstable))-psi_h2+psi_h1);
    SHF_h(unstable) = rho_atm_h(unstable).*c_pd .*u_star(unstable).*th_star;
    LHF_h(unstable) = rho_atm_h(unstable).*L_sub.*u_star(unstable).* q_star;
    L_prev = L(unstable);
    L(unstable) = u_star(unstable).^2.*(theta(unstable)).*(1+((1-eps)./eps).*q_h(unstable))./(g.*kappa.*th_star.*(1+((1-eps)./eps).*q_star));
    L_dif = abs((L_prev-L(unstable))./L_prev);
%    print,"HF iterations unstable stratification: ",i+1,sum( (L_dif > L_dif_max)),100..*sum( (L_dif > L_dif_max))./sum( (L_dif))
    if sum( (L_dif > L_dif_max)) == 1 
        break
    end
  end
end

q_h = 1000.*q_h; % from kg./kg to g./kg
% no_q =  (pres == -999 or T == -999 or RH == -999)
% if total(no_q) ne -1 then q_h(no_q) = -999
% no_HF =  (pres == -999 or T == -999 or Tsurf_h == -999 or RH == -999 or WS == -999 or HaWS == -999)
% if total(no_HF) ne -1  
%   SHF_h(no_HF) = -999
%   LHF_h(no_HF) = -999
%   end
SHF2 = SHF_h;
%----------------------------------------------------------------------------------------------------
%%
figure
hold on
plot(data.SensibleHeatFluxWm2(500:3000))
% plot(SHF1(500:1000))
plot(SHF2(500:3000))
ylabel('SHF')
legend('from file', 'calculated')


time = datetime(data.Year,data.MonthOfYear,data.DayOfMonth,data.HourOfDayUTC,0,0);

figure
plot(time,L)
hold on
plot(data_AWS.time,L_mod)
%%
figure
histogram(L,-500:10:1000)
hold on
histogram(L_mod,-500:10:1000)
xlim([-1000 1000])
