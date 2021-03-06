% m =methanol
% me=methane
% w =water
% c =carbon monoxide
% h =hydrogen
function [dx,Pout] = SimulinkReactor2(t,x,u) % where t is time (h), x is a vector of state variales, u is a vector of inputs

% molecular weights
global Mrm Mrme Mrw Mrc Mrh
Mrm = 32;
Mrme = 16;
Mrw = 18;
Mrc = 28;
Mrh = 2;

% input information
global P0
Wm0 = u(1);
Wme0 = u(2);
Ww0 = u(3);
Wc0 = u(4);
Wh0 = u(5);
P0 = u(6);
T0 = u(7);

% G=mass flux (kg/hw ); Area=reactor x-section area (m^2); Dp=catalyst size (m);
% V=reactor volume (m^3); E=voidage; Cpcat=catalyst Cp (kJ/(kg K)); M=total mass
% flowrate (kg/hr); H1=standard heat of reaction (kJ/kmol); n=number of discretised 
% steps
global G Area Dp V E Cpcat rhoc M R H1 H2 n
n=10;
M=u(8);
Area = 1;
G = M/Area;
Dp = 7.04*10^-3;
V =4.4/n;
E = 0.4;
Cpcat = 0.85;
rhoc = 1400;
R=8.314;
H1= -81000; 
H2=-185400; 
global ac bc cc dc
global am bm cm dm
global ah bh ch dh
global ame bme cme dme
global aw bw cw dw
am = 19.037; bm = 9.146.*10.^-2; cm = -1.217.*10.^-5; dm = -8.033.*10.^-9; 
ac = 28.142; bc = 0.167.*10.^-2; cc = 0.537.*10.^-5; dc = -2.22.*10.^-9; 
ah = 29.087; bh = -0.191.*10.^-2; ch = 0.4.*10.^-5; dh = -0.87.*10.^-9; 
ame = 19.874; bme = 5.021.*10.^-2; cme = 1.268.*10.^-5;dme = -11.*10.^-9; 
aw = 32.217; bw = 0.192.*10.^-2; cw = 1.055.*10.^-5; dw = -3.593.*10.^-9; 

for i=1:n
    Wm(i) = x(i);
    Wme(i) = x(i+n);
    Ww(i) = x(i+2*n);
    Wc(i) = x(i+3*n);
    Wh(i) = x(i+4*n);
    T(i) = x(i+5*n);
end

[Hr1] = Heat_reaction1 (T);
[Hr2] = Heat_reaction2 (T);
r1 = rate1 (Wm,Wme,Ww,Wc,Wh,T);
r2 = rate2 (Wm,Wme,Ww,Wc,Wh,T);
CpA = CpAverage (Wm,Wme,Ww,Wc,Wh,T);
rhog = GasDensity (Wm,Wme,Ww,Wc,Wh,T);
P = pressure (Wm,Wme,Ww,Wc,Wh,T);
Pout=P(end);

% mass and energy balance for all discretized volumes of reactor
for i=1:n
    if i==1
dWmdt(i) = -(M./(E.*rhog(i))).*(Wm(i)-Wm0)./V + Mrm.*r1(i).*rhoc.*(1-E)./(E.*rhog(i));
dWmedt(i) = -(M./(E.*rhog(i))).*(Wme(i)-Wme0)./V + Mrme.*r2(i).*rhoc.*(1-E)./(E.*rhog(i));
dWwdt(i) = -(M./(E.*rhog(i))).*(Ww(i)-Ww0)./V + Mrw.*r2(i).*rhoc.*(1-E)./(E.*rhog(i));
dWcdt(i) = -(M./(E.*rhog(i))).*(Wc(i)-Wc0)./V + Mrc.*(-r1(i)-r2(i)).*rhoc.*(1-E)./(E.*rhog(i));
dWhdt(i) = -(M./(E.*rhog(i))).*(Wh(i)-Wh0)./V + Mrh.*(-2.*r1(i)-3.*r2(i)).*rhoc.*(1-E)./(E.*rhog(i));
dTdt(i) = (-M.*CpA(i).*(T(i)-T0)./V + (r1(i).*(-Hr1(i)) + r2(i).*(-Hr2(i))).*rhoc.*(1-E))./(E.*rhog(i).*CpA(i) + (1-E).*rhoc.*Cpcat);
    else
dWmdt(i) = -(M./(E.*rhog(i))).*(Wm(i)-Wm(i-1))./V + Mrm.*r1(i).*rhoc.*(1-E)./(E.*rhog(i));
dWmedt(i) = -(M./(E.*rhog(i))).*(Wme(i)-Wme(i-1))./V + Mrme.*r2(i).*rhoc.*(1-E)./(E.*rhog(i));
dWwdt(i) = -(M./(E.*rhog(i))).*(Ww(i)-Ww(i-1))./V + Mrw.*r2(i).*rhoc.*(1-E)./(E.*rhog(i));
dWcdt(i) = -(M./(E.*rhog(i))).*(Wc(i)-Wc(i-1))./V + Mrc.*(-r1(i)-r2(i)).*rhoc.*(1-E)./(E.*rhog(i));
dWhdt(i) = -(M./(E.*rhog(i))).*(Wh(i)-Wh(i-1))./V + Mrh.*(-2.*r1(i)-3.*r2(i)).*rhoc.*(1-E)./(E.*rhog(i));
dTdt(i) = (-M.*CpA(i).*(T(i)-T(i-1))./V + (r1(i).*(-Hr1(i)) + r2(i).*(-Hr2(i))).*rhoc.*(1-E))./(E.*rhog(i).*CpA(i) + (1-E).*rhoc.*Cpcat);
end
end

% dx is a vector of ODEs
dx = [dWmdt'; dWmedt'; dWwdt'; dWcdt'; dWhdt'; dTdt'];
end

%% mole fraction
function y = MoleFraction (Wm,Wme,Ww,Wc,Wh) 
global Mrm Mrme Mrw Mrc Mrh
y.m = (Wm./Mrm)./(Wm./Mrm + Wme./Mrme + Ww./Mrw + Wc./Mrc + Wh./Mrh);
y.me = (Wme./Mrme)./(Wm./Mrm + Wme./Mrme + Ww./Mrw + Wc./Mrc + Wh./Mrh);
y.w = (Ww./Mrw)./(Wm./Mrm + Wme./Mrme + Ww./Mrw + Wc./Mrc + Wh./Mrh);
y.c = (Wc./Mrc)./(Wm./Mrm + Wme./Mrme + Ww./Mrw + Wc./Mrc + Wh./Mrh);
y.h = (Wh./Mrh)./(Wm./Mrm + Wme./Mrme + Ww./Mrw + Wc./Mrc + Wh./Mrh);
end

%% partial pressure (bar)
function PP = PartialPressure (Wm,Wme,Ww,Wc,Wh,T) 
y = MoleFraction (Wm,Wme,Ww,Wc,Wh);
P = pressure (Wm,Wme,Ww,Wc,Wh,T);
PP.c = P.*y.c;
PP.h = P.*y.h;
PP.m = P.*y.m;
end

%% equilibrium constant (bar.^(-2))
function [K] = equil_constant (T) 
K = 10.^((5332./T)-12.831); 
end

%% factors for rate equation
function [f] = factors (T) 
f.A=-1103.6-308.4.*(T./100)+171.6.*(T./100).^2-14.64.*(T./100).^3;
f.B=92.65-17.95.*(T./100)-0.5265.*(T./100).^2+0.1745.*(T./100).^3;
f.C=34.295-5.444.*(T./100)-0.6771.*(T./100).^2+0.1088.*(T./100).^3;
f.D=483.5-28.21.*(T./100)-22.8.*(T./100).^2+2.438.*(T./100).^3;
end

%% rate of reaction 1 (kmol/kgcat hr)
function r1 = rate1 (Wm,Wme,Ww,Wc,Wh,T) 
PP = PartialPressure (Wm,Wme,Ww,Wc,Wh,T);
[f] = factors (T);
[K] = equil_constant (T);
r1 = (((PP.c.*(PP.h.^2))-(PP.m./K))./((f.A+f.B.*PP.c+f.C.*PP.h+f.D.*PP.m).^3));
end

%% rate of reaction 2 (kmol/kgcat hr)
function r2 = rate2 (Wm,Wme,Ww,Wc,Wh,T)
PP = PartialPressure (Wm,Wme,Ww,Wc,Wh,T);
r2 =  0.5.*10.^15.*exp(-30000./T).*PP.c;
end

%% gas density (kg/m^3)
function rhog = GasDensity (Wm,Wme,Ww,Wc,Wh,T)
global R
P = pressure (Wm,Wme,Ww,Wc,Wh,T);
MrA = MrAverage (Wm,Wme,Ww,Wc,Wh);
rhog = (P.*(1E5)./(R.*1000)./T).*MrA;
end

%% total pressure (bar)
function P = pressure(Wm,Wme,Ww,Wc,Wh,T)
global G Area Dp V P0 E R n
MrA = MrAverage (Wm,Wme,Ww,Wc,Wh);
polynomial(1,:) = [1 -P0.*100000 1.75.*(G./3600).^2.*(1-E).*V.*R.*T(1)./(Area.*Dp.*(MrA(1)./1000).*E.^3) ];
RootVector(1,:) = roots(polynomial(1,:));
P(1) = max(RootVector(1,:))./100000;
for i=2:n
polynomial(i,:) = [1 -P(i-1).*100000 1.75.*(G./3600).^2.*(1-E).*V.*R.*T(i)./(Area.*Dp.*(MrA(i)./1000).*E.^3) ];
RootVector(i,:) = roots(polynomial(i,:));
P(i) = max(RootVector(i,:))./100000;
end
end

%% molar average molecular weight (kg/kmol)
function MrA = MrAverage (Wm,Wme,Ww,Wc,Wh) 
global Mrm Mrme Mrw Mrc Mrh
y = MoleFraction (Wm,Wme,Ww,Wc,Wh);
MrA = Mrc.*y.c + Mrh.*y.h + Mrm.*y.m + Mrme.*y.me + Mrw.*y.w;
end
 
%% average mass Cp (kJ/(K kg))
function CpA = CpAverage (Wm,Wme,Ww,Wc,Wh,T) 
CpMolarA = CpMolarAverage (Wm,Wme,Ww,Wc,Wh,T);
MrA = MrAverage (Wm,Wme,Ww,Wc,Wh);
CpA = CpMolarA./MrA;
end

%% average molar Cp (kJ/(K kmol))
function CpMolarA = CpMolarAverage (Wm,Wme,Ww,Wc,Wh,T) 
y = MoleFraction (Wm,Wme,Ww,Wc,Wh);
CpMolar_m = CpMolar_Methanol (T);
CpMolar_c = CpMolar_CO (T);
CpMolar_h = CpMolar_Hydrogen (T);
CpMolar_me = CpMolar_Methane (T);
CpMolar_w = CpMolar_Water (T);
CpMolarA = CpMolar_c.*y.c + CpMolar_h.*y.h + CpMolar_m.*y.m + CpMolar_me.*y.me + CpMolar_w.*y.w;
end

%% molar Cp of different chemicals (kJ/(K kmol))
function CpMolar_m = CpMolar_Methanol (T)
    global am bm cm dm
            CpMolar_m = (am + bm.*T + cm.*T.^2 + dm.*T.^3);
 end
function CpMolar_c = CpMolar_CO (T)
    global ac bc cc dc
        CpMolar_c = (ac + bc.*T + cc.*T.^2 + dc.*T.^3);
end
function CpMolar_h = CpMolar_Hydrogen (T)
    global ah bh ch dh
        CpMolar_h = (ah + bh.*T + ch.*T.^2 + dh.*T.^3);
end
function CpMolar_me = CpMolar_Methane (T)
    global ame bme cme dme
        CpMolar_me = (ame + bme.*T + cme.*T.^2 + dme.*T.^3);
end
function CpMolar_w = CpMolar_Water (T)
    global aw bw cw dw
    
        CpMolar_w = (aw + bw.*T + cw.*T.^2 + dw.*T.^3);
end

%% ingegral of molar Cp for different chemicals
function Int_CpMolar_m = Integrated_CpMolar_Methanol (T)
    global am bm cm dm
            Int_CpMolar_m = (am.*T + bm.*T.^2./2 + cm.*T.^3./3 + dm.*T.^4./4);
 end
function Int_CpMolar_c = Integrated_CpMolar_CO (T)
    global ac bc cc dc
        Int_CpMolar_c = (ac.*T + bc.*T.^2./2 + cc.*T.^3./3 + dc.*T.^4./4);
end
function Int_CpMolar_h = Integrated_CpMolar_Hydrogen (T)
    global ah bh ch dh
        Int_CpMolar_h = (ah.*T + bh.*T.^2./2 + ch.*T.^3./3 + dh.*T.^4./4);
end
function Int_CpMolar_me = Integrated_CpMolar_Methane (T)
    global ame bme cme dme
        Int_CpMolar_me = (ame.*T + bme.*T.^2./2 + cme.*T.^3./3 + dme.*T.^4./4);
end
function Int_CpMolar_w = Integrated_CpMolar_Water (T)
    global aw bw cw dw
    
        Int_CpMolar_w = (aw.*T + bw.*T.^2./2 + cw.*T.^3./3 + dw.*T.^4./4);
end

%% heats of reaction (kJ/kmol)
function [Hr1] = Heat_reaction1 (T) 
global H1
Hr1 = H1 + (  Integrated_CpMolar_Methanol (T)-Integrated_CpMolar_Methanol (298)  ) - (  Integrated_CpMolar_CO (T)-Integrated_CpMolar_CO (298)  ) - 2.*(  Integrated_CpMolar_Hydrogen (T)-Integrated_CpMolar_Hydrogen (298)  );
end

function [Hr2] = Heat_reaction2 (T)
global H2
Hr2 = H2 + (  Integrated_CpMolar_Methane (T)-Integrated_CpMolar_Methane (298)  ) + (  Integrated_CpMolar_Water (T)-Integrated_CpMolar_Water (298)  ) - (  Integrated_CpMolar_CO (T)-Integrated_CpMolar_CO (298)  ) - 3.*(  Integrated_CpMolar_Hydrogen (T)-Integrated_CpMolar_Hydrogen (T)  );

end
