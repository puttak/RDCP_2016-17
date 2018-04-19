function methanol_reactor_design_Group_1()
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line 15 : Specification of constant simulation parameters
% Line 35 : Initiation of feed conditions
% Line 67 : Calculation prior ode solver
% Line 77 : ode calling
% Line 81 : Calculation after ode solver
% Line 107: Graph plotting
% Line 186: ode
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
disp('Simulation started.')
global R Area Cpcon K1con R1con R2con deltaHr Catalystcon

%%Specification of constant simulation parameters
R    = 8.3144598;    %%[kJ kmol^-1 K^-1]
Area = 1;         %%[m^2]
Cpcon=[19.037   9.146    -1.217    -8.033;
       28.142   0.167     0.537    -2.220;
       29.087  -0.191     0.400    -0.870;
       19.874   5.021     1.268    -11.00;
       32.217   0.192     1.055    -3.593];
   %[From top to bottom: Methanol CO H2 Methane H2O]
   %[From left to right: a bx10^2 cx10^5 dx10^9]
   %[kJ kmol^-1]    
K1con = [5332 12.831];    %%[(unitless)]
R1con = [-1103.6    -308.4     171.6    -14.64;
           92.65    -17.95   -0.5265    0.1745;
          34.295    -5.444   -0.6771    0.1088;
           483.5    -28.21     -22.8     2.438];
R2con = 0.5*10^15;
deltaHr = [-89980 -206296];              %[J mol^-1]
Catalystcon = [5.5 16 7.04 1400 0.4];    %[mm mm mm kgm^-3 (unitless)]

%%Initiation of feed conditions
x = input('Use Group 1''s feed conditions? (Yes=1, No=0)');
if x==1
   FM0=6.3;      %%[kmol s^-1]
   CCO0=0.2;     %%[mole fraction]
   CH20=0.75;    %%[mole fraction]
   CMeOH0=0.01;  %%[mole fraction]
   CMet0=0.03;   %%[mole fraction]
   CH2O0=0.01;   %%[mole fraction]
   T0=630;       %%[K]
   P0=300;       %%[bar]
   L=7;          %%[m] %minimum length is 4.108m, extended to 7 to see full trend
elseif x==0 
   FM0   =input('Please enter feed flowrate in kmol/s:');
   CCO0  =input('Please enter mole fraction of CO in the feed:');
   CH20  =input('Please enter mole fraction of H2 in the feed:');
   CMeOH0=input('Please enter mole fraction of methanol in the feed:');
   CMet0 =input('Please enter mole fraction of methane in the feed:');
   CH2O0 =input('Please enter mole fraction of water in the feed:');
   T0    =input('Please enter feed temperature in K:');
   P0    =input('Please enter feed pressure in bar:');
   L     =input('Please enter length of reactor to be simulated in m:');
   a     =input('End of input prompt. Proceed with simulations? (Yes=1, No=0)');
   if a==0
       error('Valid input. Simulation aborted.')
   elseif a~=1
       error('Invalid input. Only allowed inputs are 0 and 1. Simulation aborted.')
   end
else
   error('Invalid input. Only allowed inputs are 0 and 1. Simulation aborted.')
end

% Calculation prior ode solver
nCO0=CCO0*FM0;       %%[kmol s^-1]
nH20=CH20*FM0;       %%[kmol s^-1]
nMeOH0=CMeOH0*FM0;   %%[kmol s^-1]
nMet0=CMet0*FM0;     %%[kmol s^-1]
nH2O0=CH2O0*FM0;     %%[kmol s^-1]
odelength=linspace(0,L,10000);
m0=[T0 P0 nCO0 nH20 nMeOH0 nMet0 nH2O0];
disp('Solving ode...')

% ode caller
[l,m]=ode45(@RCDP_ode,odelength, m0);
disp('Ode solved. Parsing and processing data...')

%Calculation after ode solver
W=Catalystcon(4)*Area*(1-Catalystcon(5))*l;
[Length,width]=size(m);
Temperature = m(:,1);
Pressure    = m(:,2);
xi1         = m(:,5)-nMeOH0;
xi2         = m(:,6)-nMet0;

for i=1:Length-1                              %Find expressions between 2 sample lengths
r1(i)=(xi1(i+1)-xi1(i))./(l(i+1)-l(i));       %Calculate rate 1 at each length intervals
r2(i)=(xi2(i+1)-xi2(i))./(l(i+1)-l(i));       %Calculate rate 2 at each length intervals
dl(i)=(l(i+1)+l(i))/2;                        %Average length at each length intervals
avexi1(i)=(xi1(i+1)+xi1(i))./2;               %Average xi1 at each length intervals
avexi2(i)=(xi2(i+1)+xi2(i))./2;               %Average xi2 at each length intervals
aveT(i)=(Temperature(i+1)+Temperature(i))./2; %Average temperature at each length intervals
end

length=interp1(xi1, l, 0.10837);
I=num2str(length);
if isnan(length)
    str=('Required methanol flowrate specification is not reached. Rerun simulation and change input.');
else
    str=strcat('Minimum length to achieve the required methanol flowrate specification is ',I,'m');
end
disp (str)
disp('Data parsed and processed. Plotting graphs...')

%Graph plotting
figure('Name','Required graphs','units','normalized','outerposition',[0 0 1 1])
subplot(3,2,1)
axis([0,0.12,0,7000])
plot(xi1,W);
title('Catalyst weight against \xi1');
xlabel('\xi1 (kmol/s)');
ylabel('Catalyst weight (kg)');

subplot(3,2,2)
axis([0,0.12,0,0.04])
plot (xi1, xi2);
title('\xi2 against \xi1')
xlabel('\xi1 (kmol/s)');
ylabel('\xi2 (kmol/s)');

subplot(3,2,3)
axis([0,0.12,500,800])
plot (xi1,Temperature);
title('Temperature against \xi1');
xlabel('\xi1 (kmol/s)');
ylabel('Temperature (K)');

subplot(3,2,4)
axis([0,0.12,0,500])
plot (xi1,Pressure);
title('Pressure against \xi1');
xlabel('\xi1 (kmol/s)');
ylabel('Pressure (bar)');

subplot(3,2,5)
plot(avexi1,r1); %plotted average extend of reaction 1 at each length interval against ROR 1 
title('Rate of reaction 1 (per unit length) against \xi1')
xlabel('\xi1 (kmol/s)')
ylabel('ROR1 (kmol/m.s)')




figure('Name','Extra graphs','units','normalized','outerposition',[0 0 1 1])
subplot(3,2,1)
plot(dl,r1);
title('Rate of reaction 1 (per unit length) against Length')
xlabel('Length (m)')
ylabel('ROR1 (kmol/m.s)')

subplot(3,2,2)
plot(dl,r2);
title('Rate of reaction 2 (per unit length) against Length')
xlabel('Length (m)')
ylabel('ROR2 (kmol/m.s)')

subplot(3,2,3)
axis([0,L,0,0.12]);
plot (l,xi1,l,xi2);
title('\xi1 and \xi2 against Length')
legend('\xi1','\xi2');
xlabel('Length (m)');
ylabel('\xi1 and \xi2 (kmol/s)');

subplot(3,2,4)
plot(l,Temperature);
axis([0,L,600,750]);
title('Temperature against Length')
xlabel('Length (m)');
ylabel('Temperature (K)');

subplot(3,2,5)
plot(aveT,r1);
title('Rate of reaction 1 (per unit length) against Temperature')
xlabel('Temperature (K)')
ylabel('ROR1 (kmol/m.s)')

subplot(3,2,6)
plot(avexi1,r2);
title('Rate of reaction 2 (per unit length) against \xi1')
xlabel('\xi1 (kmol/s)')
ylabel('ROR2 (kmol/m.s)')

disp('Run completed.')
end


function dm = RCDP_ode(~, m)

global Area Cpcon K1con R1con R2con deltaHr Catalystcon R

%m = [T P nCO nH2 nMeOH nMet nH2O];
dm = zeros(7,1);    % a column vector
nCO=m(3);   %%[kmol s^-1]
nH2=m(4);   %%[kmol s^-1]
nMeOH=m(5); %%[kmol s^-1]
nMet=m(6);  %%[kmol s^-1]
nH2O=m(7);  %%[kmol s^-1]
n = [nCO nH2 nMeOH nMet nH2O];
ntot=n(1)+n(2)+n(3)+n(4)+n(5);


%catalyst weight
W=Catalystcon(4)*Area*(1-Catalystcon(5));

%Mass flowrate calculation
%define molecular masses
M = [28 2 32 16 18];
%ma=[mCO mH2 mMeOH mMet mH2O];
ma = n.*M;
mtot=ma(1)+ma(2)+ma(3)+ma(4)+ma(5); %%[kg s^-1]

%Volumetric flowrate calculation using ideal gas law
%molar flowrate vector
%n = [nCO nH2 nMeOH nMet nH2O];
%v=[vCO vH2 vMeOH vMet vH2O]
v = (n*1000)*R*m(1)/(m(2)*10^5);
vtot=v(1)+v(2)+ v(3)+ v(4)+ v(5);   %%[m^3]

%Partial Pressure calculation
%n = [nCO nH2 nMeOH nMet nH2O];
%Pa=[PCO PH2 PMeOH]
Pa = m(2)*(n/(ntot)); %%[bar]

%Cp calculation
Cpcon(:,2)=Cpcon(:,2).*(10^-2);
Cpcon(:,3)=Cpcon(:,3).*(10^-5);
Cpcon(:,4)=Cpcon(:,4).*(10^-9);
CpMeOH=(Cpcon(1,1)+Cpcon(1,2).*m(1)+Cpcon(1,3).*m(1).^2+Cpcon(1,4).*m(1).^3);
CpCO=(Cpcon(2,1)+Cpcon(2,2).*m(1)+Cpcon(2,3).*m(1).^2+Cpcon(2,4).*m(1).^3);
CpH2=(Cpcon(3,1)+Cpcon(3,2).*m(1)+Cpcon(3,3).*m(1).^2+Cpcon(3,4).*m(1).^3);
CpMet=(Cpcon(4,1)+Cpcon(4,2).*m(1)+Cpcon(4,3).*m(1).^2+Cpcon(4,4).*m(1).^3);
CpH2O=(Cpcon(5,1)+Cpcon(5,2).*m(1)+Cpcon(5,3).*m(1).^2+Cpcon(5,4).*m(1).^3);   %%[J K^-1 mol^-1]
Cptot=((nMeOH*CpMeOH+nCO*CpCO+nH2*CpH2+nMet*CpMet+nH2O*CpH2O).*1000);          %%[J K^-1 s^-1]

%Reaction 1 calculations
K=10^(K1con(1)/m(1)-K1con(2));    %[bar^-2]
A=(R1con(1,1)+R1con(1,2).*(m(1)/100)+R1con(1,3).*(m(1)/100)^2+R1con(1,4).*(m(1)/100)^3);
B=(R1con(2,1)+R1con(2,2).*(m(1)/100)+R1con(2,3).*(m(1)/100)^2+R1con(2,4).*(m(1)/100)^3);
C=(R1con(3,1)+R1con(3,2).*(m(1)/100)+R1con(3,3).*(m(1)/100)^2+R1con(3,4).*(m(1)/100)^3);
D=(R1con(4,1)+R1con(4,2).*(m(1)/100)+R1con(4,3).*(m(1)/100)^2+R1con(4,4).*(m(1)/100)^3);
r1=((Pa(1)).*((Pa(2)).^2)-((Pa(3))./K)).*W/(((A+B.*(Pa(1))+C.*(Pa(2))+D.*(Pa(3))).^3).*60.*60); %%[kmol s^-1 m^-1]

%Reaction 2 calculation
r2=R2con.*exp(-30000./m(1)).*(Pa(1)).*W./(60*60);   %%[kmol s^-1 m^-1]

%Pressure drop calculation
G=mtot./Area;    %%[kg m^-2 s^-1]
ro=mtot/vtot;    %%[kg m^-3]

%Cp calculations
funMeOH=@(T) (Cpcon(1,1).*T+Cpcon(1,2).*(T.^2)./2+Cpcon(1,3).*(T.^3)./3+Cpcon(1,4).*(T.^4)./4);
funCO=@(T) (Cpcon(2,1).*T+Cpcon(2,2).*(T.^2)./2+Cpcon(2,3).*(T.^3)./3+Cpcon(2,4).*(T.^4)./4);
funH2=@(T) (Cpcon(3,1).*T+Cpcon(3,2).*(T.^2)./2+Cpcon(3,3).*(T.^3)./3+Cpcon(3,4).*(T.^4)./4);
funMet=@(T) (Cpcon(4,1).*T+Cpcon(4,2).*(T.^2)./2+Cpcon(4,3).*(T.^3)./3+Cpcon(4,4).*(T.^4)./4);
funH2O=@(T) (Cpcon(5,1).*T+Cpcon(5,2).*(T.^2)./2+Cpcon(5,3).*(T.^3)./3+Cpcon(5,4).*(T.^4)./4);
dHr1=(funCO(298)-funCO(m(1)))+2*(funH2(298)-funH2(m(1)))+deltaHr(1)+(funMeOH(m(1))-funMeOH(298));
dHr2=(funCO(298)-funCO(m(1)))+3*(funH2(298)-funH2(m(1)))+deltaHr(2)+(funMet(m(1))-funMet(298))+(funH2O(m(1))-funH2O(298));

%mass, energy and pressure derivative eqns
%dT/dl
dm(1)=-(dHr1*1000*r1+dHr2*1000*r2)/(Cptot);
%dP/dl
dm(2)=-(1.75*(G^2)*(1-Catalystcon(5))/((Catalystcon(3)*10^-3)*(Catalystcon(5)^3)*(ro)))*(1e-5);
%convert all rates to dl
%dC/dl
dm(3)=-r1-r2;
dm(4)=-2*r1-3*r2;
dm(5)=r1;                    
dm(6)=r2;
dm(7)=r2;
end
