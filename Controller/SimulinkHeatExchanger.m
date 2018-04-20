function [dx]=SimulinkHeatExchanger (t,x,u)
global A L Circ U0 Cpw rhow AO n2  % Tin 

Fiwatertotal=u(1)/18; %Total cooling water flow rate (kmol/hr)
Fiwater=Fiwatertotal/50; %Cooling Water Flowrate per double pipe (splits into 50) (kmol/hr)
Wm=u(2); %Inlet weight fractions
Wme=u(3);
Ww=u(4);
Wc=u(5);
Wh=u(6);
Pin=u(7); %Inlet pressure (bar)
Tin = u(8); % temperature (K)
FiMass=u(9)/50; %Inlet process stream mass flowrate (kg/hr)

global FiCH3OH FiCH4 FiH2O FiCO FiH2 FiTotal
FiCH3OH=(FiMass)*Wm*(1/32); %Molar Flowrates (kmol/h)
FiCH4=(FiMass)*Wme*(1/16);
FiH2O=(FiMass)*Ww*(1/18);
FiCO=(FiMass)*Wc*(1/28);
FiH2=(FiMass)*Wh*(1/2);
FiTotal=FiCH3OH+FiCH4+FiH2O+FiCO+FiH2;%Total Molar Flowrate (kmol/h)

Cpw=75.3; %Cooling Water Heat Capacity (kJ/kmol.K)
rhow=54.848; %Cooling Water Density (kmol/m^3)

n2=5;

A=pi*0.1024.^2*(1/4); %Pipe Inner Area (m^2)
AO=(1/4)*pi()*(0.1542.^2-0.114.^2); %Outside Annulus Area (m^2)
L=3.8/n2; %Length of heat exchanger pipe divided by number of discrete sections it is split up into
Circ=pi*0.114; %Outside pipe circumference (m)
U0=(5.6783*3600/1000)*(0.0081*Fiwatertotal+5.7); %Heat Transfer Coefficient (kJ/K.m^2.hr)

for i=1:n2
    Th(i)=x(i);
    Tc(i)=x(n2+i);
end

for i=1:n2
    Cp=heatcapacities(Th(i)); 
    rhog=(Pin*10^5)/(8.314*Th(i)*1000); %Gas Density (kmol/m^3)
    if i==1
        DT=Th(i)-Tc(i+1);
        dThdt(i)=((FiTotal.*(Tin-Th(i)))./(rhog*A*L))-((Circ*U0*DT)./(rhog*A*Cp)); %Hot Stream Energy Balance 
        dTcdt(i)=((Fiwater.*(Tc(i+1)-Tc(i)))./(rhow*AO*L))+((Circ*U0*DT)./(rhow*AO*Cpw)); %Cooling Water Engergy Balance
    else
        if i==n2
        DT=Th(n2)-293;
        dThdt(i)=((FiTotal.*(Th(i-1)-Th(i)))./(rhog*A*L))-((Circ*U0*DT)./(rhog*A*Cp)); %Hot Stream Energy Balance
        dTcdt(i)=((Fiwater.*(293-Tc(i)))./(rhow*AO*L))+((Circ*U0*DT)./(rhow*AO*Cpw)); %Cooling Water Engergy Balance 
        else
        DT=Th(i)-Tc(i+1);
        dThdt(i)=((FiTotal.*(Th(i-1)-Th(i)))./(rhog*A*L))-((Circ*U0*DT)./(rhog*A*Cp)); %Hot Stream Energy Balance
        dTcdt(i)=((Fiwater.*(Tc(i+1)-Tc(i)))./(rhow*AO*L))+((Circ*U0*DT)./(rhow*AO*Cpw)); %Cooling Water Engergy Balance 
        end
        
    end 
end


[dx]=[  dThdt' ; dTcdt'  ];



end

function Cp=heatcapacities(Tout)
global FiCH3OH FiCH4 FiH2O FiCO FiH2 FiTotal
CpCH3OH=19.037+(9.146*10.^-2*Tout)-(1.217*10.^-5*Tout.^2)-(8.033*10.^-9*Tout.^3); %kJ/K.kmol
CpCH4=19.874+(5.021*10.^-2*Tout)+(1.268*10.^-5*Tout.^2)-(11*10.^-9*Tout.^3); %kJ/K.kmol
CpH2O=32.217+(0.192*10.^-2*Tout)+(1.055*10.^-5*Tout.^2)-(3.593*10.^-9*Tout.^3); %kJ/K.kmol
CpCO=28.142+(0.167*10.^-2*Tout)+(0.537*10.^-5*Tout.^2)-(2.22*10.^-9*Tout.^3); %kJ/K.kmol
CpH2=29.087-(0.191*10.^-2*Tout)+(0.4*10.^-5*Tout.^2)-(0.87*10.^-9*Tout.^3); %kJ/K.kmol
Cp=(CpCH3OH*FiCH3OH+CpCH4*FiCH4+CpH2O*FiH2O+CpCO*FiCO+CpH2*FiH2)/FiTotal; %Average Cp (kJ/kmol.K)
end

