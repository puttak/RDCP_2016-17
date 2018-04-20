function [sys,x0,str,ts]= HeatExchanger_sfcn(t,x,u,flag, Th1, Th2, Th3, Th4, Th5,    Tc1, Tc2, Tc3, Tc4, Tc5)

switch flag
 case 0 % initialize
 str=[] ;
 ts = [0 0] ;
 s = simsizes ;
 s.NumContStates = 10 ;
 s.NumDiscStates = 0 ;
 s.NumOutputs = 10 ;
 s.NumInputs = 9 ;
 s.DirFeedthrough = 0 ;
 s.NumSampleTimes = 1 ;
 sys = simsizes(s) ;
 x0 = [  Th1, Th2, Th3, Th4, Th5,    Tc1, Tc2, Tc3, Tc4, Tc5  ] ;
 case 1 % derivatives
 u = u ;
 sys = SimulinkHeatExchanger(t,x,u) ;
 case 3 % output
 sys = x;
 case {2 4 9} % 2:discrete
 % 4:calcTimeHit
 % 9:termination
 sys =[];
 otherwise
 error(['unhandled flag =',num2str(flag)]) ;
 
end
end