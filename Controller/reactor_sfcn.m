function [sys,x0,str,ts]= reactor_sfcn(t,x,u,  flag,    Wm1, Wm2, Wm3, Wm4, Wm5, Wm6, Wm7, Wm8, Wm9, Wm10,    Wme1, Wme2, Wme3, Wme4, Wme5, Wme6, Wme7, Wme8, Wme9, Wme10,     Ww1, Ww2, Ww3, Ww4, Ww5, Ww6, Ww7, Ww8, Ww9, Ww10,    Wc1, Wc2, Wc3, Wc4, Wc5, Wc6, Wc7, Wc8, Wc9, Wc10,    Wh1, Wh2, Wh3, Wh4, Wh5, Wh6, Wh7, Wh8, Wh9, Wh10,    T1, T2, T3, T4, T5, T6, T7, T8, T9, T10)

n=10;
switch flag
 case 0 % initialize
 str=[] ;
 ts = [0 0] ;
 s = simsizes ;
 s.NumContStates = 6*n ;
 s.NumDiscStates = 0 ;
 s.NumOutputs = 60 ;
 s.NumInputs = 3 ;
 s.DirFeedthrough = 0 ;
 s.NumSampleTimes = 1 ;
 sys = simsizes(s) ;
 x0 = [Wm1 Wm2 Wm3 Wm4 Wm5 Wm6 Wm7 Wm8 Wm9 Wm10,    Wme1 Wme2 Wme3 Wme4 Wme5 Wme6 Wme7 Wme8 Wme9 Wme10,     Ww1 Ww2 Ww3 Ww4 Ww5 Ww6 Ww7 Ww8 Ww9 Ww10,    Wc1 Wc2 Wc3 Wc4 Wc5 Wc6 Wc7 Wc8 Wc9 Wc10,    Wh1 Wh2 Wh3 Wh4 Wh5 Wh6 Wh7 Wh8 Wh9 Wh10,    T1 T2 T3 T4 T5 T6 T7 T8 T9 T10] ;
 
 
 case 1 % derivatives
 u = u ;
 [dx] = SimulinkReactor(t,x,u) ;

 sys = dx;
 
 %  sys= SimulinkReactor(t,x,u) ;
 
 case 3 % output
  sys = x;
%   for i=1:n
%     Wm(i) = x(i);
%     Wme(i) = x(i+n);
%     Ww(i) = x(i+2*n);
%     Wc(i) = x(i+3*n);
%     Wh(i) = x(i+4*n);
%     T(i) = x(i+5*n);
%   end
%    P0=300;
%    P = pressure(Wm,Wme,Ww,Wc,Wh,T,P0);
% % %   
%   sys = [  Wm(end) ; Wme(end) ; Ww(end) ; Wc(end) ; Wh(end) ; T(end); P(end)  ];% 
% %  ; P(end)
 
 case {2 4 9} % 2:discrete
 % 4:calcTimeHit
 % 9:termination
 sys =[];
 otherwise
 error(['unhandled flag =',num2str(flag)]) ;
 
end
end

% initial conditions is 0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,    600,600,600,600,600,600,600,600,600,600