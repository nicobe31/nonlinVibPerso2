function [IC, hvalue] = shootingNLFR(period,X0,M,K,C,a,Ampl,phaseF)
%%
% [IC, hvalue] = shootingNLFR(period,X0,M,K,C,a,Ampl,phaseF)
%
% Shooting method to determine the initial conditions for a periodic
% soltuion with the given excitation
%
% period, the period of integration, 1/frequency
% X0 initial guess for the initial conditions
% M,K,C the matrices of the linear system
% a the nonlinear coefficients
% Ampl and phaseF: amplitude and phase of the force
%
% IC the initial conditions satsfyng the periodicity condition
% hvalue periodicity condition for the solution
%
%%

unit2fminunc = 1000; % pass from m to mm
x0 = X0*unit2fminunc;

optionsOpti = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','TolX',1e-5,'TolFun',1e-5);
[dxOpti, hvalue, ~,~] = fminunc(@(dx) h2minNLFR(M,K,C,a,dx,period,unit2fminunc,x0,Ampl,phaseF),[0;0;0],optionsOpti);

dxOpti = [dxOpti(1:2); 0; dxOpti(3)];
dxOpti = dxOpti/unit2fminunc;
signX0 = sign(X0);
signX0(signX0==0) = 1;
IC = X0 + signX0.*dxOpti;

end
