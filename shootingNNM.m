function [IC, hvalue] = shootingNNM(period,X0,M,K,a)
%%
% [IC, hvalue] = shootingNNM(period,X0,M,K,a)
%
% Shooting method to determine the Nonlinear normal mode.
%
% period, the period of integration, 1/frequency
% X0 initial guess for the initial conditions
% M,K the matrices of the linear system
% a the nonlinear coefficients
%
% IC the initial conditions satsfyng the periodicity condition
% hvalue periodicity condition for the solution
%
%%

C = zeros(2); % conservative system
unit2fminunc = 1000; % pass from m to mm
x0 = X0*unit2fminunc;

optionsOpti = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','TolX',1e-5,'TolFun',1e-5);
[dxOpti, hvalue, ~,~] = fminunc(@(dx) h2minNNM(M,K,C,a,dx,period,unit2fminunc,x0,0),[0;0;0;0],optionsOpti);

dxOpti = dxOpti/unit2fminunc;
signX0 = sign(X0);
signX0(signX0==0) = 1;
IC = X0 + signX0.*dxOpti;

end

