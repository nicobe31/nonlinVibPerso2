function [a,MSE] = RFS(res,loc,M,K,C)
%%
% [a,MSE] = RFS(res,loc,M,K,C)
% 
% This function compute the coeffficients of the nonlinearity and the mean
% square error.
% res is a structure of the results, loc is index of the measured force,
% M,K,C are the matrices of the linear system
% 
% ex: [a,MSE] = RFS(resR_40_1_22_40,2,M,K,C)
%%

xx = res.x(2,:) - res.x(1,:);
xxd = res.xd(2,:) - res.xd(1,:);
% xxdd = res.xdd(loc(1),:);
m = M(loc(1),loc(1));

linAcc = K*res.x + C*res.xd;
nonLinAcc = res.xdd + linAcc;
Z = -m*nonLinAcc(loc(1),:);

eqFit = 'a2*x.^2 + a3*x.^3';
[fitresult, ~] = createFit2(xx, Z, eqFit);
a =  coeffvalues(fitresult);

N = length(xx);
sigma = var(Z);
xddhat = feval(fitresult,xx);
MSE = 100/(N*sigma)*sum((Z- xddhat').^2);
end
