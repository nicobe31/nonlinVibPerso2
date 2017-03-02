function fNL = nonlinearForce(x,xd,a)
%%
% fNL = nonlinearForce(x,xd,a)
%
% This function use the defined mathematical model of the nonlinearity, the
% matrices of the displacement and velocities and the coefficients of the
% nonlinearity to compute the nonlinear force applied on all dof
%%

xx = x(2,:) - x(1,:);
xxd = xd(2,:) - xd(1,:);

f1 = a(1)*xx.^2 + a(2)*xx.^3;

fNL = [-1;1]*f1;

end
