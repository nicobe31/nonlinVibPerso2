function [fitresult, gof] = createFit2(x, z, eqFit)
%%
% [fitresult, gof] = createFit2(x, z, eqFit)
%
% create the 2d fit of the variable z in function of x with the mathemtical
% model eqFit
%%

% Set up fittype and options.
ft = fittype(eqFit);

% Fit model to data.
[fitresult, gof] = fit( x', z', ft );

% Plot
figure( 'Name', 'RSM' );
h = plot( fitresult, x, z );
xlabel x
ylabel xdd
grid on


