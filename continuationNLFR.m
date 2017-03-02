function [AmplX,f,IC,hvalue] = continuationNLFR(K,M,C,a,f0,step,df_end,x0,Ampl)
%%
% [AmplX,f,IC,hvalue] = continuationNLFR(K,M,C,a,f0,step,df_end,x0,Ampl)
%
% sequential continuation with shooting method for the NLFR computation
% K,M, C matrices of the linear system
% a coefficients of the nonlinearity
% f0, step, df_end respectively, the starting frequency, the step and the
% final increment of the frequency
% x0 inital guess for the first shooting method
% Ampl amplitude of the excitation force
%
% AmplX is the maximum dislacement along one period of the reference degree 
% of freedom , f the frequencies, IC the initial conditions satisfaying the
% periodicity condition, and hvalue the value of the periodicity condition 
% for the solutions.
%
% ex: [AmplX,f,IC,hvalue] = continuationNLFR(K,M,C,a,24,0.01,10,[1e-4 -1e-4 0 0]',25)
%%

% determine the initiale condition for the first shooting
[t,x] = integration(M,K,C,a,x0,10,Ampl,f0,0);
indx1 = find(abs(x(:,3))<=5e-5,1,'last');
x0 = x(indx1,:)';
% determine the phase of the force
phaseF = mod(2*pi*f0*t(indx1),2*pi);

dofMeas = 2;
nf = floor(df_end/step)+1;
f = f0:step:(f0+df_end);
IC = zeros(4,nf);
AmplX = zeros(1,nf);
hvalue = zeros(1,nf);

period = 1/f0;          
[IC(:,1), hvalue(1)] = shootingNLFR(period,x0,M,K,C,a,Ampl,phaseF);
[t,x] = integration(M,K,C,a,IC(:,1),6*period,Ampl,f0,phaseF);
AmplX(1) = max(x(:,dofMeas));

h1 = figure;
plot(t,x)
grid on
h2 = figure;
title('Nonlinear frequency response')
xlabel('frequency [Hz]')
ylabel('displacement [m]')
hold on
plot(f(1),AmplX(1),'bo');
pause(0.1)

for i = 2:nf
    period = 1/f(i);
                           
    [IC(:,i), hvalue(i)] = shootingNLFR(period,IC(:,i-1),M,K,C,a,Ampl,phaseF);
    
    [t,x] = integration(M,K,C,a,IC(:,i),6*period,Ampl,f(i),phaseF);
    AmplX(i) = max(x(:,dofMeas));
    
    figure(h1)
    plot(t,x(:,1:2))
    grid on
    figure(h2);
    plot(f(i),AmplX(i),'bo');
    pause(0.1)
end

figure(h2);
hold off

end
