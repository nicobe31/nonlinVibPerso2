function [E,f,IC,hvalue] = continuationNNM(K,M,a,f0,step,df_end,x0)
%%
% [E,f,IC,hvalue] = continuationNNM(K,M,a,f0,step,df_end,x0)
%
% sequential continuation with shooting method for the NNM computations
% K,M matrices of the linear system
% a coefficients of the nonlinearity
% f0, step, df_end respectively, the starting frequency, the step and the
% final increment of the frequency
% x0 inital guess for the first shooting method
%
% E is the energy in the system, f the frequencies, IC the initial
% conditions satisfaying the periodicity condition, and hvalue the value of
% the periodicity condition for the solutions.
%
% ex: [E,f,IC,hvalue] = continuationNNM(K,M,a,27.5664,0.1,4,[0.015,-0.015,0,0]')
%%

nf = floor(df_end/step)+1;
f = f0:step:(f0+df_end);
IC = zeros(4,nf);
E = zeros(1,nf);
hvalue = zeros(1,nf);

period = 1/f0;          
[IC(:,1), hvalue(1)] = shootingNNM(period,x0,M,K,a);
E(1) = 0.5*IC(1:2,1)'*K*IC(1:2,1) + 0.5*IC(3:4,1)'*M*IC(3:4,1);
h1 = figure;
title('Frequency - Energy plot')
xlabel('log10(E) [J]')
ylabel('frequency [Hz]')
hold on
plot(log10(E(1)),f(1),'bo');
pause(0.1)
for i = 2:nf
    period = 1/f(i);
    
    [IC(:,i), hvalue(i)] = shootingNNM(period,IC(:,i-1),M,K,a);
    E(i) =  0.5*IC(1:2,i)'*K*IC(1:2,i) + 0.5*IC(3:4,i)'*M*IC(3:4,i);
    figure(h1);
    plot(log10(E(i)),f(i),'bo');
    pause(0.1)
end
figure(h1);
hold off

end
