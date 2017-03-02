% Second version of the continuation method,
% with a pseudo- arc-length method

function [AmplX,f,IC,hvalue] = continuationNLFR_02(K,M,C,a,f0,step,df_end,x0,Ampl)

[t,x] = integration(M,K,C,a,x0,10,Ampl,f0,0);
plot(t,x);
indx1 = find(abs(x(:,3))<=5e-5,1,'last');
x0 = x(indx1,:)';
x0(3) = 0;
phaseF = mod(2*pi*f0*t(indx1),2*pi);

dofMeas = 2;
nf = floor(df_end/step)+1;
% f = f0:step:(f0+df_end);
IC = zeros(4,nf);
AmplX = zeros(1,nf);
hvalue = zeros(1,nf);

f(1) = f0;
period = 1/f(1);
[IC(:,1), hvalue(1),iter, flag] = shootingNLFR(period,x0,M,K,C,a,Ampl,phaseF);
[t,x] = integration(M,K,C,a,IC(:,1),6*period,Ampl,f0,phaseF);
AmplX(1) = max(x(:,dofMeas));

iterOpti = 5;

h1 = figure;
plot(t,x)
grid on
h2 = figure;
hold on
plot(f(1),AmplX(1),'bo');
pause(0.1)

f(2) = f0 + step;
period = 1/f(2);          
[IC(:,2), hvalue(2),iter,flag] = shootingNLFR(period,IC(:,1),M,K,C,a,Ampl,phaseF);
[t,x] = integration(M,K,C,a,IC(:,1),6*period,Ampl,f(2),phaseF);
AmplX(2) = max(x(:,dofMeas));

figure(h1);
plot(t,x)
grid on
figure(h2)
plot(f(2),AmplX(2),'bo');
pause(0.1)

stepf_prec = step;
stepf = step;
dx = IC(:,2) - IC(:,1);
stepX = norm(dx);
stepff(1:2) = [step, step];
stepXX(1:2) = [stepX,stepX];

for i = 3:nf
    
    stepf_prec = stepf;
    stepf = min(stepf*(iterOpti+2)/(iter+2), 0.1);
    stepf = max(stepf,0.01);
    stepX = stepX*stepf./stepf_prec;
    stepff(i) = stepf;
    stepXX(i) = stepX;
    df = f(i-1) - f(i-2);
    dx = IC(:,i-1) - IC(:,i-2);
    
    f(i) = f(i-1) + stepf*df/norm(df);
    x0 = IC(:,i-1) + stepX*dx/norm(dx);
    
    period = 1/f(i);
                           
    [IC(:,i), hvalue(i),iter,flag] = shootingNLFR(period,x0,M,K,C,a,Ampl,phaseF);
    
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
