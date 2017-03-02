% Third version of the continuation function. 
% Try of a kind of arc-length method

function [E,f,IC,hvalue] = continuationNNM_03(K,M,a,f0,step,df_end,x0)

iterOpti = 1;
maxIter = 200;
C = zeros(2);
Ampl = 0; 
freqF = 0;

% 1
period = 1/f0;  
f(1) = f0;
[IC(1:4,1), hvalue(1)] = shootingNNM_03(period,x0,M,K,a);
IC(5,1) = period;
E(1) = 0.5*IC(1:2,1)'*K*IC(1:2,1) + 0.5*IC(3:4,1)'*M*IC(3:4,1);
h1 = figure;
hold on
plot(log(E(1)),f(1),'bo');
pause(0.1)

% 2
f(2) = f0+step;
period = 1/f(2);          
[IC(1:4,2), hvalue(2)] = shootingNNM_03(period,IC(1:4,1),M,K,a);
IC(5,2) = period;
E(2) = 0.5*IC(1:2,2)'*K*IC(1:2,2) + 0.5*IC(3:4,2)'*M*IC(3:4,2);
plot(log(E(2)),f(2),'bo');
pause(0.1);

i =3;
while f(i-1) <= (f0 + df_end) && i < maxIter
    
    difft = IC(:,i-1) - IC(:,i-2);
    t = difft/norm(difft);
    x0 = IC(:,i-1) + step*t;
    
    [IC(:,i), hvalue(i),iter] = shooting11(x0,M,K,a);
    E(i) =  0.5*IC(1:2,i)'*K*IC(1:2,i) + 0.5*IC(3:4,i)'*M*IC(3:4,i);
    f(i) = 1/IC(5,i);
    figure(h1);
    plot(log(E(i)),f(i),'bo');
    pause(0.1)
    if iter == 0;
        iter = 1;
    end
    
    step = step*iterOpti/iter;
    
    i = i+1;
end

figure(h1);
hold off

end
