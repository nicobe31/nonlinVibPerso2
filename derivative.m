function [hz, ht] = derivative(x0, stepX, period, M, C, K, a, AmplF,freqF,phaseF)

if nargin <= 7
    AmplF = 0;
    freqF = 0;
    phaseF = 0;
end

nx = length(x0);
hz = zeros(nx);

h0 = periodicity(x0,period,M,K,C,a,AmplF,freqF,phaseF);

for i=1:nx
    dx = zeros(nx,1);
    dx(i) = stepX;
    IC = x0 + dx;
    
    h = periodicity(IC,period,M,K,C,a,AmplF,freqF,phaseF);
    hz(:,i) = (h - h0)/stepX;
end

dofFext = [1;0];
fext = exteriorFroce(AmplF,freqF,period,dofFext,phaseF,x0);
fnl = nonlinearForce(x0(1:2),x0(3:4),a);

L = [zeros(nx/2), eye(nx/2);
    -M\K, -M\C];
gnl = [zeros(nx/2,1); M\fnl];
gext = [zeros(nx/2,1); M\fext];

ht = L*x0 -gnl + gext;

end
