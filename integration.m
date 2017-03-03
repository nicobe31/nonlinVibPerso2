function [t,x] = integration(M,K,C,a,x0,tend,AmplF,freqF,phaseF,absTol,relTol)

if nargin <=6
    AmplF = 0;
    freqF = 0;
    phaseF = 0;
end
if nargin <= 9
    absTol = 1e-8;
    relTol = 1e-8;
end

optionsOdeo = odeset('RelTol',relTol, 'AbsTol',absTol);

[t,x] = ode45(@(t,x) eqMot(t,x,M,K,C,a,AmplF,freqF,x0,phaseF),[0, tend],x0,optionsOdeo);

end

function xd = eqMot(t,x,M,K,C,a,Ampl,freq,x0,phase)

dof = 2;

L = [zeros(dof), eye(dof);
     -M\K, -M\C];

fnl = nonlinearForce(x(1:2),x(3:4),a);

fext = exteriorFroce(Ampl,freq,t,[1;0],phase,x0);

gnl = [zeros(dof,1);
       M\fnl];
gext = [zeros(dof,1);
        M\fext];
    
xd = L*x - gnl + gext;

end
