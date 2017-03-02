function [t,x] = integration(M,K,C,a,x0,tend,Ampl,freq,phase)

optionsOdeo = odeset('RelTol',1e-8, 'AbsTol',1e-8);

[t,x] = ode45(@(t,x) eqMot(t,x,M,K,C,a,Ampl,freq,x0,phase),[0, tend],x0,optionsOdeo);

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
