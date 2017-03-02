function h_norm = h2minNNM(M,K,C,a,dx,tend,unit,X0,Ampl)

%scaling
sf1 = [[1,1]*10^(floor(log10(max(abs(X0(1:2)))))),[1,1]*10^(floor(log10(max(abs(X0(3:4))))))]';
dx = dx.*sf1;

% centering
signX0 = sign(X0);
signX0(signX0==0) = 1;
x0 = X0 +signX0.*dx;

x0 = x0/unit; % from mm to m

freq = 1/tend;

[~,x1] = integration(M,K,C,a,x0,tend,Ampl,freq,0);

h = x1(end,:)'- x0;
h_norm = norm(h);
end
