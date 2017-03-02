function h_norm = h2minNLFR(M,K,C,a,dx,tend,unit,X0,Ampl,phaseF) 

%scaling
sf1 = [[1,1]*10^(floor(log10(max(abs(X0(1:2)))))),[1,1]*10^(floor(log10(max(abs(X0(3:4))))))]';
dx = [dx(1:2); 0; dx(3)];
dx = dx.*sf1;

% centering
signX0 = sign(X0);
signX0(signX0==0) = 1;
x0 = X0 +signX0.*dx;

x0 = x0/unit; %repass from mm to m


freq = 1/tend;

[~,x1] = integration(M,K,C,a,x0,tend,Ampl,freq,phaseF);
[~,x2] = integration(M,K,C,a,x1(end,:),tend,Ampl,freq,phaseF);
[~,x3] = integration(M,K,C,a,x2(end,:),tend,Ampl,freq,phaseF);
[~,x4] = integration(M,K,C,a,x3(end,:),tend,Ampl,freq,phaseF);
[~,x5] = integration(M,K,C,a,x4(end,:),tend,Ampl,freq,phaseF);
[~,x6] = integration(M,K,C,a,x5(end,:),tend,Ampl,freq,phaseF);

h = (x1(end,:)'+x2(end,:)'+x3(end,:)'+x4(end,:)'+x5(end,:)'+x6(end,:)' - 6*x0)/6;
% h1 = norm(x1(end,:)'-x0);
% h2 = norm(x2(end,:)'-x0);
% h3 = norm(x3(end,:)'-x0);
% h4 = norm(x4(end,:)'-x0);
% h5 = norm(x5(end,:)'-x0);
% h6 = norm(x6(end,:)'-x0);
% h_max = max([h1,h2,h3,h4,h5,h6]);
h_norm = norm(h);
end
