function fext = exteriorFroce(Ampl,freq,t,dofCoefs,phase,x0)
%%
% fext = exteriorFroce(Ampl,freq,t,dofCoefs,phase,x0)
%
% Compute a sinusoidale cosinusoidale force of amplitude Ampl N, of frequency
% freq Hz, for the time(s) t (s), scaled by dofCoefs (use to place the
% force for example), an initial phase in rad, and the initial condition of
% the time integration in m and m/s x0 (can be replaced by 1)
% 
%% 

si = sign(x0(1));
si(si==0) = 1;
fext = dofCoefs *-si*Ampl*cos((2*pi*freq*t) + phase);

end

