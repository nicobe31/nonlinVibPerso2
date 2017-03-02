function [F,t] = forces(type,param1,param2, param3, param4, param5)
%%
% [F,t] = forces(type,param1,param2)
% 
% This function determine the force et the corresponding time vector.
% 
% Arguments: type: type of force: 'sineSweep' or 'multisine'
% if type= sineSweep: param1 = amplitude [N] 
%                     param2 = starting frequency [Hz]
%                     param3 = ending frequency [Hz]
%                     param4 = sweep rate [Hz/s]
%                     param5: not used
% if type= multisine: param1 = amplitude [N] 
%                     param2 = starting frequency [Hz]
%                     param3 = ending frequency [Hz]
%                     param4 = number of point for the time vector
%                     param5 = number of frequencies used to represent the
%                     door function
%%

switch type
    case 'sineSweep'
    %% linear sine sweep
    ampl = param1; %N
    f0 = param2; %Hz
    f1 = param3; %Hz
    r = param4; %rate Hz/s
    
        % parameters
        phi0 = 0; %initial phase angle degree
        samplingFreq = 300; %Hz

    phi0 = phi0*pi/180; % °-> rad
    t0 = 0;
    t1 = (f1-f0)/r;
    seriet = t0:1/samplingFreq:t1;

    arg1 = 2*pi*f0*(seriet-t0);
    arg2 = 2*pi*r*(seriet-t0).^2/2;
    arg = arg1 + arg2 + phi0;

    F = ampl*sin(arg);
    
    t = seriet;

    case 'multisine' 
    %% multisine with random phase
    ampl = param1; 
    f0 = param2;
    f1 = param3;
    Nt = param4;
    Neff = param5;
    
        % parameters 
        t0 = 0;
        t1 = 20;
        repetition = 4; % power of 2
       
    serieNt = Nt/repetition;
    mod(serieNt,1)
    seriet = linspace(t0,t1/repetition,serieNt);
    dt = seriet(2)-seriet(1);
    samplingFreq = 1/dt; %Hz

    N = Neff*samplingFreq/(f1-f0);
    k0 = floor(f0*N/samplingFreq);
    k1 = ceil(f1*N/samplingFreq);

    phik = (2*pi*rand(k1-k0+1,1))-pi;

    arg1 = 1i*2*pi*[k0:k1]'*samplingFreq/N*seriet;
    arg2 = 1i*phik;
    arg = bsxfun(@plus, arg1, arg2);

    sines = exp(arg);

    u  = real(N^(-1/2)*ampl*sum(sines,1));
    
    F = zeros(1,Nt);
    t = zeros(1,Nt);
    F(1:serieNt) = u;
    t(1:serieNt) = seriet;
    for rep = 2:repetition
        F((rep-1)*serieNt + 1: rep*serieNt) = u;
        t((rep-1)*serieNt + 1: rep*serieNt) = seriet + t((rep-1)*serieNt) + dt;
    end
% 
y = fft(F); 
fd = abs(y/Nt);
p2 = fd(1:Nt/2+1);
p2(2:end-1) = 2*p2(2:end-1);
f = samplingFreq*(0:(Nt/2))/Nt;
plot(f,p2)

end


