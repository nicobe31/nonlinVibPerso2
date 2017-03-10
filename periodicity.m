function h = periodicity(x0,period,M,K,C,a,AmplF,freqF,phaseF)
%%
% h = periodicity(x0,period,M,K,C,a,AmplF,freqF,phaseF)
%
% AmpF, freqF, phaseF are the parameters of the force, if omited set to 0
%
% If do NNM C must be zeros(n)
%% 

if nargin <= 6
    AmplF = 0;
    freqF = 0;
    phaseF = 0;
end

[~,x] = integration(M,K,C,a,x0,period,AmplF,freqF,phaseF);

h = x(end,:)' - x0;

end
