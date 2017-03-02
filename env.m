function [xUpper, xLower, freqUpper, freqLower] = env(resMatrix,grpDisplay)
%% 
% function [xUpper, xLower, freqUpper, freqLower] = env(resMatrix,grpDisplay)
%
% Arguments:
% resMatrix: matrix containing the res structure from simulations
% grpDisplay: 'dof', 'simu': decide how to group the display, by dof's or by simulation. If omitted default is by dof.
% 
% Outputs: 
% xUpper, freqUpper: upper envelope of the displacement and corresponding time. Given in an array form of dimension nbreSimulaiton x nDof.
% xLower, freqLower: lower envelope of the displacement and corresponding time. Given in an array form of dimension nbreSimulaiton x nDof.


if nargin < 2
     grpDisplay = 'dof';
end
 
nSimu = length(resMatrix);
[ndof, nSample] = size(resMatrix{1}.x);

xUpper = cell(nSimu,ndof);
xLower = cell(nSimu,ndof);
freqUpper = cell(nSimu,ndof);
freqLower = cell(nSimu,ndof);
color = ['b'; 'g'; 'r'; 'c'; 'm'; 'y'];

%% Envelope computation

for i = 1:nSimu
	for j = 1:ndof
	[xUpper{i,j},indUpper] = findpeaks(resMatrix{i}.x(j,:));
	[xLower{i,j},indLower] = findpeaks(-resMatrix{i}.x(j,:));
    xLower{i,j} = -xLower{i,j};
    freqUpper{i,j} = resMatrix{i}.f(indUpper);
    freqLower{i,j} = resMatrix{i}.f(indLower);
	end
end

%% Display

if strcmp(grpDisplay, 'dof')

    hp = zeros(nSimu,1);
    leg = cell(nSimu,1);
    for j = 1:ndof
        name = sprintf('dof %i',j);
        figure('name',name);
        hold on
        for i = 1:nSimu
            hp(i) = plot(freqUpper{i,j},xUpper{i,j},color(mod(i-1,6)+1,:));
            plot(freqLower{i,j},xLower{i,j},color(mod(i-1,6)+1,:));
            leg{i} = sprintf('simu %i',i);
        end
        title(sprintf('envelope of the displacement'),'FontSize', 16);
        legend(hp,leg);
        xlabel('Frequency [Hz]','FontSize', 14);
        ylabel('Displacement [m]','FontSize', 14);
        
    end

elseif strcmp(grpDisplay,'simu')

    hp = zeros(ndof,1);
    leg = cell(ndof,1);
    for i = 1:nSimu
        name = sprintf('simu %i',i);
        figure('name',name);
        hold on
        for j = 1:ndof
            hp(j) = plot(freqUpper{i,j},xUpper{i,j},color(mod(i-1,6)+1,:));
            plot(freqLower{i,j},xLower{i,j},color(mod(i-1,6)+1,:));
            leg{j} = sprintf('dof %i',j);
        end
        title(sprintf('envelope of the displacement, simu %i',i));
        legend(hp,leg);
        xlabel('frequency [Hz]');
        
    end       

end
