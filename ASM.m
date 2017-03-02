 function ASM(res,index,M,K,C,displayDensity,a)
%%
% ASM(res,index,M,K,C,displayDensity,a)
%
% Display the Acceleration Surface Method for the data structure res given in
% argument. index is a vector containing the indices between which the
% acceleration surface is computed. Set one of the index to 0 to compute
% the AS between the dof and its fixation. M,L,C are the matrices of the
% linear system.
% displayDensity can either be 'normal','light' or 'extralight': this
% parameter change the density of point represented in the 3d view. If
% omited is set to 'light'
% 'a' is the coefficients of the nonlinearity (as nedeed in the function
% nonlinearForces). If omitted is setted to [0,0];

%%
if nargin < 6
    displayDensity = 'light';
    a = zeros(1,2);
end
if nargin < 7
    a = zeros(1,2);
end
ratioWidth = 0.02; %ratio of the width of the slice
colorSlice = ['bx'; 'gx'; 'rx'; 'cx'; 'mx'; 'yx'];

switch displayDensity
    case 'normal'
        nskip = 1;
    case 'light'
        nskip = 10;
    case 'extralight'
        nskip = 50;
    otherwise 
        nskip = 10;
end

%% data preparation
lengthQ = length(res.t);
deltaQ = zeros(3,lengthQ);
testIndex = find(index==0,1);

if isempty(testIndex);
    deltaQ(1,:) = res.x(index(1),:) - res.x(index(2),:);
    deltaQ(2,:) = res.xd(index(1),:) - res.xd(index(2),:);
else
    index(testIndex) = [];
    deltaQ(1,:) = res.x(index(1),:);
    deltaQ(2,:) = res.xd(index(1),:);
end    

linAcc = K*res.x + C*res.xd;
fNL = nonlinearForce(res.x,res.xd,a);
nonLinAcc = res.xdd + linAcc + fNL;

deltaQ(3,:) = -nonLinAcc(index(1),:);

%% ASM 3D plot
fig3d = figure('Name','ASM 3D');
plot3(deltaQ(1,1:nskip:end),deltaQ(2,1:nskip:end),deltaQ(3,1:nskip:end),'ko');
title('ASM','FontSize', 16);
xlabel('x: displacement','FontSize', 14);
ylabel('y: velocity','FontSize', 14);
zlabel('z: acceleration','FontSize', 14);

%% Slice
numberSlice = 1;
promp = sprintf('What slice do you want? axes-angleD-offset (ex: x-045-0.5) or no \n');
sliceStrg = input(promp,'s');

while ~strcmpi(sliceStrg,'no')
    hyphenInd = find(sliceStrg == '-');
        
    % axe
    sliceAxeStrg = sliceStrg(1);
    switch sliceAxeStrg
        case 'x'
            sliceAxe = 1;
        case 'y'
            sliceAxe = 2;
        case 'z'
            sliceAxe = 3;
        otherwise
            sliceAxe = 2;
    end

    % value
    sliceValue = str2double(sliceStrg(hyphenInd+1:end));
    widthSliceHalf = ratioWidth*(max(deltaQ(sliceAxe,:)) - min(deltaQ(sliceAxe,:)))/2;

    % find indice in deltaQ that are in the slice
    iSlice = find((deltaQ(sliceAxe,:)>(sliceValue-widthSliceHalf))&(deltaQ(sliceAxe,:)<(sliceValue+widthSliceHalf)));
        
    % plot
    figure(fig3d);
    hold on
    plot3(deltaQ(1,iSlice),deltaQ(2,iSlice),deltaQ(3,iSlice),colorSlice(mod(numberSlice,6),:));
     
    axePlot2d = 1:3;
    axePlot2d(sliceAxe) = [];
    axeName = {'Displacement';'Velocity';'Acceleration'};
    fig2d = figure('Name','ASM slice');%,'FontSize', 16);
    plot(deltaQ(axePlot2d(1),iSlice),deltaQ(axePlot2d(2),iSlice),colorSlice(mod(numberSlice,6),:),'LineWidth',2);
    titleStrg = sprintf('Slice number %i', numberSlice);
    title(titleStrg);
    xlabel(axeName{axePlot2d(1)},'FontSize', 14);
    ylabel(axeName{axePlot2d(2)},'FontSize', 14);
    
    numberSlice = numberSlice + 1;
    sliceStrg = input(promp,'s');
end
end
        
