addpath '../../Interpolation/'
addpath '../../Visualize/'
addpath '../../Models/'

verbose_Plot = false;
verbose_Patch = true;
verbose_Patch_XY = true;

%% Load mesh
%% Vessel models
% load ShortMicrocapillary_16El.mat
% load RefinedConstrictedVessel_16El.mat
load LongConstrictedVessel_16El.mat
name = nameVessel;

%%
CellVolume = 93.78; % micron^3
RefLength = ((3*CellVolume)/(4*pi))^(1/3); % micron

if verbose_Plot
    figure('Color','white')
    Plot_Mesh(coord, connect, 0.25, true, true)
    axis off
end
if verbose_Patch
    figure('Color','white')
    Patch_Mesh(coord, connect)
end

%%
set(gcf,'PaperPositionMode','auto')
if verbose_Plot
    view([0 90])
    print([name,'_Plot'],'-dsvg') 
end
if verbose_Patch
    print([name,'_Patch'],'-dsvg') 
end
if verbose_Patch_XY
    view([0 0])
    print([name,'_Patch_XY'],'-dsvg') 
end
