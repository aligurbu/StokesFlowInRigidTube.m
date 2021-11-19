%% Load mesh
%% Vessel models
load ShortMicrocapillary_16El.mat
name = nameVessel;

%% Plotting the geometry of the model
if verbose_Plot
    figure('Color','white')
    Plot_Mesh(coord, connect, 0.25, true, false)
    axis off
end
if verbose_Patch
    figure('Color','white')
    Patch_Mesh(coord, connect)
end

%% Reference parameters
ReferenceParameters

%% Parameters 
ParametersForTheAnalysis