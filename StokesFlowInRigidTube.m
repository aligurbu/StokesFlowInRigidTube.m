%% Direct boundary element method to analyze 
%% the extracellular plasma flow in a rigid vessel
%% 
addpath 'Integration/'
addpath 'Interpolation/'
addpath 'Visualize/'
addpath 'Models/'

verbose_Plot = false;
verbose_Patch = false;
Starttime = tic;

%% Input the model and parameters for the analysis from Models folder
Load_ShortMicrocapillary_16El


%%
ModelTime = toc(Starttime)