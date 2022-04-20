%% Direct boundary element method to analyze
%% the extracellular plasma flow in a rigid vessel
%%
addpath 'Integration/'
addpath 'Interpolation/'
addpath 'Visualize/'
addpath 'Models/'
addpath 'Tests/'
addpath 'Utilities/'

verbose_Plot = false;
verbose_Patch = false;
verbose_PrintProfiles = false;
Starttime = tic;

%% Input the model and parameters for the analysis from Models folder
%% Vessel models: Load mesh
load ShortMicrocapillary_16El.mat
% load RefinedConstrictedVessel_16El.mat
% load LongConstrictedVessel_16El.mat
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

%% Field points (Gauss quadrature nodes)
[FieldPts, BasisFn, NormalV, Weights] = ...
                   FieldProperties(coord, connect, numElem, ...
                                   numDofPerNode, numNodesPerElem, gx, gw);

%% Solving
Solution = SolveBEM_GMRES(coord, connect, inletelem, outletelem, ...
                          elemDofNum, ...
                          NeumannDofs, NeumannNode, DirichletElem, ...
                          FieldPts, NormalV, Weights, BasisFn, ...
                          Telem, ...
                          grx, grw, gtx, gtw, mu, numGaussPoints, ...
                          numNodes, numDofPerNode, ToleranceGMRES);

%% Collect the solutions
%% Collect velocities
unodal(NeumannDofs) = Solution(NeumannDofs);

%% Collect tractions
Telem(:,DirichletElem) = Solution(elemDofNum(:,DirichletElem));

%%
ModelTime = toc(Starttime)

%% Analyze and post-process the results
%% Analyze the results
AnalyzeResults

%% Post-processing: the flow in the vessel
PostProcessingResults