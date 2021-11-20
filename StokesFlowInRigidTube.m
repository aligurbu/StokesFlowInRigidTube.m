%% Direct boundary element method to analyze 
%% the extracellular plasma flow in a rigid vessel
%% 
addpath 'Integration/'
addpath 'Interpolation/'
addpath 'Visualize/'
addpath 'Models/'
addpath 'Tests/'

verbose_Plot = false;
verbose_Patch = false;
Starttime = tic;

%% Input the model and parameters for the analysis from Models folder
Load_ShortMicrocapillary_16El
% Load_RefinedConstrictedVessel_16El
% Load_LongConstrictedVessel_16El

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
verbose_PrintProfiles = false;

%% Analyze the results 
AnalyzeResults