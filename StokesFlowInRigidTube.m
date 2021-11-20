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


%%
ModelTime = toc(Starttime)