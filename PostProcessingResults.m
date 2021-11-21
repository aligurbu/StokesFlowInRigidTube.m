%% Post-Processing
StartPostProcessingTime = tic;

%% Set-up the evaluation points
SetUp_EvaluationPoints

%% Compute the velocity field inside the vessel
VelocityInsideDomain = PostProcessing(EvaluationPts, coord, connect, ...
                                      inletelem, elemDofNum, ...
                                      DirichletElem, NeumannElem, ...
                                      FieldPts, NormalV, Weights, BasisFn, ...
                                      unodal, Telem, ...
                                      grx, grw, gtx, gtw, mu, numGaussPoints, ...
                                      numDofPerNode, numDofPerElem);
PostProcessingTime = toc(StartPostProcessingTime)

%% Dimensionalization of the parameters
coord = coord*RefLength; % micron
EvaluationPts = EvaluationPts*RefLength; % micron
LinePts = LinePts*RefLength; % micron
CirclePtsmid = CirclePtsmid*RefLength; % micron
CirclePtsinlet = CirclePtsinlet*RefLength; % micron
CirclePtsoutlet = CirclePtsoutlet*RefLength; % micron
VelocityInsideDomain = VelocityInsideDomain*RefVelocity*10^(-3); %mm/sec

%% Plot Settings 
PlotSettings

%% Contour graph of speed of the fluid inside the vessel
%% Compute the speed of the fluid flow
[LineSpeed, CrossSectionSpeedmid, CrossSectionSpeedinlet, ...
            CrossSectionSpeedoutlet, MinCaxis, MaxCaxis] = ...
            SpeedComputations(LinePts, LinePtsX, ...
                              CirclePtsmid, CirclePtsXmid, ...
                              CirclePtsinlet, CirclePtsXinlet, ...
                              CirclePtsoutlet, CirclePtsXoutlet, ...
                              VelocityInsideDomain);

%% Contour graph 
ContourSpeedInsideVessel

%% 
% %% Maximum velocity
% InletPressure = InletPressure*RefPressure; % Pa = N/m^2
% TubeRadius = TubeRadius*RefLength; % micron
% TubeLength = TubeLength*RefLength; % micron
% mu = mu*RefViscosity; % Pa.s = N.s/m^2
% Vmax = ((InletPressure * TubeRadius^2) / (4*mu*TubeLength))*10^(-3) %mm/s