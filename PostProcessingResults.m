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