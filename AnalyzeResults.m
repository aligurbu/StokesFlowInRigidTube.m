%% Analyze Results of StokesFlowInRigidTube

%% Check conservation of mass and forces
[inflow, outflow, inletForce, SumofForces, ...
 SumofMoments, SumofNormofMoments] = ...
          BalanceofMassForcesMoments(FieldPts, BasisFn, NormalV, Weights, ...
                                     inletelem, outletelem, wallelem, ...
                                     elemDofNum, ...                            
                                     unodal, Telem, ...
                                     numGaussPoints);
%% 
MassBalance = (inflow-outflow)*RefVelocity*RefLength^2;
MassBalanceError = ((inflow-outflow)/(inflow+outflow))*100;
ForceBalanceError = (norm(SumofForces)/norm(inletForce))*100;
MomentBalanceError = (norm(SumofMoments)/SumofNormofMoments)*100;
fprintf('\nMass balance; \ndifference of inlet and outlet flow: %8.6f micron^3/s \n', ...
                                                               MassBalance)
fprintf('Mass balance error in tube: %8.6f%%\n',MassBalanceError)
fprintf('Sum of forces error in tube: %8.6f%%\n',ForceBalanceError)
fprintf('Sum of moment error in tube: %8.6f%%\n',MomentBalanceError)

%% Plotting the velocity and traction profiles on the vessel boundary
%% Settings 
PlotSettings

%% A line along the tube length
numNodesOnCircumference = 2*16; % for 16 elements around the circimference
if strcmp(name,'ShortMicrocapillary_16El')
    numNodesAlongTubeLength = 2*16+1; % for 16 element along the length
    wallnodeArray = reshape(wallnode, numNodesOnCircumference, ...
                                      numNodesAlongTubeLength);
    wallnodeOnlyWallLine = wallnodeArray(9,[1 3:numNodesAlongTubeLength 2]);
elseif strcmp(name,'RefinedConstrictedVessel_16El')
    numNodesAlongTubeLength = 2*48+1; % for 48 element along the length
    wallnodeArray = reshape(wallnode, numNodesOnCircumference, ...
                                      numNodesAlongTubeLength);
    wallnodeOnlyWallLine = wallnodeArray(9,[1 3:numNodesAlongTubeLength 2]);
elseif strcmp(name,'LongConstrictedVessel_16El')
    numNodesAlongTubeLength = 2*48+1; % for 48 element along the length
    wallnodeArray = reshape(wallnode, numNodesOnCircumference, ...
                                      numNodesAlongTubeLength);
    wallnodeOnlyWallLine = wallnodeArray(9,[1 3:numNodesAlongTubeLength 2]);
end
% % Depict wallnodeOnlyWallLine on vessel
% figure('Color','white')
% hold on
% Patch_Mesh(coord, connect)
% plot3(coord(1,wallnodeOnlyWallLine), ...
%       coord(2,wallnodeOnlyWallLine), ...
%       coord(3,wallnodeOnlyWallLine), 'k.','MarkerSize',MarkerSizeind)

%% Inlet and Outlet velocity profile
Plot_VectorFields(coord, connect, ...
                  [coord(1,inletnode) coord(1,outletnode)], ...
                  [coord(2,inletnode) coord(2,outletnode)], ...
                  [coord(3,inletnode) coord(3,outletnode)], ...
                  [unodal(1,inletnode) unodal(1,outletnode)], ...
                  [unodal(2,inletnode) unodal(2,outletnode)], ...
                  [unodal(3,inletnode) unodal(3,outletnode)], ...
                  Scaleind, LineWidthind, MaxHeadSizeind)
view([0, 90])
fig = gcf;
fig.NumberTitle = 'off';
fig.Name = 'Inlet and Outlet velocity profile';
if verbose_PrintProfiles
    print('InletOutletVelocityProfile','-dpng','-r0')
end

%% Arrangement of connectivity matrix
[connectArranged,connectArrangedInd] = unique(connect);

%% Set-up for patch plotting
[connectIndArrange, coordPatchx, coordPatchy, coordPatchz] = ...
                                          SetUp_CoordPatch(coord, connect);

%% Computing the unit normal vectors at element nodes
[NormalVNodal, NormalVNodalPatch] = ...
          UnitNormalVectorAtElementNodes(coord, connect, connectIndArrange, ...
                                         numElem, numNodesPerElem, ...
                                         numDofPerNode);
% % Depict NormalVNodalPatch
% Plot_VectorFields(coord, connect, ...
%                   coordPatchx, coordPatchy, coordPatchz, ...
%                   NormalVNodalPatch(:,:,1), ...
%                   NormalVNodalPatch(:,:,2), ...
%                   NormalVNodalPatch(:,:,3), ...
%                   Scaleind, LineWidthind, MaxHeadSizeind)

%% Nodal traction, pressure and shear fields
[TractionNodalPatch, ...
 NormTractionNodalPatch, NormTractionNodalArranged, ...
 MinNormTractionNodalPatch, MaxNormTractionNodalPatch, ...
 PressureNodalPatch, ...
 NormPressureNodalPatch, NormPressureNodalArranged, ...
 MinNormPressureNodalPatch, MaxNormPressureNodalPatch, ...
 ShearNodalPatch, ...
 NormShearNodalPatch, NormShearNodalArranged, ...
 MinNormShearNodalPatch, MaxNormShearNodalPatch] = ...
 TractionPressureShearFields_Nodal_Patch(Telem, ...
                                         NormalVNodal, ...
                                         NormalVNodalPatch, ...
                                         connectIndArrange, ...
                                         connectArrangedInd, ...
                                         numElem, ...
                                         numDofPerNode, ...
                                         numDofPerElem, ...
                                         numNodesPerElem);

%% Patch graph for norm of tractions
Plot_PatchColorFields(coord, connect, wallnodeOnlyWallLine, ...
                      coordPatchx, coordPatchy, coordPatchz, ...
                      NormTractionNodalPatch, ...
                      MinNormTractionNodalPatch, ...
                      MaxNormTractionNodalPatch, ...
                      MarkerSizeind, TransparencyInd)
fig = gcf;
fig.NumberTitle = 'off';
fig.Name = 'NormTractionProfile';
if verbose_PrintProfiles
    print('NormTractionProfile','-dpng','-r0')
end

%% Line graph for norm of tractions
Plot_LineFieldsOnVessel(coord, NormTractionNodalArranged, ...
                               wallnodeOnlyWallLine, ...
                               MarkerSizeind, LineWidthind, ...
                               'Norm of traction, Pa', ...
                               'A line along the vessel wall ($\mu$m)')
fig = gcf;
fig.NumberTitle = 'off';
fig.Name = 'NormTractionLine';
if verbose_PrintProfiles
    print('NormTractionLine','-dpng','-r0')
end

