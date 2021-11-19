%% Fluid properties
%% Viscosity of exterior fluid
mu_out = 1.2*10^(-3); % Pa.s 
mu = mu_out/RefViscosity; 

%% Applied inlet pressure 
PressureGradient = 4; % mmHg/mm
% 1 mmHg = 133.322 Pa
PressureGradient = PressureGradient*0.133322; % Pa/micron
PressureDrop = PressureGradient*(TubeLength*RefLength); % Pa
InletPressure = PressureDrop/RefPressure; % Inlet pressure

%% Gauss quadrature for the regular integrals
%% $\chi \not\in \Gamma_e$ 
numGaussPoints = 10;
[gx, gw] = lgwt(numGaussPoints, -1, 1);

%% Gauss quadrature for the weakly- and nearly-singular integrals 
%% $\chi \in \Gamma_e$ or dmin/LengE < 1
numGaussPointsPolar = 20;
%% Gauss point and weights along the radius and theta in polar coordinates
[grx, grw] = lgwt(numGaussPointsPolar, -1, 1); 
[gtx, gtw] = lgwt(numGaussPointsPolar, -1, 1); 

%% The tolerance of the GMRES method
ToleranceGMRES = 10^(-6);

%% Set-up
numNodes = size(coord,2); % Total number of nodes in the model

numElem = size(connect,2); % Total number of elements in the model

numDofPerNode = size(coord,1); % Number of velocity components

numNodesPerElem = size(connect,1); % Number of nodes per element

numDofPerElem = numNodesPerElem * numDofPerNode; 
                                 % Number of DOF associated with an element

%% Set the nodal boundary conditions
% bcs == 1; Dirichlet B.C; velocity is known
% bcs == 0; Neumann B.C; traction is known
bcs = zeros(size(coord));
bcs(:) = 1;
bcs(:,setdiff(inletnode,wallnode)) = 0;
bcs(:,setdiff(outletnode,wallnode)) = 0;

DirichletNode = wallnode;
NeumannNode = union(setdiff(inletnode,wallnode), ...
                    setdiff(outletnode,wallnode));

DirichletDofs = find(bcs == 1); % Dofs of the wall nodes
NeumannDofs = find(bcs == 0);   % Dofs of the inlet/outlet nodes 
                                % excluding dofs related to the edge nodes.

%% Set element boundary conditions
bcselem = zeros(1, size(connect,2));
bcselem(1,wallelem) = 1; % Dirichlet B.C; velocity is known

DirichletElem = wallelem; % = find(bcselem);
NeumannElem = union(inletelem,outletelem); % = find(~bcselem);

%% Set prescribed (inlet) traction
Telem = zeros(numDofPerElem, numElem);
%% Prescribe the inlet pressure on the inlet element nodes in x-direction
Telem(1:numDofPerNode:numDofPerElem,inletelem) = InletPressure; 

%% Set prescribed velocity (all zero for now)
unodal = zeros(numDofPerNode, numNodes);

%% Array to store element DOF numbers
elemDofNum = zeros(numDofPerElem, numElem);
for m = 1:numElem
    for n = 1:numNodesPerElem
        elemDofNum((n-1)*numDofPerNode+(1:numDofPerNode),m) = ...
                          (connect(n,m)-1)*numDofPerNode+(1:numDofPerNode);
    end
end