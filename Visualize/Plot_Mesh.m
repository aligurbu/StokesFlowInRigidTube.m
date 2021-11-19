function Plot_Mesh(coord, connect, Fineness, PlotNodeCoordinates, ...
                                             PlotUnitNormalVector, ...
                                             LineWidth, LineColor)
%%
% Plot the mesh from given coordinates of nodes and connectivity matrix 
% of elements. 
% Fineness of mesh \in (0,1]. Smaller value is corresponding to better 
% resolution. 
% if PlotNodeCoordinates is true then plot the node coordinates.
% if PlotUnitNormalVector is true then plot the unit normal vector of
% element at the 9th node.

%% Settings 
if nargin >5
    LineWidthind = LineWidth; % Line width
else
    LineWidthind = 1.1; % Line width
end
if nargin > 6
    LineColorind = LineColor; % Line color
else
    LineColorind = 'b'; % Line color
end
MarkerSizeind = 10; % Marker size
Scaleind = 0.25; % Scaling index for quiver3
MaxHeadSizeind = 0.7; 

%% Drawing the elements 
numEl = size(connect,2);
s = -1:Fineness:1;
s1 = [s ones(length(s),1)' flip(s) -ones(length(s),1)'];
s2 = [-ones(length(s),1)' s ones(length(s),1)' flip(s)];
hold on
for el = 1:numEl
    xie = coord(:,connect(:,el));
    X = zeros(size(xie,1),length(s1));
    for k = 1:length(s1)
        X(:,k) = xie*interpolate_9nodequad(s1(k), s2(k));
    end
    plot3(X(1,:),X(2,:),X(3,:),[LineColorind,'-'],'LineWidth',LineWidthind)
end
axis equal
xlabel('X'), ylabel('Y'), zlabel('Z')
view(3)

if PlotNodeCoordinates
    plot3(coord(1,:),coord(2,:),coord(3,:),'r.','MarkerSize',MarkerSizeind)
    axis equal
    xlabel('X'), ylabel('Y'), zlabel('Z')
end

if PlotUnitNormalVector
    X = zeros(1,numEl);
    Y = zeros(1,numEl);
    Z = zeros(1,numEl);
    nhatX = zeros(1,numEl);
    nhatY = zeros(1,numEl);
    nhatZ = zeros(1,numEl);
    for n = 1:numEl
        [nhat, ~] = jacobian_9nodequad(0, 0, coord(:,connect(:,n)));
        nhatX(n) = nhat(1);
        nhatY(n) = nhat(2);
        nhatZ(n) = nhat(3);
        X(n) = coord(1,connect(9,n));
        Y(n) = coord(2,connect(9,n));
        Z(n) = coord(3,connect(9,n));
    end
    quiver3(X,Y,Z,nhatX,nhatY,nhatZ,Scaleind,'k', ...
            'LineWidth',LineWidthind,'MaxHeadSize',MaxHeadSizeind),
    axis equal,
    xlabel('X'), ylabel('Y'), zlabel('Z')
end
axis equal
rotate3d