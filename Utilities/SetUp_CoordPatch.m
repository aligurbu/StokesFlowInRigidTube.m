function [connectIndArrange, coordPatchx, coordPatchy, coordPatchz] = ...
          SetUp_CoordPatch(coord, connect)
%% Set-up for patch plotting
% connectIndArrange triangularizes the 9 node quadrilateral element for 
% patch ploting so that the 9th node would be used in the interpolation of 
% field on the element. Basically, quadrilateral is subdivided to two 
% triangles by using connectIndArrange.
% Coordinates of element nodes 
connectIndArrange = [1 5 2 6 3 9 1 3 7 4 8 1 9 3];
connect_rearrange = connect(connectIndArrange,:);
coordPatchx = reshape(coord(1,connect_rearrange),size(connect_rearrange));
coordPatchy = reshape(coord(2,connect_rearrange),size(connect_rearrange));
coordPatchz = reshape(coord(3,connect_rearrange),size(connect_rearrange));