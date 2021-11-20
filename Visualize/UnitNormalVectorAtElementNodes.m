function [NormalVNodal, NormalVNodalPatch] = ...
          UnitNormalVectorAtElementNodes(coord, connect, connectIndArrange, ...
                                         numElem, numNodesPerElem, ...
                                         numDofPerNode)
%% Computing the unit normal vectors at element nodes
NormalVNodal = zeros(numNodesPerElem,numElem,numDofPerNode);
for mm = 1:numElem
    xi_e = coord(:,connect(:,mm));
    for kk = 1:numNodesPerElem
        [zeta1,zeta2] = getZeta(kk);
        [nhat, ~] = jacobian_9nodequad(zeta1, zeta2, xi_e);
        NormalVNodal(kk,mm,:) = nhat;
    end
end
NormalVNodalPatch = NormalVNodal(connectIndArrange,:,:);
