function Plot_PatchColorFields(coord, connect, wallnodeOnlyWallLine, ...
                               coordPatchx, coordPatchy, coordPatchz, ...
                               FieldNodalPatch, ...
                               MinFieldNodalPatch, ...
                               MaxFieldNodalPatch, ...
                               MarkerSizeind, TransparencyInd)
                               
%% Patch graph for force field
figure('Color','white')
hold on
Plot_Mesh(coord, connect, 0.25, false, false, 0.1, 'k')
Vessel = patch(coordPatchx,coordPatchy,coordPatchz,FieldNodalPatch);
plot3(coord(1,wallnodeOnlyWallLine), ...
      coord(2,wallnodeOnlyWallLine), ...
      coord(3,wallnodeOnlyWallLine), 'k.','MarkerSize',MarkerSizeind)
alpha(Vessel, TransparencyInd) % to set transparency
set(Vessel,'FaceColor','interp','EdgeColor','none')
Vessel.CDataMapping = 'scaled';
caxis([min(MinFieldNodalPatch) max(MaxFieldNodalPatch)])
colormap(jet)
cbar = colorbar('southoutside');
set(cbar,'FontSize',12)
set(get(cbar,'title'),'string','Pa');
set(gca,'FontName','cambria math','FontSize',12)
material Dull
axis off
view(3)
rotate3d
axis equal
box off