%% Contour graph of speed of the fluid inside the vessel
figure('Color','white')
hold on
Patch_Mesh(coord, connect, 0.1)
% LineSpeed
ContourLineSpeed(LinePts(:,:,1), LinePts(:,:,2), LineSpeed);
% Middle section speed
ContourCrossSectionSpeed(CirclePtsmid(:,:,2), CirclePtsmid(:,:,3), ...
                         CrossSectionSpeedmid, CirclePtsmid(1,1,1));
% inlet section speed
ContourCrossSectionSpeed(CirclePtsinlet(:,:,2), CirclePtsinlet(:,:,3), ...
                         CrossSectionSpeedinlet, CirclePtsinlet(1,1,1));
% outlet section speed
ContourCrossSectionSpeed(CirclePtsoutlet(:,:,2), CirclePtsoutlet(:,:,3), ...
                         CrossSectionSpeedoutlet, CirclePtsoutlet(1,1,1));
caxis([MinCaxis MaxCaxis])
colormap(jet)
cbar = colorbar('south');
set(cbar,'FontSize',12)
set(get(cbar,'title'),'string','Speed (mm/s)');
set(gca,'FontName','cambria math','FontSize',12)
axis off
view([45 10])
rotate3d
axis equal
box off

if verbose_PrintProfiles
    print('SpeedProfileInsideVessel','-dpng','-r0')
end

%% Change to the side view 
figure('Color','white')
hold on
Patch_Mesh(coord, connect, 0.1)
% LineSpeed
ContourLineSpeed(LinePts(:,:,1), LinePts(:,:,2), LineSpeed);
% Middle section speed
ContourCrossSectionSpeed(CirclePtsmid(:,:,2), CirclePtsmid(:,:,3), ...
                         CrossSectionSpeedmid, CirclePtsmid(1,1,1));
% inlet section speed
ContourCrossSectionSpeed(CirclePtsinlet(:,:,2), CirclePtsinlet(:,:,3), ...
                         CrossSectionSpeedinlet, CirclePtsinlet(1,1,1));
% outlet section speed
ContourCrossSectionSpeed(CirclePtsoutlet(:,:,2), CirclePtsoutlet(:,:,3), ...
                         CrossSectionSpeedoutlet, CirclePtsoutlet(1,1,1));
caxis([MinCaxis MaxCaxis])
colormap(jet)
cbar = colorbar('southoutside');
set(cbar,'FontSize',12)
set(get(cbar,'title'),'string','Speed (mm/s)');
set(gca,'FontName','cambria math','FontSize',12)
axis off
view([0 90])
rotate3d
axis equal
box off

if verbose_PrintProfiles
    print('SpeedProfileInsideVesselSideView','-dpng','-r0')
end

%% Depict VelocityInsideDomain
figure('Color','white')
hold on
Patch_Mesh(coord, connect, 0.1)
quiver3(EvaluationPts(1,1:numel(LinePtsX)), ...
        EvaluationPts(2,1:numel(LinePtsX)), ...
        EvaluationPts(3,1:numel(LinePtsX)), ...
        VelocityInsideDomain(1,1:numel(LinePtsX)), ...
        VelocityInsideDomain(2,1:numel(LinePtsX)), ...
        VelocityInsideDomain(3,1:numel(LinePtsX)), ...
        0, 'b', 'LineWidth', 1.1)
axis off
view([0 90])
rotate3d
axis equal
box off

if verbose_PrintProfiles
    print('VelocityProfileInsideVesselSideView','-dpng','-r0')
end

%% Plot speed along the centerline of vessel
figure('Color','white')
plot(LinePts(11,:,1), LineSpeed(11,:), ...
     'k-o','MarkerSize',0.25*MarkerSizeind,'LineWidth',LineWidthind)
ylabel('Speed (mm/s)')
xlabel('The centerline along the vessel ($\mu$m)')
xlim([min(LinePts(11,:,1)) max(LinePts(11,:,1))])
set(gca,'FontName','cambria math','FontSize',12)
% axis equal
box off

if verbose_PrintProfiles
    print('SpeedCenterlineVessel','-dpng','-r0')
end