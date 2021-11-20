function Plot_LineFieldsOnVessel(coord, NormField, WallLine, ...
                                 MarkerSizeind, LineWidthind, ...
                                 ylabelText, xlabelText)

figure('Color','white')
hold on
plot(coord(1,WallLine), NormField(WallLine), ...
     'k-o','MarkerSize',0.25*MarkerSizeind,'LineWidth',LineWidthind)
ylabel(ylabelText)
xlabel(xlabelText)
xlim([min(coord(1,WallLine)) max(coord(1,WallLine))])
set(gca,'FontName','cambria math','FontSize',12)