function[hxy] = plot_deviationSii(gridSiiX,gridSiiY,Sii_dev,mapsize,diffmap)
% the function calculates the stress/strain-rate distribution over an
% indexed image and plot the results as 2D color-coded map
%--------------------------------------------------------------------------
hxy = axes('InnerPosition',[0.25 0.042 0.46 0.825],'PositionConstraint','innerposition');
xlabel('distance [mm]');
ylabel('distance [mm]');
hcb = colorbar(hxy);
set(hxy,'Colormap',diffmap.map);
caxis([-100 +100]);
hold on
set(hxy,'TickDir','out','XAxisLocation','top');
set(hxy,'XLim',[0 mapsize(2)]);  
set(hxy,'YLim',[0 mapsize(1)]);
set(hcb,'Visible','on','Location','SouthOutSide','TickDirection','out');
pcolor(hxy,gridSiiX,gridSiiY,Sii_dev);
shading interp
axis ij image
box on
set(gca,'FontSize',8);
end