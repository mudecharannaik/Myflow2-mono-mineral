function[hxy] = plot_dmm(AX,Daxis,mech,mechcmap,T,type)
% the function calculates the stress/strain-rate distribution over an
% indexed image and plot the results as 2D color-coded map
%--------------------------------------------------------------------------
pos_ax = [0.25 0.042 0.46 0.825];
hxy = axes('InnerPosition',pos_ax,'PositionConstraint','innerposition');
switch type
    case 'voigt'
        % isostrain condition
        xlabel([char(963),' [MPa]']);
    case 'reuss'
        % isostress condition
        xlabel(['Log_{10}(',char(941),') [s^{-1}]']);
end      
ylabel('log_{10}(d) [mm]');
hold on
xticks('auto'); yticks('auto');
set(hxy,'TickDir','out','XAxisLocation','top');
hcb = colorbar(hxy);
set(hcb,'Visible','off','Location','SouthOutSide','TickDirection','out');
[X,Y] = meshgrid(AX,Daxis); 
pcolor(hxy,X,Y,mech);
set(hxy,'YScale','log')
set(hxy,'XScale','log');
xlim([AX(1) AX(end)]);
axis ij
shading interp
box on
colormap(mechcmap);
Le = (AX(end)-AX(1)); 
He = (Daxis(end))-(Daxis(1)); 
maxH = (Daxis(end));
minL = AX(1); 
pos_rec = [minL+0.08*Le maxH-0.9*He 0.8*Le 0.5*He];
pos_dc = [minL+0.085*Le maxH-0.88*He 0.07*Le 0.03*He];
pos_diff = [minL+0.085*Le maxH-0.81*He 0.07*Le 0.05*He];
pos_gbs = [minL+0.085*Le maxH-0.7*He 0.07*Le 0.07*He];
% show legend
rectangle(hxy,'Position',pos_rec,...
    'Curvature',0.1,...
    'FaceColor','w',...
    'EdgeColor','k',...
    'LineWidth',0.1);
rectangle(hxy,'Position',pos_dc,...
    'FaceColor',mechcmap(1,:),...
    'EdgeColor','k',...
    'LineWidth',0.1);
text(minL+0.17*Le,maxH-0.87*He,'Dislocation creep');
rectangle(hxy,'Position',pos_diff,...
    'FaceColor',mechcmap(2,:),...
    'EdgeColor','k',...
    'LineWidth',0.1);
text(minL+0.17*Le,maxH-0.79*He,'Diffusion creep');
rectangle(hxy,'Position',pos_gbs,...
    'FaceColor',mechcmap(3,:),...
    'EdgeColor','k',...
    'LineWidth',0.1);
text(minL+0.17*Le,maxH-0.68*He,'GBS');
text(minL+0.088*Le,maxH-0.55*He,['temperature = ',num2str(T-273.15,'%.1f'),'°C']);
end