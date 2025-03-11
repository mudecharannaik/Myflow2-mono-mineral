function[hxy] = plot1d(lab1,lab2,   var,        vary,varx,mech,Xl,Yl,model,varargin)
% the function calculates the stress/strain-rate distribution over an
% indexed image and plot the results as 2D color-coded map
%--------------------------------------------------------------------------
pos_ax = [0.25 0.042 0.46 0.825];
hxy = axes('InnerPosition',pos_ax,'PositionConstraint','innerposition');
xlabel(' ');ylabel(' ');
xlabel(Xl);ylabel(Yl);
hold on
xticks(hxy,[]);yticks(hxy,[]);
xticks('auto'); yticks('auto');
set(hxy,'TickDir','out','XAxisLocation','top');
Xr = [min(varx) max(varx)]; 
Yr = [min(vary) max(vary)];
hcb = colorbar(hxy);
set(hcb,'Visible','off','Location','SouthOutSide','TickDirection','out');
switch nargin 
    case 10
        % evaluate one scale (to go logarithmic)
        axdir = varargin{1};
        if strcmp(axdir,'logY')==1
            vary = log10(vary);
            Yr = log10(Yr);
        else
            var = log10(var);
            Xr = log10(Xr);
        end
    case 11
    % change both scales to logarithmic  
    vary = log10(vary);
    Yr = log10(Yr);
    var = log10(var);
    Xr = log10(Xr);
end
set(hxy,'YLim',Yr)
set(hxy,'XLim',Xr) 
axis ij 
box on
set(gca,'FontSize',8);
if model==0 
    % monophase models show different deformation mechanisms
    plot(hxy,var(mech==1),vary(mech==1),'-b');
    plot(hxy,var(mech==2),vary(mech==2),'-r');
    plot(hxy,var(mech==3),vary(mech==3),'-k');
    % add labels to the plot (monophase models only)
    Le = (Xr(2)-Xr(1)); He = (Yr(2)-Yr(1)); minH = Yr(1); minL = Xr(1); 
    pos_rec = [minL+0.03*Le minH+0.05*He 0.3*Le 0.18*He];
    pos_dc = [minL+0.05*Le minH+0.06*He 0.03*Le 0.02*He];
    pos_diff = [minL+0.05*Le minH+0.09*He 0.03*Le 0.02*He];
    pos_gbs = [minL+0.05*Le minH+0.12*He 0.03*Le 0.02*He];
    rectangle(hxy,'Position',pos_rec,...
        'Curvature',[0.1 0.1],...
        'FaceColor','w',...
        'EdgeColor','k',...
        'LineWidth',0.1);
    rectangle(hxy,'Position',pos_dc,...
        'FaceColor','b',...
        'EdgeColor','k',...
        'LineWidth',0.1);
    text(minL+0.1*Le,minH+0.07*He,'Dislocation creep');
    rectangle(hxy,'Position',pos_diff,...
        'FaceColor','r',...
        'EdgeColor','k',...
        'LineWidth',0.1);
    text(minL+0.1*Le,minH+0.1*He,'Diffusion creep');
    rectangle(hxy,'Position',pos_gbs,...
        'FaceColor','k',...
        'EdgeColor','k',...
        'LineWidth',0.1);
    text(minL+0.1*Le,minH+0.13*He,'grain boundary sliding');
    text(minL+0.05*Le,minH+0.17*He,' ')
    text(minL+0.05*Le,minH+0.2*He,' ');
    text(minL+0.05*Le,minH+0.17*He,lab1);
    text(minL+0.05*Le,minH+0.2*He,lab2);
else
    % composite do not show variations of deformation mechanisms in 1d
    % plots
    plot(hxy,var,vary,'-m','LineWidth',2);
end
end