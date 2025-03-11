function[h2]=initialize_plot
% the function initialize_plot.m creates an empty plot
%--------------------------------------------------------------------------
h2 = axes('InnerPosition',[0.25 0.042 0.46 0.825],'PositionConstraint','innerposition','NextPlot','replacechildren');
x_label = 'variable (x)';
y_label = 'variable (y)';
hcb = colorbar(h2);
colormap(parula(128));
x_range = [0.1 600];  
y_range = [273.15 1473.15]; 
set(h2,'TickDir','out','XAxisLocation','top');
set(h2,'XLim',x_range);  
set(h2,'YLim',y_range);
set(hcb,'Visible','off','Location','SouthOutSide','TickDirection','out');
set(get(gca,'XLabel'),'string',x_label);
set(get(gca,'YLabel'),'string',y_label);
box on
set(gca,'FontSize',8);
axis(h2,'ij'); 
end