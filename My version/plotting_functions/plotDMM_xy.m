function[hdmm] = plotDMM_xy(w,isize,bound,grain_D,grain_ID,eff_mech_xy,mapsize,mechcmap)
% the function calculates the distribution of stable deformation mechanisms over an
% indexed image and plot the results as 2D color-coded map
%--------------------------------------------------------------------------
MECH_map = zeros(isize(1),isize(2));  
for i = 1:numel(grain_D)
    MECH_map(grain_ID{i}) = eff_mech_xy(i);
end
% hide grain-boundaries
for i = 1:isize(1)*isize(2)
    if bound(i)==0
        MECH_map(i) = nan;
    end
end
% plot phasemaps and deformation mechanism map
SiiX = linspace(0,w,isize(2));
SiiY = linspace(0,(w/isize(2))*isize(1),isize(1));
[gridSiiX,gridSiiY] = meshgrid(SiiX,SiiY);

% plot deformation mechanism distribution map in the second figure
hdmm = axes('InnerPosition',[0.25 0.042 0.46 0.825],'PositionConstraint','innerposition');
set(gca,'FontSize',8);
xlabel('distance [mm]');
ylabel('distance [mm]');
hcb = colorbar(hdmm);
%-----------------------------
% evaluate active deformation mechanism and adjust the colormap if needed
defm = unique(eff_mech_xy);
range = min(defm):1:max(defm);
if numel(range)<3
    if numel(range)==1
        switch range
            case 1
                % dislocation creep only
                mechcmap = mechcmap(1,:);
            case 2
                % diffusion creep only
                mechcmap = mechcmap(2,:);
            case 3
                % gbs only
                mechcmap = mechcmap(3,:);
        end
    else
        % current map has only two deformation mechanisms
        if min(defm)==1 && max(defm)==2
            % dislocation +diffusion creep 
            mechcmap = mechcmap(1:2,:);
        elseif min(defm)==2 && max(defm)==3
            % diffusion creep + gbs
            mechcmap = mechcmap(2:3,:);
        else
            % dislocation creep + gbs
            mechcmap = mechcmap([1,3],:);
        end
    end
end   
%---------------------
colormap(mechcmap);
hold on
set(hdmm,'TickDir','out','XAxisLocation','top');
set(hdmm,'XLim',[0 mapsize(2)]);  
set(hdmm,'YLim',[0 mapsize(1)]);
set(hcb,'Visible','off','Location','SouthOutSide','TickDirection','out');
pcolor(hdmm,gridSiiX,gridSiiY,MECH_map);
shading interp
axis ij image
box on
set(gca,'FontSize',8);
end