function txt = cursorXYmap_deviationSii(~,info)
% the function format data to be displayed as the cursor moves on the
% plot
%----------------------------------------------------------------------
x = info.Position(1);
y = info.Position(2);
z = info.Target(1);
xgrid = z.XData;
ygrid = z.YData;
val = z.CData;
% interpolate to find the value of current cursor position
cursor_val = interp2(xgrid,ygrid,val,x,y);
% display current value on the dynamic window
txt = ['X  = ',num2str(x,'%.1f'),' [mm], Y  = ',num2str(y,'%.1f'),' [mm], dev. from average Sii = ',num2str(cursor_val),'%.1f',' [%]'];
end