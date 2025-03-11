function txt = cursorXYmap_dmm(~,info)
% the function format data to be displayed as the cursor moves on the
% plot
%----------------------------------------------------------------------
x = info.Position(1);
y = info.Position(2);
z = info.Target(1);
% display current value on the dynamic window
txt = ['X  = ',num2str(x,'%.2f'),' [mm], Y  = ',num2str(y,'%.2f'),' [mm]'];
end