function[hcmenu] = context_menu
% the function creates context menu that allow to change the aspect of
% graphic objects plotted within the mainplot
%--------------------------------------------------------------------------
hcmenu = uicontextmenu;

% Define callbacks 
hdel = 'set(gco,''LineStyle'',''none'')';
hdashed = 'set(gco,''LineStyle'',''--'')';
hdotted = 'set(gco,''LineStyle'','':'')';
hsolid = 'set(gco,''LineStyle'',''-'')';
hthin = 'set(gco,''LineWidth'',1)';
hmid = 'set(gco,''LineWidth'',2)';
hthick = 'set(gco,''LineWidth'',3)';
hblack = 'set(gco,''Color'',''k'')';
hblue = 'set(gco,''Color'',''b'')';
hred = 'set(gco,''Color'',''r'')';
hgreen = 'set(gco,''Color'',''g'')';
hcyan = 'set(gco,''Color'',''c'')';
hmagenta = 'set(gco,''Color'',''m'')';

% Define context menu items and install callbacks
uimenu(hcmenu,'Label','Delete','Callback',hdel);
uimenu(hcmenu,'Label','Style:');
uimenu(hcmenu,'Label','---','Callback',hdashed);
uimenu(hcmenu,'Label','...','Callback',hdotted);
uimenu(hcmenu,'Label','___','Callback',hsolid);
uimenu(hcmenu,'Label','Size:','Separator','on');
uimenu(hcmenu,'Label','thin','Callback',hthin);
uimenu(hcmenu,'Label','medium','Callback',hmid);
uimenu(hcmenu,'Label','thick','Callback',hthick);
uimenu(hcmenu,'Label','Color:','Separator','on');
uimenu(hcmenu,'Label','black','Callback',hblack);
uimenu(hcmenu,'Label','blue','Callback',hblue);
uimenu(hcmenu,'Label','red','Callback',hred);
uimenu(hcmenu,'Label','green','Callback',hgreen);
uimenu(hcmenu,'Label','cyan','Callback',hcyan);
uimenu(hcmenu,'Label','magenta','Callback',hmagenta);

% Define context menu for MARKERS
hcmenup = uicontextmenu;
% Define callbacks for context menu 
mdel = 'set(gco,''Marker'',''none'')';
mshape0 = 'set(gco,''Marker'',''o'')';
mshape1 = 'set(gco,''Marker'',''s'')';
mshape2 = 'set(gco,''Marker'',''d'')';
mshape3 = 'set(gco,''Marker'',''p'')';
mshape4 = 'set(gco,''Marker'',''h'')';
msize0 = 'set(gco,''MarkerSize'',2)';
msize1 = 'set(gco,''MarkerSize'',3)';
msize2 = 'set(gco,''MarkerSize'',5)';
msize3 = 'set(gco,''MarkerSize'',7)';
msize4 = 'set(gco,''MarkerSize'',9)';
msize5 = 'set(gco,''MarkerSize'',11)';
msize6 = 'set(gco,''MarkerSize'',13)';
msize7 = 'set(gco,''MarkerSize'',15)';
mcolor1 = 'set(gco,''MarkerFaceColor'',''w'')';
mcolor2 = 'set(gco,''MarkerFaceColor'',''b'')';
mcolor3 = 'set(gco,''MarkerFaceColor'',''g'')';
mcolor4 = 'set(gco,''MarkerFaceColor'',''y'')';
mcolor5 = 'set(gco,''MarkerFaceColor'',''c'')';
mcolor6 = 'set(gco,''MarkerFaceColor'',''m'')';
% Define context menu items and install callbacks
uimenu(hcmenup,'Label','Delete','Callback',mdel);
uimenu(hcmenup,'Label','Shape:');
uimenu(hcmenup,'Label','circle','Callback',mshape0);
uimenu(hcmenup,'Label','square','Callback',mshape1);
uimenu(hcmenup,'Label','diamond','Callback',mshape2);
uimenu(hcmenup,'Label','pentagon','Callback',mshape3);
uimenu(hcmenup,'Label','hexagon','Callback',mshape4);
uimenu(hcmenup,'Label','Size:','Separator','on');
uimenu(hcmenup,'Label','1','Callback',msize0);
uimenu(hcmenup,'Label','2','Callback',msize1);
uimenu(hcmenup,'Label','3','Callback',msize2);
uimenu(hcmenup,'Label','4','Callback',msize3);
uimenu(hcmenup,'Label','5','Callback',msize4);
uimenu(hcmenup,'Label','6','Callback',msize5);
uimenu(hcmenup,'Label','7','Callback',msize6);
uimenu(hcmenup,'Label','8','Callback',msize7);
uimenu(hcmenup,'Label','Color:','Separator','on');
uimenu(hcmenup,'Label','white','Callback',mcolor1);
uimenu(hcmenup,'Label','blue','Callback',mcolor2);
uimenu(hcmenup,'Label','green','Callback',mcolor3);
uimenu(hcmenup,'Label','yellow','Callback',mcolor4);
uimenu(hcmenup,'Label','cyan','Callback',mcolor5);
uimenu(hcmenup,'Label','magenta','Callback',mcolor6);
end