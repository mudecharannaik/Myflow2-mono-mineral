function MYflow2
% The software allow calculating the flow strength of composite materials 
% from end-member flow laws using several mixture models. MYflow also 
% compute and displays deformation mechanism maps of single-phase materials, 
% and paleopiezometric data.
% 
% developers: Leonardo Casini & Matteo Maino
% programmed by L. Casini, University of Sassari, 07100 Sassari (Italy)
% contact: casini@uniss.it
%-------------------------------------------------------------------------
% DEBUG HISTORY
% 2023-09-19 fixed wrong call to function plot1d.m in Eii-T maps
% 2023-09-19 modified sliders for setting grainsize, strainrate, stress 
% 2023-09-16 wrong call to mehl & hirth 2008 piezometer (calcPP function)
% 2023-05-02 fixed bug in the grainsize.m function (wrong formula for
%            calculating the partial area fractions of mineral phases)
% 2023-05-01 fixed bug in Huet et al (2014) mixing rule
% 2023-04-11 fixed rounding error in the expression for median
%            stress/strain rate
% 2023-04-06 fixed error in text labels of median strain-rate in Eii_xy
%            maps
% 2023-03-09 fixed error in mixing rules (Sii_T_voigt.m function)
% 2023-02-07 corrected wrong phase indexing of monophase rocks in 
%            'grainsize.m' function.
% 2023-01-30 modified function median_dev_Sii - now the median stress value
%            is computed from the cumulative sum of grain stresses 
%            weighted by the molar volume proportion of phases in the composite. 
% 2023-01-30 modifed plotting functions for Sii_xy and Eii_xy maps, now the
%            grains of non identified material are excluded from the plot
% 2023-01-30 fixed handle visibility issue of figure 'grain size statistics'
% 2023-01-28 fixed bug in the voigt, reuss functions - now both mechanical
%            constraints are enabled after previous selection is cancelled
% 2023-01-26 modified grainsize statistics in the function 'grainsize' -
%            now, each phase has its own statistics 
% 2022-11-16 added references curves of viscous creep of quartzite and
%            olivine (menu 'Plot')
% 2022-11-14 fixed unwanted change of radiobutton group properties and position
%            after 'Reset_rheology' callback execution
% 2022-11-09 modified behavior of paleopiezometry table - now data can be
%            displayed independently of the model type selected (if any), &
%            the whole dataset can mix subsets of different mineral phases  
% 2022-11-06 modified function caclPP
% 2022-11-05 added database of piezometric calibrations
% 2022-10-30 added stress/strainrate deviations maps from a median value.
%            For the moment, the median value is used, to be replaced with
%            a median value weighted by molar proportion
% 2022-10-27 added dynamic datatip to display stress/strainrate value
%            (in xy-maps) on cursor click 
% 2022-10-26 added dynamic text to display current value of intensive
%            variables modified by action on sliders
% 2022-10-26 set fixed ranges of stress/strainrate in Sii/Eii-xy colormaps
% 2022-10-25 fixed T settings (now appears as °C in the GUI and is
%            thoroughly computed as K degree by all rheological functions)
% 2022-10-25 fixed dimensional problem in flow law of Sii/Eii-xy maps
% 2022-10-24 corrected re-initialization of the main plot after resetting
% 2022-10-21 modified calculation of mean stress/strainrate for composite
%            based on different mixing rules - now composite strength is
%            evaluated based on selected mechanical constraint, mixing rule
%            and projection space (involves changes of functions
%            Sii_T_voigt.m, Eii_T_voigt.m)
% 2022-10-20 fixed reset_rheology.m function
% 2022-10-18 fixed wrong handle numbering of main plot
% 2022-10-17 added slider to change dynamically grainsize in Sii/T and
%            Eii/T plots
% 2022-10-14 added sliders to allow modifying dynamically T,Eii and Sii 
% 2022-10-13 changed order of message boxes appearing while setting up
%            compositional maps
% 2022-09-39 added radiobutton that allow switching between different plots
%            (Sii-Eii/xy maps only) 
% 2022-09-27 fixed bug in Sii_xy and Eii_xy functions (maps were plotted 
%            without requesting T/Sii-Eii conditions)
% 2022-09-27 fixed unwanted behavior of 'reset rheology' menu  
% 2022-09-27 the deformation mechanism map now shows the correct labels 
% 2022-09-25 function calcPP.m unwanted change of the color of latest paleostress 
%            when asking for non-existent piezometric calibration
% 2022-09-25 fixed wrong behaviour of function 'calcPP.m' - piezometric 
%            data were not shown while selecting 'add single grain size'  
%==========================================================================
SZ           = get(0,'ScreenSize');
pos_f        = [5,5,SZ(3)-10,SZ(4)-200];
f            = figure('Visible','off','Position',pos_f,'Resize','off');
% set path to current folder
currentFolder = pwd;
warning('off','all');
addpath(genpath(currentFolder));

% set colors and colormaps
inact_txt    = [0.6 0.6 0.6];           % inactive text
active_txt   = [0.3 0.3 0.3];           % active text
chkcol       = [0.94 0.94 0.94];        % background of text boxes
dccol   = [0.1137    0.7725    0.8863]; % blue
diffcol = [0.7333    0.4235    0.8784]; % violet
gbscol  = [0.9647    0.3333    0.4196]; % red
mechcmap  = [dccol; diffcol; gbscol];
cmapdiff = load('diffmap.mat');
siimap = load('siimap.mat');

%% DATABASES

% extract mineral database
listm  = dir('DB_mineral_parameters\*.mat');
[r,~]  = size(listm);
listDB = cell(r,1);
for mineralID = 1:r
    namewithext       = listm(mineralID).name;
    namewithoutext    = strsplit(namewithext,'.');
    currname          = namewithoutext(1);
    listDB(mineralID) = currname;
end

% extract composite database
listc  = dir('DB_composite\*.mat');
[r1,~] = size(listc);
listCO = cell(r1,1);
for compID = 1:r1
    cnamewithext    = listc(compID).name;
    cnamewithoutext = strsplit(cnamewithext,'.');
    ccurrname       = cnamewithoutext(1);
    listCO(compID)  = ccurrname;
end

% extract the piezometric calibration database 
listCAL = dir('DB_calibrations\*.mat');
[r2,~]  = size(listCAL);
listcal = cell(r2,1);
mincal = cell(1,r2);
indV = 1:r2;
for calID = 1:r2
    cnamewithext    = listCAL(calID).name;
    cnamewithoutext = strsplit(cnamewithext,'.');
    ccurrname       = cnamewithoutext(1);
    listcal(calID)   = ccurrname;
    calpar = load(['DB_calibrations\',cnamewithext]);
    mincal(calID) = calpar.par(1);
end

%% (0)    - Mineral Database
m0  = uimenu(f,'Label','DB materials');
m0h = gobjects(5,1);
m0h(1) = uimenu(m0,'Label','Edit mineral DB');    
         % 1 sub-menu, single minerals DB
         m0h(2) = uimenu(m0h(1),'Label','Add new phase','Callback',{@new_phase});
         m0h(3) = uimenu(m0h(1),'Label','Edit existent phase','Callback',{@edit_phase});
         m0h(4) = uimenu(m0h(1),'Label','Delete phase','Callback',{@del_phase});
m0h(5) = uimenu(m0,'Label','Edit composites DB');
         % 2 sub-menu, composites
       uimenu(m0h(5),'Label','NEW composite','Callback',{@new_comp});
       uimenu(m0h(5),'Label','Delete composite','Callback',{@del_comp});

%% (1)    - select type of model as single phase/composite
m1   = uimenu(f,'Label','Rheological Model'); 
m1h  = gobjects(22,1);
     uimenu(m1,'Label','Model type');
     % I sub-menu, model type (single phase/composite) & mechanical
     % constraint
     m1h(1) = uimenu(m1,'Label','Single mineral'); 
            m1h(2) = uimenu(m1h(1),'Label','Load mineral','Callback',{@select_phase});
     m1h(3) = uimenu(m1,'Label','Composite');
            m1h(4) = uimenu(m1h(3),'Label','Load composite','Callback',{@load_comp});
     m1h(5) = uimenu(m1,'Label','Mechanical constraints','Separator','on');
            m1h(6) = uimenu(m1h(5),'Label','Iso strain (Voigt,1928)','Callback',{@voigt});
            m1h(7) = uimenu(m1h(5),'Label','Iso stress (Reuss, 1929)','Callback',{@reuss});
     m1h(8) = uimenu(m1,'Label','Rule of mixture','Separator','on');
            m1h(9)= uimenu(m1h(8),'Label','Geometric mean (Matthies & Humber, 1993)','Callback',{@geomean});
            m1h(10) = uimenu(m1h(8),'Label','Arithmetic mean (Voigt, 1928)','Callback',{@arithmetic});
            m1h(11) = uimenu(m1h(8),'Label','Harmonic mean (Reuss, 1929)','Callback',{@harmonic});
            m1h(12) = uimenu(m1h(8),'Label','Minimized Power Geometric mean (Huet et al., 2014)','Callback',{@huet});
            m1h(13) = uimenu(m1h(8),'Label','Thermodynamic mixing (Hobbs et al., 2019)','Callback',{@hobbs});
     m1h(14) = uimenu(m1,'Label','Projection Space','Separator','on');
            m1h(15) = uimenu(m1h(14),'Label','Stress/T','Callback',{@Sii_T});
            m1h(16) = uimenu(m1h(14),'Label','Strain-rate/T','Callback',{@Eii_T});
            m1h(17) = uimenu(m1h(14),'Label','Stress/grain-size','Callback',{@Sii_d});    
            m1h(18) = uimenu(m1h(14),'Label','Strain-rate/grain-size','Callback',{@Eii_d}); 
            m1h(19) = uimenu(m1h(14),'Label','Stress/x,y','Callback',{@Sii_xy},'Separator','on');
            m1h(20) = uimenu(m1h(14),'Label','Strain-rate/x,y','Callback',{@Eii_xy});
            m1h(21) = uimenu(m1h(14),'Label','Deformation Mechanism Map','Callback',{@dmm},'Separator','on');
     m1h(22) = uimenu(m1,'Label','RESET model','Separator','on','Callback',{@reset_rheology});
     m1h(23) = uimenu(m1,'Label','RESET mineral','Callback',{@reset_mineral});
set(m1h(6:7),'Enable','off');   % Voigt/Reuss
set(m1h(9:13),'Enable','off');  % mixture rule 7=median,8=aritmetic,9=harmonic,10=huet,11=hobbs
set(m1h(15:23),'Enable','off'); % projection spaces 13=Sii/T,14=Eii/T,15=Sii/d,16=Eii/d,17=Sii/xy.18=Eii/xy,19=dmm
 
%% (2)  PLOTTING COMMANDS
m2      = uimenu(f,'Label','Plot');
m2h     = gobjects(7,1);
        uimenu(m2,'Label','Reference envelopes');
        % I sub-menu, Byerlee envelope and other reference curves
        m2h(1) = uimenu(m2,'Label','Byerlee normal fault','Callback',{@byen});
        m2h(2) = uimenu(m2,'Label','Byerlee strike-slip fault','Callback',{@byess});
        m2h(3) = uimenu(m2,'Label','Byerlee reverse fault','Callback',{@byer});
        m2h(4) = uimenu(m2,'Label','Quartzite disl.creep','Callback',{@qtzcreep});
        m2h(5) = uimenu(m2,'Label','Olivine disl.creep','Callback',{@olcreep});
        uimenu(m2,'Label','Display data','Separator','on');
        % II sub-menu, plot data
        m2h(6) = uimenu(m2,'Label','Export figure','Callback',{@export});
        m2h(7) = uimenu(m2,'Label','RESET plot ','Callback',{@cancel_plot},'Separator','on');
set(m2h(6:7),'Enable','off');

%% (3)  paleopiezometry
m3h     = gobjects(5,1);
m3      = uimenu(f,'Label','Paleopiezometry');
        uimenu(m3,'Label','Set new calibration','Callback',{@new_cal});
        uimenu(m3,'Label','Modify calibration','Callback',{@mod_cal});
        m3h(1) = uimenu(m3,'Label','Import grain-size','Callback',{@d_imp},'Separator','on');
        m3h(2) = uimenu(m3,'Label','Add single data','Callback',{@d_sing});        
        m3h(3) = uimenu(m3,'Label','CANCEL dataset','Callback',{@cancel_pp},'Enable','off');
        m3h(4) = uimenu(m3,'Label','calculate','Callback',{@calcPP},'Enable','off','Separator','on');
        m3h(5) = uimenu(m3,'Label','Hide data','Callback',{@hidePP},'Enable','off');
%% (4)  HELP
        uimenu(f,'Label','HELP','Callback',{@help});
%% OTHER GUI COMPONENTS---------------------------------------------------------------------------
% datazone left of the main plot 
datapan = uicontrol('Style','frame',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[0.0115*SZ(3) SZ(4)*0.379 SZ(3)*0.197 SZ(4)*0.295],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
dzone_text = uicontrol('Style','text',...
    'String','MODEL PARAMETERS ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.02 SZ(4)*0.64 SZ(3)*0.1042 SZ(4)*0.0243],...   
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'FontWeight','bold');
modtxt1          = uicontrol('Style','text',...
    'String','Model type: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.03 SZ(4)*0.61 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
mdtxt1_sel      = uicontrol('Style','text',...
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.61 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);
modname      = uicontrol('Style','text',...
    'String','Material name: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.024 SZ(4)*0.59 SZ(3)*0.07 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
modname_sel      = uicontrol('Style','text',...
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.59 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);
modbehavior = uicontrol('Style','text',...
    'String','Mechanical constraint:',...           
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.012 SZ(4)*0.57 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
modbehavior_sel = uicontrol('Style','text',...
    'String','  ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.57 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);
mixname      = uicontrol('Style','text',...
    'String','Mixture rule: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.0295 SZ(4)*0.55 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
mixname_sel    = uicontrol('Style','text',...
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.55 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);
piezo_name   = uicontrol('Style','text',...
    'String','Piezometer: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.0498 SZ(4)*0.53 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8,'FontWeight','bold');
piezo_name_sel    = uicontrol('Style','text',...    
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.53 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);  
cal_name     = uicontrol('Style','text',...
    'String','Calibration: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.029 SZ(4)*0.51 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,'FontWeight','bold');
cal_name_sel     = uicontrol('Style','text',...
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.51 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);
pmtxt1          = uicontrol('Style','text',...
    'String','Projection space: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.029 SZ(4)*0.49 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
pmtxt1_sel      = uicontrol('Style','text',...
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.49 SZ(3)*0.08 SZ(4)*0.0243],...
    'HorizontalAlignment','left',...
    'FontSize',8);
pmphase      = uicontrol('Style','text',...
    'String','minerals: ',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.029 SZ(4)*0.47 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
pmphase_sel      = uicontrol('Style','text',...
    'String',' ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.43 SZ(3)*0.08 SZ(4)*0.0643],...
    'HorizontalAlignment','left',...
    'FontSize',8);
% create sliders 
slidepan = uicontrol('Style','frame',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[0.0115*SZ(3) SZ(4)*0.15 SZ(3)*0.197 SZ(4)*0.228],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
slid_text = uicontrol('Style','text',...
    'String','VARIABLES ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.02 SZ(4)*0.344 SZ(3)*0.1042 SZ(4)*0.0243],...  
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'FontWeight','bold');
tempsld3 = uicontrol('Style','text',...
    'String','T [°C]:',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.019 SZ(4)*0.295 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
tempsld2 = uicontrol('Style','text',...
    'String','0                                        1500',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.31 SZ(3)*0.1 SZ(4)*0.0243],...
    'HorizontalAlignment','center',...
    'FontSize',7,...
    'FontWeight','bold');
tempsld1 = uicontrol('Style','slider',...
    'Max',1500,...
    'Min',0,...
    'SliderStep',[5/1500,50/1500],...
    'BackGroundColor',chkcol,...
    'Enable','off',...
    'Position',[SZ(3)*0.1 SZ(4)*0.3 SZ(3)*0.0842 SZ(4)*0.02],...
    'Callback',{@slidT},...
    'value',500);
% stress
stresssld3 = uicontrol('Style','text',...
    'String','Stress [MPa]:',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.019 SZ(4)*0.255 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
stresssld2 = uicontrol('Style','text',...
    'String','0                                         500',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.095 SZ(4)*0.27 SZ(3)*0.1 SZ(4)*0.0243],...
    'HorizontalAlignment','center',...
    'FontSize',7,...
    'FontWeight','bold');
stresssld1 = uicontrol('Style','slider',...
    'Max',500,...
    'Min',0,...
    'SliderStep',[5/500,25/500],...
    'BackGroundColor',chkcol,...
    'Enable','off',...
    'Position',[SZ(3)*0.1 SZ(4)*0.26 SZ(3)*0.0842 SZ(4)*0.02],...
    'Callback',{@slidSii},...
    'value',100);
% strainrate
strainsld3 = uicontrol('Style','text',...
    'String','Strain-rate [1/s]:',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.019 SZ(4)*0.215 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
strainsld2 = uicontrol('Style','text',...
    'String','1e-15                                    1e-5',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.094 SZ(4)*0.23 SZ(3)*0.1 SZ(4)*0.0243],...
    'HorizontalAlignment','center',...
    'FontSize',7,...
    'FontWeight','bold');
strainsld1 = uicontrol('Style','slider',...
    'Max',-5,...
    'Min',-15,...
    'SliderStep',[0.1/10,1/10],...
    'BackGroundColor',chkcol,...
    'Enable','off',...
    'Position',[SZ(3)*0.1 SZ(4)*0.22 SZ(3)*0.0842 SZ(4)*0.02],...
    'Callback',{@slidEii},...
    'value',-12);
% grainsize
dsld3 = uicontrol('Style','text',...
    'String','grain size [mm]:',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.019 SZ(4)*0.175 SZ(3)*0.064 SZ(4)*0.0243],...
    'HorizontalAlignment','right',...
    'FontSize',8,...
    'FontWeight','bold');
dsld2 = uicontrol('Style','text',...
    'String','0.001                                     10',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.094 SZ(4)*0.19 SZ(3)*0.1 SZ(4)*0.0243],...
    'HorizontalAlignment','center',...
    'FontSize',7,...
    'FontWeight','bold');
dsld1 = uicontrol('Style','slider',...
    'Max',10,...
    'Min',0.001,...
    'SliderStep',[0.005/9.999,0.1/9.999],...
    'BackGroundColor',chkcol,...
    'Enable','off',...
    'Position',[SZ(3)*0.1 SZ(4)*0.18 SZ(3)*0.0842 SZ(4)*0.02],...
    'Callback',{@slidd},...
    'value',0.1);

% creates radiobutton group
bg = uibuttongroup(f,'Position',[0.01 0.055 0.2 0.139],...
    'SelectionChangedFcn',@plotselection); 
radio1 = uicontrol(bg,'Style','radiobutton',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'String','stress/strain-rate map',...
    'Position',[SZ(3)*0.05 SZ(4)*0.045 SZ(3)*0.0842 SZ(4)*0.02],...
    'Enable','off');         
radio2 = uicontrol(bg,'Style','radiobutton',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'String','deformation mechanism map',...
    'Position',[SZ(3)*0.05 SZ(4)*0.025 SZ(3)*0.1 SZ(4)*0.02],...
    'Enable','off');
radio3 = uicontrol(bg,'Style','radiobutton',...
    'ForegroundColor',inact_txt,...
    'BackGroundColor',chkcol,...
    'String','standard deviation map',...
    'Position',[SZ(3)*0.05 SZ(4)*0.005 SZ(3)*0.1 SZ(4)*0.02],...
    'Enable','off');
% creates radiobuttons
radiotxt = uicontrol('Style','text',...
    'String','SWITCH RESULTS ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.02 SZ(4)*0.113 SZ(3)*0.1042 SZ(4)*0.0243],...   
    'HorizontalAlignment','left',...
    'FontSize',10,...
    'FontWeight','bold');

% creates legend for deformation mechanism maps
legend_dmm = gobjects(1,9);
legend_dmm(1) = uicontrol('Style','text','ForegroundColor',[0 0 0],...
    'BackGroundColor',[0 0 0],'String','-----',...
    'Position',[SZ(3)*0.2845 SZ(4)*0.069 SZ(3)*0.041 SZ(4)*0.022],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(2) = uicontrol('Style','text','ForegroundColor',mechcmap(1,:),...
    'BackGroundColor',mechcmap(1,:),'String','-----',...
    'Position',[SZ(3)*0.285 SZ(4)*0.07 SZ(3)*0.04 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(3) = uicontrol('Style','text','String','dislocation creep',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.33 SZ(4)*0.07 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(4) = uicontrol('Style','text','ForegroundColor',[0 0 0],...
    'BackGroundColor',[0 0 0],'String','-----',...
    'Position',[SZ(3)*0.3995 SZ(4)*0.069 SZ(3)*0.041 SZ(4)*0.022],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(5) = uicontrol('Style','text','ForegroundColor',mechcmap(2,:),...
    'BackGroundColor',mechcmap(2,:),'String','-----',...
    'Position',[SZ(3)*0.4 SZ(4)*0.07 SZ(3)*0.04 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(6) = uicontrol('Style','text','String','diffusion creep',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.444 SZ(4)*0.07 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(7) = uicontrol('Style','text','ForegroundColor',[0 0 0],...
    'BackGroundColor',[0 0 0],'String','-----',...
    'Position',[SZ(3)*0.5015 SZ(4)*0.069 SZ(3)*0.041 SZ(4)*0.022],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(8) = uicontrol('Style','text','ForegroundColor',mechcmap(3,:),...
    'BackGroundColor',mechcmap(3,:),'String','-----',...
    'Position',[SZ(3)*0.502 SZ(4)*0.07 SZ(3)*0.04 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_dmm(9) = uicontrol('Style','text','String','grain-boundary sliding',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.548 SZ(4)*0.07 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_Siixy = uicontrol('Style','text','String','MPa',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.650 SZ(4)*0.015 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_Eiixy = uicontrol('Style','text','String',['log(',char(941),')'],...
    'HorizontalAlignment','left','Position',[SZ(3)*0.650 SZ(4)*0.015 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_deviation_Sii = uicontrol('Style','text','String',[char(963),' % deviation'],...
    'HorizontalAlignment','left','Position',[SZ(3)*0.650 SZ(4)*0.015 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_deviation_Eii = uicontrol('Style','text','String',['Log10',char(941),' % deviation'],...
    'HorizontalAlignment','left','Position',[SZ(3)*0.650 SZ(4)*0.015 SZ(3)*0.07 SZ(4)*0.02],...
    'Visible','off',...
    'FontSize',8);
legend_intvar1 = uicontrol('Style','text','String','T = 0.00[°C]',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.25 SZ(4)*0.73 SZ(3)*0.15 SZ(4)*0.015],...
    'Visible','off',...
    'FontSize',8);
legend_intvar2 = uicontrol('Style','text','String','d = 0.00[mm]',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.25 SZ(4)*0.71 SZ(3)*0.15 SZ(4)*0.015],...
    'Visible','off',...
    'FontSize',8);
legend_medianVal = uicontrol('Style','text','String',' ',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.25 SZ(4)*0.69 SZ(3)*0.15 SZ(4)*0.015],...
    'Visible','off',...
    'FontSize',8);
legend_peakVal = uicontrol('Style','text','String','peak value = ',...
    'HorizontalAlignment','left','Position',[SZ(3)*0.25 SZ(4)*0.67 SZ(3)*0.15 SZ(4)*0.015],...
    'Visible','off',...
    'FontSize',8);
% create table for paleopiezometric data
tabtext      = uicontrol('Style','text',...
    'String','Samples ',...
    'ForegroundColor',active_txt,...
    'BackGroundColor',chkcol,...
    'Position',[SZ(3)*0.835 SZ(4)*0.68 SZ(3)*0.1042 SZ(4)*0.0243],...    
    'HorizontalAlignment','left',...
    'FontSize',9,...
    'FontWeight','bold');
colname      = {'ON','mineral','color',['d [',char(956),'m]'],' err.d ','T [°C]','err.T',[char(963),' [MPa]']};
coleditable  = [false,false,false,false,false,false,false,false];
colformat    = {'logical','char', 'numeric','numeric', 'numeric','numeric','numeric'};
dat          = cell(999,8);
postable     = [SZ(3)*0.72,SZ(4)*0.05,SZ(3)*0.26,SZ(4)*0.615];
widthtable   = {30,50,43,46,39,40,39,51};
tab1         = uitable('Parent',f,'Position',postable,...
            'Data',dat,...
            'ColumnName', colname,...
            'ColumnFormat', colformat,...
            'ColumnEditable', coleditable,...          
            'CellEditCallback',@piezo_data,...
            'Visible','on');  
            set(tab1,'ColumnWidth',widthtable);  

%% set context menus (abilitate changing aspect of objects in the plot)
[hcmenu] = context_menu;

%% INITIALIZE VARIABLES----------------------------------------------------
% rheology
mixr                 = '';             % mixture rule
flowp                = cell(3,1);      % flow law parameters
mech                 = [];             % index of nodal effective deformation mechanism 
molp                 = [];             % molar proportions of composite
model                = [];             % composite (1), mono-phase model (0)
VR                   = [];             % initialize choice of Reuss/Voigt model
cgs                  = [];             % median grain-size of phases in composite
comp                 = [];             % name of composite
phase_name           = {};             % name of phases in current model
phase_ID             = [];             % ID of selected phase(s), >1 for composites
cvol                 = [];             % phase proportions in composite
% map variables
T                    = 600;            % mean temperature
Eii                  = 1e-12;          % strainrate
Sii                  = 100;            % stress [MPa]
plottype             = [];             % ID of currently active plot 1=Sii/T, 2=Eii/T, 3=Sii/D, 4=Eii/D, 5=Sii/xy, 6=Eii/xy, 7=DMM
mapgrain_ID          = [];
grain_comp           = [];
mapgrain_D           = [];
mapsize              = [];
wmap                 = [];
mapbound             = [];
eff_mech_xy          = [];
Eii_median           = [];
Sii_median           = [];
Sii_dev              = [];
Eii_dev              = [];
Sii_map              = [];
Eii_map              = [];
Sii_xymap            = [];
Eii_xymap            = [];
SiiEiiMap            = [];
gridSiiX             = [];
gridSiiY             = [];
gridEiiX             = [];
gridEiiY             = [];
% reference curves
bye_n                = 0;              % handle to Byerlee envelope for normal faults
bye_ss               = 0;              % Byerlee envelope for strike-slip faults
bye_r                = 0;              % Byerlee envelope for reverse faults
qtz_dc               = 0;
ol_dc                = 0;
t1                   = 0;
t2                   = 0;
t3                   = 0;
t4                   = 0;
t5                   = 0;
% paleostress data
hmdata               = [];             % handles to paleostress data added to current plot
pstress_mean         = [];             % mean paleostress values computed from grainsize dataset (REQUIRED!)
pstress_min          = [];             % minimum paleostress values (if available)
pstress_max          = [];             % maximum paleostress values (if available)
pTmin                = [];             % minimum temperature estimate associated to grainsize data (if available)
pTmax                = [];             % maximum temperature associated to grainsize data (if available)
datin                = 0;              % counter of uploaded grainsize/T data
datcal               = 0;              % number of piezometric data (to be calculated at the single call of calc_piezo)
datcalALL            = 0;              % number of total calculated piezometric data (before actual call to calc_piezo)
phase_pp             = [];             % ID of phase currently selected for piezometric analysis
calibration          = 0;              % ID to currently selected piezometric calibration

%% PHYSICAL PARAMETERS and AXIS---------------------------------------------
% set physical constraints for all models
% compute pressure profile from linear sectioned geotherm
nXY     = 1000;                                  % number of gridpoints (equal X,Y,Z directions)
Siiaxis = linspace(0.1,400,nXY);                 % stress axis
Eiiaxis = logspace(-15,-6,nXY);                  % strain-rate axis
Daxis   = linspace(0.001,10,nXY);                % grainsize axis [mm]
Taxis   = linspace(273.15,273.15+1200,nXY);      % T axis
P0      = 0.1;                                   % atmospheric pressure at z=0 [MPa]
alpha   = (2700*9.81*1000)/1e6;                  % linear pressure increment/km [MPa]
tgrad   = 30;                                    % linear geotherm [°C/km]
P       = P0+alpha.*(Taxis-273.15)./tgrad;       % linear pressure gradient [GPa]
lambda             = 0.4;                        % pore fluid-pressure factor
byerlee_norm       = (1-lambda)*0.75.*P;         % Byerlee envelope for normal fault 
byerlee_strikeslip = (1-lambda)*1.2.*P;          % Byerlee for strike-slip faults
byerlee_reverse    = (1-lambda)*3.0.*P;          % Byerlee for thrusts

%% set axes (default uses 2d plots)
h2 = initialize_plot;
% enable dynamic data display on moving of the cursor 
dcm = datacursormode(f);
dcm.Enable ='on';
dcm.updateDataCursors;
%% Initialize the GUI.
% Change units to normalized so components resize 
% automatically.
set([f,h2,...
    tabtext,dzone_text,slidepan,slid_text,...
    datapan,radiotxt,radio1,radio2,radio3,...
    legend_dmm,legend_Siixy,legend_Eiixy,legend_deviation_Sii,legend_deviation_Eii,legend_intvar1,legend_intvar2,legend_medianVal,legend_peakVal...
    tempsld1,stresssld1,strainsld1,dsld1,...
    tempsld2,stresssld2,strainsld2,dsld2,...
    tempsld3,stresssld3,strainsld3,dsld3,...
    pmtxt1,pmtxt1_sel,pmphase,pmphase_sel,...
    modtxt1,...
    modname,...
    mixname,...
    mixname_sel,...
    mdtxt1_sel,...
    modname_sel,...
    modbehavior,...
    modbehavior_sel,...
    piezo_name_sel,piezo_name,cal_name,cal_name_sel],'Units','normalized');
% initialize the GUI
set(f,'Name','MYflow2.0')
set(f,'NumberTitle','off')
% Move the GUI to the center of the screen.
movegui(f,'center')               % center the GUI
set(f,'Visible','on');            % Make the GUI visible.


%% (0) EDIT DATABSE---------------------------------------------------------
% this menu allow the following actions to be done:
% 1. create new mineral                 
% 2. delete existent mineral          
% 3. create new composite 
% 4. delete existent composite

function new_phase(~,~)
    % add new phase to current database
    % the function open a modal dialog box that allow the user to
    % upload a new mineral phase by specifying its: 1) molar volume 
    % [kJ/kbar], 2)Burgers vector length [mm], 3) shear modulus [MPa], 
    % 4) pre-exponential constant of dislocation creep [MPa^-n mm^-m],
    % 5) activation energy for lattice diffusion [kJ/mol], 6)
    % activation energy of grain-boundary diffusion [kJ/mol]
    %------------------------------------------------------------------
    name_min = inputdlg('Set the name of new mineral phase being added to the database:','New mineral');
    par = nan(2,3);
    options.Interpreter = 'tex';
    options.Resize      = 'on';
    options.WindowStyle = 'normal';
    prompt = {['Set molar volume of  ',cell2mat(name_min),', V_m, as  [kJ/kBar]'],'Burgers vector, b [mm]','shear modulus, G [MPa]'};
    inpdlgtitle = cell2mat(name_min);
    defans = {'0.00','1e-7','1e+4'};
    par(1,:) = str2double(inputdlg(prompt,inpdlgtitle,1,defans,options));                               
    prompt1 = {'{A}_{dc} dislocation creep constant [MPa^-n mm^-m 1/s]',...
        '{Q}_l activation energy of lattice diffusion [kJ/mol]',...
        '{Q}_{gb} activation energy of grain-boundary diffusion [kJ/mol]'};
    inpdlgtitle1 = cell2mat(name_min);
    defans1 = {'1e-10','150','75'};
    par(2,:) = str2double(inputdlg(prompt1,inpdlgtitle1,1,defans1,options));   
    save(['DB_mineral_parameters\',char(name_min)],'par');
    % update mineral database as new mineral has been created
    listm = dir('DB_mineral_parameters\*.mat');
    [r,~]           = size(listm);
    listDB          = cell(r,1);
    for IDmineral = 1:r
        namewithext = listm(IDmineral).name;
        namewithoutext = strsplit(namewithext,'.');
        currname = namewithoutext(1);
        listDB(IDmineral) = currname;
    end
end

function edit_phase(~,~)
    % edit existent phase
    % the function allow to select an existent phase from DB and edit
    % one or more of its parameters
    %------------------------------------------------------------------
    [selmin,~] = listdlg('ListString',listDB,'Name','Select phase to be modified: ','SelectionMode','single');
    if isempty(selmin)==0
        options.Interpreter = 'tex';
        matmin     = load(['DB_mineral_parameters\',listDB{selmin},'.mat']);
        valedit = matmin.par;                   % extract flow law parameters of selected mineral
        promptmin = {'V, molar volume [kJ/kBar]','b, burgers vector [mm]','Shear modulus [MPa]'};
        titleditmin = ['Extensive properties of: ',selmin];
        defanseditmin = cell(1,3);
        for i = 1:3
            defanseditmin(i) = {num2str(valedit(1,i))};
        end
        par(1,:) = str2double(inputdlg(promptmin,titleditmin,1,defanseditmin,options));
        promptmin = {'DC pre-exponential constant [MPa^-n mm^-m 1/s]',...
            'Ql activation energy of lattice diffusion [kJ/mol]',...
            'Qgb activation energy of grain-boundary diffusion [kJ/mol]'};
        titleditmin = ['Creep parameters of: ',selmin];
        defanseditmin = cell(1,3);
        for i = 1:3
            defanseditmin(i) = {num2str(valedit(2,i))};
        end
        par(2,:) = str2double(inputdlg(promptmin,titleditmin,1,defanseditmin));
        save(['DB_mineral_parameters\',listDB{selmin}],'par');
    end
end

function del_phase(~,~)
    % the function allow deleting one of the existing phases from DB (permanent action)
    %-----------------------------------------------------------------------------------
    delcheck = questdlg('Are you sure you want to delete a phase from DB?','Delete phase','yes','no','no');
    switch delcheck
        case 'yes'
            [phasedeleted,~] = listdlg('ListString',listDB,'Name','Select phase being deleted: ','SelectionMode','single');
            delete(['DB_mineral_parameters\',listDB{phasedeleted},'.mat'])
        case 'no'    
    end                
end
    
function new_comp(~,~)
    % the function allow to create a new composite material composed of
    % n>=2 phases with specific proportions. The phases must be already included within the DB
    %------------------------------------------------------------------------------------------
    Nphase = inputdlg('Specify the number of phases (min = 1): ','Set Composite ',[1 50],{'1'});
    Nphase = str2double(Nphase);         % convert to double
    bulk = 0;                            % initialize total volume of the composite (must be 1.0 at the end of initialization)
    count = 0;
    phase_name = cell(1,Nphase);         % name of phases of composite material
    phase_vol = zeros(1,Nphase);
    phase_d = zeros(1,Nphase);
    while count<Nphase
        count = count+1;
        [comp_x,~] = listdlg('ListString',listDB,'Name','Select phase','SelectionMode','single');
        phase_name(count) = listDB(comp_x);
        pstring = {['Specify vol.% [min 0.0, max 1.0] of phase ',cell2mat(listDB(comp_x)),':'],'Set average grain-size [um]: '};
        val = str2double(inputdlg(pstring,cell2mat(listDB(comp_x)),1,{'1.0','100.0'}));
        if count<Nphase 
            if bulk<1.0
                bulk = bulk+val(1);
            else
                errordlg('The cumulative volume fraction of phases should not exceed 1.0!')
            end
        else
            if bulk+val(1)>1.0 || bulk+val(1)<1.0
                val(1) = 1-bulk;
                lname = length(cell2mat(listDB(comp_x)));
                rest1 = 19-lname;                         % spaces to be added for consistency
                newval = length(num2str(val(1)));
                rest2 = 64-newval;
                spaceadd = ' ';
                spaceadd1 = ' ';
                for i = 1:rest1
                    spaceadd(i) = ' ';
                end
                for i = 1:rest2
                    spaceadd1(i) = ' ';
                end
                errordlg(['The volume of the composite should not exceeds 1.0!                   ';...
                    'the volume fraction of phase ',cell2mat(listDB(comp_x)),' will be set to       ',spaceadd;...
                    'phi = ',num2str(val(1)),spaceadd1])    
            end
        end
        phase_vol(count) = val(1);
        phase_d(count) = val(2);
    end
    new_comp = inputdlg('Specify name of new composite material: ','New Composite',1,{'composite_X'});
    save(['DB_composite\',cell2mat(new_comp)],'phase_name','phase_vol','phase_d');
end

function del_comp(~,~)
    % the function allow to delete one of the existent composite
    % material from DB (permanent action!)
    %------------------------------------------------------------------
    filetodel = uigetfile('DB_composite');
    delete(['DB_composite\',filetodel]);
end

%% (A) RHEOLOGY----------------------------------------------------------------
%  this menu allow selection of the following settings:
%  1. elastic bound (Voigt/Reuss) -> functional to composite models only                                  
%  2. model type (single phase/composite)                               
%  3. projection space (2D/3D)                                  
%  4. authomatically visualize rheological profiles and DMM
%--------------------------------------------------------------------------
function select_phase(~,~)
    % select mineral, single-phase model
    % this switch off composite model and composite rheologies 
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(6:7),'Checked','off');
        set(m1h(6:7),'Enable','off');
        set(m1h(3),'Enable','on');      % enable composite models
        % reset variables
        model                = [];             % composite (1), mono-phase model (0)
        phase_name           = {};             % name of phases in composite (this is = 1 for single phase materials!)
        phase_ID             = [];             % ID of phase currently selected
        phase_pp             = [];             % ID of phase currently selected for piezometric analysis
        % change Data Zone properties
        set(modtxt1,'ForeGroundColor',inact_txt);
        set(mdtxt1_sel,'String',' ');
        set(modname,'ForeGroundColor',inact_txt);
        set(modname_sel,'String',' ');
        set(pmphase,'ForeGroundColor',inact_txt);
        set(pmphase_sel,'String',' ');
    else
        set(gcbo, 'Checked', 'on');
        [phase_ID,~] = listdlg('ListString',listDB,'Name','Select phase: ','SelectionMode','single');
        phase_name = listDB(phase_ID);
        
        % show selection message
        msgbox({'Single phase model selected: ';'Set mechanical constraint';'to continue'},'Model setup');
        
        % change menu properties
        set(m1h(3),'Enable','off');      % disabilitate composite models
        set(m1h(6:7),'Enable','on');     % mechanical constraints 5=V, 6=R
        set(m1h(22),'Enable','on');      % reset rheology
       
        % change datazone properties
        set(modtxt1,'ForeGroundColor',active_txt);
        set(mdtxt1_sel,'String','single mineral','ForeGroundColor',active_txt);
        set(modname,'ForeGroundColor',active_txt);
        set(modname_sel,'String',cell2mat(listDB(phase_ID)),'ForeGroundColor',active_txt);
        set(pmphase,'ForeGroundColor',active_txt);
        set(pmphase_sel,'String',phase_name);
        model = 0;
        phase_pp = phase_name;
    end
end

function load_comp(~,~)
    % load composite from DB.
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(6:7),'Enable','off');   % mechanical constraints 5=V, 6=R
        set(m1h(1),'Enable','on');      % enable single phase models
        % reset variables
        model                = [];      % composite (1), mono-phase model (0)
        cgs                  = [];      % median grain-size of phases in composite
        comp                 = [];      % name of composite
        phase_name           = {};      % name of phases in composite (this is = 1 for single phase materials!)
        cvol                 = [];      % phase proportions in composite
        % reset datazone
        set(modtxt1,'ForeGroundColor',inact_txt);
        set(mdtxt1_sel,'String','  ');
        set(modname,'ForeGroundColor',inact_txt);
        set(modname_sel,'String',' ');
        set(pmphase,'ForeGroundColor',inact_txt);
        set(pmphase_sel,'String',' ');
    else
        set(gcbo, 'Checked', 'on');
        % load composite material
        listc  = dir('DB_composite\*.mat');
        [r1,~] = size(listc);
        listCO = cell(r1,1);
        for compID1 = 1:r1
            cnamewithext    = listc(compID1).name;
            cnamewithoutext = strsplit(cnamewithext,'.');
            ccurrname       = cnamewithoutext(1);
            listCO(compID1)  = ccurrname;
        end
        [compname,~] = listdlg('ListString',listCO,'Name','Select composite: ','SelectionMode','single');
        comp = load(['DB_composite\',cell2mat(listCO(compname)),'.mat']);
        phase_name = comp.phase_name;
        
        % get ID of phases in current composite model
        phase_ID = zeros(1,numel(phase_name));
        ind = 1:numel(listDB);
        for i = 1:numel(phase_name)
            phase_ID(i) = ind(strcmp(phase_name(i),listDB));
        end 
        
        % get volume proportion, molar proportion and mean grainsize of all
        % phases in composite materials
        molp = molar_proportion(comp);
        cvol = comp.phase_vol;
        cgs = (comp.phase_d).*1e-3;              % grain-size in mm
        model = 1;
        
        % show selection message
        msgbox({'Composite model: ';'select mechanical constraint ';'to continue'},'Model setup');
        
        % change menu properties
        set(m1h(6:7),'Enable','on');     % mechanical constraints 5=V, 6=R
        set(m1h(1),'Enable','off');      % enable single phase models
        set(m1h(22),'Enable','on');      % reset rheology
        
        % change Data Zone properties
        set(modtxt1,'ForeGroundColor',active_txt);
        set(mdtxt1_sel,'String','composite','ForeGroundColor',active_txt);
        set(modname,'ForeGroundColor',active_txt);
        set(modname_sel,'String',cell2mat(listCO(compname)),'ForeGroundColor',active_txt);
        set(pmphase,'ForeGroundColor',active_txt);
        set(pmphase_sel,'String',phase_name);
    end
end

function voigt(~,~)
    % the function allow to select iso-strainrate condition, i.e., the
    % model assumes strain is uniform everywhere within the material.
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(9:13),'Enable','off');     % disable mixture rules 
        set(m1h(15:21),'Enable','off');    % disable projection spaces
        set(m1h(6:7),'Enable','on');       % enable selecting mechanical constraint
        % reset variables
        VR = [];
        mech = [];
        flowp = [];
        % reset Data Zone 
        set(modbehavior,'ForeGroundColor',inact_txt);
        set(modbehavior_sel,'String',' ');
    else
        set(gcbo, 'Checked', 'on');
        set(m1h(7),'Enable','off');
        VR = 1;
        % change Data Zone properties
        set(modbehavior,'ForeGroundColor',active_txt);
        set(modbehavior_sel,'String','iso-strain (Voigt)');
        if model==1
            % composite models
            msgbox({'Mechanical constraint selected: ';'set mixture rule to continue, or ';...
                ['plot the results in ',char(963),'-XY projection space']},'Model setup');
            % change menu properties
            set(m1h(9:13),'Enable','on');         % mixture rule 9=geometric,10=aritmetic,11=harmonic,12=huet,13=hobbs
            set(m1h(19),'Enable','on');           % Sii/xy projection space does not requires selecting a mixture rule
            set(m1h(22),'Enable','on');           % reset rheology
        else
            % single phase models
            msgbox({'Mechanical constraint selected: ';'set projection space to visualize results'},'Model setup');
            % change menu properties
            set(m1h([15,17,20,21]),'Enable','on'); % voigt's spaces + dmm
            set(m1h(22),'Enable','on');            % reset rheology
        end
    end
end

function reuss(~,~)
    % select isostress condition, i.e.,  the
    % model assumes stress is uniform everywhere within the material.
    % The function requires composite model selection (model = 1)
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu propertiesì
        set(m1h(9:13),'Enable','off');     % mixture rule 7=median,8=aritmetic,9=harmonic,10=huet,11=hobbs
        set(m1h(14:20),'Enable','off');    % disable projection spaces 
        set(m1h(6:7),'Enable','on');       % enable selecting mechanical constraint
        % reset variables
        VR = [];
        mech = [];
        flowp = [];
        % reset datazone
        set(modbehavior,'ForeGroundColor',inact_txt);
        set(modbehavior_sel,'String',' ');
    else
        set(gcbo, 'Checked', 'on');
        set(m1h(6),'Enable','off');
        VR = 0;
        % change Data Zone properties
        set(modbehavior,'ForeGroundColor',active_txt);
        set(modbehavior_sel,'String','iso-stress (Reuss)');
        if model==1
            % composite models
            msgbox({'Mechanical constraint selected: ';'set mixture rule to continue, or ';...
                ['plot the results in ',char(941),'-XY projection space']},'Model setup');
            % change menu properties
            set(m1h(9:13),'Enable','on');          % mixture rule 8=geometric,9=aritmetic,10=harmonic,11=huet,12=hobbs
            set(m1h(20),'Enable','on');            % Eii/xy projection space does not requires selecting a mixture rule
            set(m1h(22),'Enable','on');            % reset rheology
        else
            % single phase models
            msgbox({'Mechanical constraint selected: ';'set projection space to visualize results'},'Model setup');
            % change menu properties
            set(m1h([16,18,19,21]),'Enable','on'); % voigt's spaces + dmm
            set(m1h(22),'Enable','on');            % reset rheology
        end
    end  
end

function geomean(~,~) 
    % the  function computes composite rheology geometrically-averaging
    % the physical properties of the constituents. 
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(15:21),'Enable','off'); 
        % change Data Zone properties
        set(mixname,'ForeGroundColor',inact_txt);
        set(mixname_sel,'string',' ');
        mixr = '';
    else
        set(gcbo, 'Checked', 'on');
        % show selection message
        msgbox({'Mixture rule: Median ';'select one of the available projection space '; 'to continue'},'Model setup');
        mixr = 'median';
        % set flow law parameter formulations:
        % change menu properties
        set(m1h(9:13),'Checked','off');  % cancel previous selection if any
        set(m1h(22),'Enable','on');      % reset rheology
        if VR==0 
            % Reuss's projection spaces 15=Eii/T,17=Eii/d,19=Eii/xy
            set(m1h([16,18,20]),'Enable','on'); 
        else
           % Voigts's projection spaces 14=Sii/T,16=Sii/d,18=Sii/xy
           set(m1h([15,17,19]),'Enable','on');
        end
        % change Data Zone properties
        set(mixname,'ForeGroundColor',active_txt);
        set(mixname_sel,'string','Geometric mean');
    end
end

function arithmetic(~,~) 
     % (2) mixture model
     % select artihmetic mean (Voigt) mixture model 
     % it requires Voigt mechanical constraint
     %---------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(15:21),'Enable','off'); % all projection spaces 
        % change Data Zone properties
        set(mixname,'ForeGroundColor',inact_txt);
        set(mixname_sel,'String',' ');
        % reset variables
        mixr = ''; 
    else
        set(gcbo, 'Checked', 'on');
        % show selection message
        msgbox({'Mixture rule: arithmetic mean (Voigt) ';'select one of the available projection space '; 'to continue'},'Model setup');
        mixr = 'mean';
        % change menu properties
        set(m1h(22),'Enable','on');            % reset rheology
        set(m1h(9:13),'Checked','off');        % cancel previous selection if any
        if VR==0 
            % Reuss's projection spaces 15=Eii/T,17=Eii/d,19=Eii/xy
            set(m1h([16,18,20]),'Enable','on'); 
        else
           % Voigts's projection spaces 14=Sii/T,16=Sii/d,18=Sii/xy
           set(m1h([15,17,19]),'Enable','on');
        end
        % change Data Zone properties
        set(mixname,'ForeGroundColor',active_txt);
        set(mixname_sel,'String','Arithmetic');
    end
end

function harmonic(~,~) 
     % (2) mixture model
     % select harmonic mean (Reuss) mixture model 
     % it requires Reuss mechanical constraint
     %---------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(15:21),'Enable','off'); % all projection spaces 
        % change Data Zone properties
        set(mixname,'ForeGroundColor',inact_txt);
        set(mixname_sel,'String',' ');
        % reset variables
        mixr = '';   
    else
        set(gcbo, 'Checked', 'on');
        % show selection message
        msgbox({'Mixture rule: harmonic mean (Reuss) ';'select one of the available projection space '; 'to continue'},'Model setup');
        mixr = 'harmonic';
        % change menu properties
        set(m1h(22),'Enable','on');            % reset rheology
        set(m1h(9:13),'Checked','off');        % cancel previous selection if any
        if VR==0 
            % Reuss's projection spaces 15=Eii/T,17=Eii/d,19=Eii/xy
            set(m1h([16,18,20]),'Enable','on'); 
        else
           % Voigts's projection spaces 14=Sii/T,16=Sii/d,18=Sii/xy
           set(m1h([15,17,19]),'Enable','on');
        end
        % change Data Zone properties
        set(mixname,'ForeGroundColor',active_txt);
        set(mixname_sel,'String','Harmonic');  
    end
end

function huet(~,~)
    % (3) mixture model
    % select minimized power geometric mean model. Could be either Reuss-
    % or Voigt-based depending on the selected mechanical constraint
    %----------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(15:21),'Enable','off');
        % change datazone
        set(mixname,'ForeGroundColor',inact_txt);
        set(mixname_sel,'String',' ');
        % reset variables
        mixr = '';
    else
        set(gcbo, 'Checked', 'on');
        % show selection message
        msgbox({'Mixture rule: Minimized Power model ';'select one of the available projection space '; 'to continue'},'Model setup');
        % change menu properties
        set(m1h(22),'Enable','on');      % reset rheology
        set(m1h(9:13),'Checked','off');  % cancel previous selection if any
        mixr = 'huet';
        if VR==0
            set( m1h([16,20]),'Enable','on');  
        else
            set(m1h([15,19]),'Enable','on');
        end
        % change Data Zone properties
        set(mixname,'ForeGroundColor',active_txt);
        set(mixname_sel,'String','Huet et al. (2014)');
    end
end

function hobbs(~,~)
    % (4) mixture model
    % select thermodynamic mixing model 
    %----------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        % change menu properties
        set(m1h(15:21),'Enable','off');
        % change datazone
        set(mixname,'ForeGroundColor',inact_txt);
        set(mixname_sel,'String',' ');
        % reset variables
        mixr = '';
    else
        set(gcbo, 'Checked', 'on');
        % show selection message
        msgbox({'Mixture rule: Thermodynamic mixing model ';'select one of the available projection space '; 'to continue'},'Model setup');
        % change menu properties
        set(m1h(22),'Enable','on');      % reset rheology
        set(m1h(9:13),'Checked','off');  % cancel previous selection if any
        if VR==0
            set( m1h([16,20]),'Enable','on');  
        else
            set(m1h([15,19]),'Enable','on');
        end
        % change Data Zone properties
        set(mixname,'ForeGroundColor',active_txt);
        set(mixname_sel,'String','Hobbs et al. (2019)');
        % set mixture rule ID
        mixr = 'hobbs';
    end 
end

function Sii_T(~,~) 
% Set stress/T as projection space - requires specifying a reference strain
% rate.
% 1. For single-phase models, this requires specifying
% grain-size as the fixed variable. 
% 2. composite have fixed grain-size thus specifying grain-size 
% is not required
% -----------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(pmtxt1,'ForeGroundColor',inact_txt);
        set(pmtxt1_sel,'String',' ');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 1;
        set(m1h(15:21),'Checked','off');  % cancel previous selection 
        set(gcbo, 'Checked', 'on');
        % change menu properties
        set(m1h(22),'Enable','on');       % reset rheology
        set(m2h(7),'Enable','on');        % reset graph on
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Stress/temperature');  
        if model==0
            % single phase model require specifying both strainrate and
            % grainsize
            options.Interpreter = 'tex';
            defans = {'1e-15','0.01'};
            prompt = {'Set a reference strain-rate of deformation [s^{-1}]: ','Set grain size [mm]: '};
            fvar = str2double(inputdlg(prompt,'Model setup',1,defans,options));
            Eii = fvar(1);
            cgs = fvar(2);
            [mech,Siimin] = Sii_T_voigt(model,phase_name,Taxis,nXY,cgs,Eii);
            % change grainsize slider properties
            set(dsld3,'ForeGroundColor',active_txt);
            set(dsld2,'ForeGroundColor',active_txt);
            set(dsld1,'Enable','on','BackGroundColor',[1 1 1]);
            set(dsld1,'Value',cgs);
            % change dynamic legend properties
            set(legend_intvar2,'String',['grain size = ',num2str(cgs),' [mm]']);      
            set(legend_intvar2,'Visible','on');
        else
            % composite models require only strainrate
            options.Interpreter = 'tex';
            defans = {'1e-15'};
            prompt = {'Set a reference strain-rate of deformation [s^{-1}]: '};
            Eii = str2double(inputdlg(prompt,'Model setup',1,defans,options));
            [mech,Siimin] = Sii_T_voigt(model,phase_name,Taxis,nXY,cgs,Eii,molp,cvol,mixr);
        end 
        % change slider properties
        set(strainsld3,'ForeGroundColor',active_txt);
        set(strainsld2,'ForeGroundColor',active_txt);
        set(strainsld1,'Enable','on','BackGroundColor',[1 1 1]);
        set(strainsld1,'Value',log10(Eii));
        % change dynamic legend properties
        set(legend_intvar1,'String',['Log(',char(941),') = ',num2str(log10(Eii)),' [1/s]']);      
        set(legend_intvar1,'Visible','on');
        % show selection message
        msg = msgbox(['Visualizing results in ',char(963),'/T projection space!'],'results');
        uiwait(msg);
        cla(h2,'reset');
        set(h2,'Visible','off')
        Xl = [char(963),' [MPa]'];
        Yl = 'T [°C]';  
        lab1 = ['mean grainsize = ',num2str(mean(cgs),'%.3f'),' mm'];
        lab2 = ['Log_{10}(',char(941),')',num2str(log10(Eii),'%.1f'),' s^{-1}'];
        [h2] = plot1d(lab1,lab2,Siimin,Taxis-273.15,Siiaxis,mech,Xl,Yl,model);
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end
end

function Eii_T(~,~) 
% 1. For single-phase models, this requires specifying
% grain-size as the fixed variable. 
% 2. composite have fixed grain-size thus specifying grain-size 
% is not required
%------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(pmtxt1,'ForeGroundColor',inact_txt);
        set(pmtxt1_sel,'String',' ');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 2;
        set(m1h(15:21),'Checked','off');  % cancel previous selection
        set(gcbo, 'Checked', 'on');
        % change menu properties
        set(m2h(7),'Enable','on');        % reset graph on
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Strain-rate/temperature');
        if model==0
            % single phase model require specifying both stress and grainsize
            options.Interpreter = 'tex';
            defans = {'100','0.01'};
            prompt = {'Set a reference stress [MPa]: ','Set grain size [mm]: '};
            fvar = str2double(inputdlg(prompt,'Model setup',1,defans,options));
            Sii = fvar(1);
            cgs = fvar(2);  
            [mech,Eiimax] = Eii_T_reuss(model,phase_name,Taxis,nXY,cgs,Sii);
            % change grainsize slider properties
            set(dsld3,'ForeGroundColor',active_txt);
            set(dsld2,'ForeGroundColor',active_txt);
            set(dsld1,'Enable','on','BackGroundColor',[1 1 1]);
            set(dsld1,'Value',cgs);
            % change dynamic legend properties
            set(legend_intvar2,'String',['grain size = ',num2str(cgs),' [mm]']);      
            set(legend_intvar2,'Visible','on');
        else
            % composite model require only mean stress
            options.Interpreter = 'tex';
            defans = {'100'};
            prompt = {'Set reference stress [MPa]: '};
            Sii = str2double(inputdlg(prompt,'Model setup',1,defans,options));
            [mech,Eiimax] = Eii_T_reuss(model,phase_name,Taxis,nXY,cgs,Sii,molp,cvol,mixr); 
        end
        % change stress slider properties
        set(stresssld3,'ForeGroundColor',active_txt);
        set(stresssld2,'ForeGroundColor',active_txt);
        set(stresssld1,'Enable','on','BackGroundColor',[1 1 1]);
        set(stresssld1,'Value',Sii);
        % change dynamic legend properties
        set(legend_intvar1,'String',[char(963),' = ',num2str(Sii),' [MPa]']);      
        set(legend_intvar1,'Visible','on');
        % show selection message
        msg = msgbox(['Visualizing results in ',char(941),'/T projection space!'],'results');
        uiwait(msg);   
        cla(h2,'reset');
        set(h2,'Visible','off')
        Xl = [char(941),' [s^{-1}]'];
        Yl = 'T [°C]';
        lab1 = ['grainsize = ',num2str(mean(cgs),'%.3f'),' mm'];
        lab2 = ['Mean stress ',char(941),')',num2str(Sii,'%.1f'),' MPa'];
        [h2] = plot1d(lab1,lab2,Eiimax,Taxis-273.15,Eiiaxis,mech,Xl,Yl,model,'logX');        
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end
end

function Sii_d(~,~)
% Set stress-grainsize as projection space.
% For single-phase models only, this requires specifying T and Eii as the 
% fixed variables. NOT SUPPORTED FOR COMPOSITE MODELS!
%------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(pmtxt1,'ForeGroundColor',inact_txt);
        set(pmtxt1_sel,'String',' ');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 3;
        set(m1h(15:21),'Checked','off');  % cancel previous selection
        set(gcbo, 'Checked', 'on');
        % change menu properties
        set(m2h(7),'Enable','on');        % reset graph on
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Stress/grain size');
        % single phase models requires specifying temperature and
        % strain-rate
        options.Interpreter = 'tex';
        defans = {'1e-15','600'};
        prompt = {'Set a reference strain-rate of deformation [s^{-1}]: ','Set temperature [°C]: '};
        fvar = str2double(inputdlg(prompt,'Model setup',1,defans,options));
        Eii = fvar(1);
        T = fvar(2);
        % change slider properties
        set(tempsld3,'ForeGroundColor',active_txt);
        set(tempsld2,'ForeGroundColor',active_txt);
        set(tempsld1,'Enable','on','BackGroundColor',[1 1 1])
        set(tempsld1,'Value',T);
        set(strainsld3,'ForeGroundColor',active_txt);
        set(strainsld2,'ForeGroundColor',active_txt);
        set(strainsld1,'Enable','on','BackGroundColor',[1 1 1]);
        set(strainsld1,'Value',log10(Eii));
        [mech,Siimin] = Sii_d_voigt(phase_name,Daxis,nXY,T,Eii);
        % change dynamic legend properties
        set(legend_intvar1,'String',['Log(',char(941),') = ',num2str(log10(Eii)),' [1/s]']);      
        set(legend_intvar1,'Visible','on');
        set(legend_intvar2,'String',['Temperature = ',num2str(T,'%.1f'),' [°C]']);      
        set(legend_intvar2,'Visible','on');
        % show selection message
        msg = msgbox(['Visualizing results in ',char(963),'/d projection space!'],'results');
        uiwait(msg);
        cla(h2,'reset');
        set(h2,'Visible','off')
        Xl = [char(963),' [MPa]'];
        Yl = 'Log_{10}d [mm]';
        lab1 = ['Temperature = ',num2str(T,'%.1f'),'°C'];
        lab2 = ['Log_{10}(',char(941),')',num2str(log10(Eii),'%.1f'),'s^{-1}'];
        [h2] = plot1d(lab1,lab2,Siimin,Daxis,Siiaxis,mech,Xl,Yl,model,'logY');
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end
end

function Eii_d(~,~)
% Set strainrate grainsize as the projection space.
% For single-phase models, this requires specifying T and Sii as the 
% fixed variables. NOT SUPPORTED FOR COMPOSITE MODELS!
%------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(pmtxt1,'ForeGroundColor',inact_txt);
        set(pmtxt1_sel,'String',' ');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 4;
        set(m1h(15:21),'Checked','off');  % cancel previous selection
        set(gcbo, 'Checked', 'on');
        % change menu properties
        set(m2h(7),'Enable','on');        % reset graph on
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Strain-rate/grain size');
        % single phase models requires specifying temperature and stress
        options.Interpreter = 'tex';
        defans = {'50','600'};
        prompt = {'Set reference stress [MPa]: ','Set temperature [°C]: '};
        fvar = str2double(inputdlg(prompt,'Model setup',1,defans,options));
        Sii = fvar(1);
        T = fvar(2);
        % change slider properties
        set(tempsld3,'ForeGroundColor',active_txt);
        set(tempsld2,'ForeGroundColor',active_txt);
        set(tempsld1,'Enable','on','BackGroundColor',[1 1 1])
        set(tempsld1,'Value',T);
        set(stresssld3,'ForeGroundColor',active_txt);
        set(stresssld2,'ForeGroundColor',active_txt);
        set(stresssld1,'Enable','on','BackGroundColor',[1 1 1]);
        set(stresssld1,'Value',Sii);
        % change dynamic legend properties
        set(legend_intvar1,'String',[char(963),' = ',num2str(Sii),' [MPa]']);      
        set(legend_intvar1,'Visible','on');
        set(legend_intvar2,'String',['Temperature = ',num2str(T,'%.1f'),' [°C]']);      
        set(legend_intvar2,'Visible','on');
        [mech,Eiimin] = Eii_d_reuss(phase_name,Daxis,nXY,T,Sii);
        % show selection message
        msg = msgbox(['Visualizing results in ',char(941),'/d projection space!'],'results');
        uiwait(msg);
        cla(h2,'reset');
        set(h2,'Visible','off')
        Xl = ['Log_{10}(',char(941),') [s^{-1}]'];
        Yl = 'Log_{10}(d) [mm]';
        lab1 = ['Temperature = ',num2str(T,'%.1f'),'°C'];
        lab2 = [char(963),' = ',num2str(Sii,'%.1f'),'MPa'];
        [h2] = plot1d(lab1,lab2,Eiimin,Daxis,Eiiaxis,mech,Xl,Yl,model,'logY','logX');
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end
end

function Sii_xy(~,~)
% Colormap of stress pattern on 2D phase map. Stress values are
% computed for each grain based on their composition and grainsize,
% assuming that the stable deformation mechanism is the one that
% locally minimizes stress (i.e., Voigt mechanical constraint
% is required)
% this projection requires: 1) specifying a temperature of deformation
% and mean strain rate, 2) uploading a phase-map showing the boundaries between
% different grains
%------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 5;
        SiiEiiMap = 1;
        % change menu properties
        set(gcbo, 'Checked', 'on');
        set(m1h(15:21),'Checked','off');  % cancel previous selection 
        set(m1h(22),'Enable','on');       % reset rheology
        set(m2h(7),'Enable','on');        % reset graph on
        set(radio1,'Enable','on','ForeGroundColor',active_txt);
        set(radio2,'Enable','on','ForeGroundColor',active_txt);
        set(radio3,'Enable','on','ForeGroundColor',active_txt);
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Stress field');
        % show selection message
        msg = msgbox(['Visualizing results in ',char(963),'/XY projection space!'],'results');
        uiwait(msg);     
        % upload image to construct phasemap and analyze grain diameters
        [mapgrain_ID,grain_comp,mapgrain_D,mapsize,wmap,mapbound] = grainsize(model,listDB,currentFolder,phase_name,phase_ID);
        options.Interpreter = 'tex';
        options.Resize      = 'on';
        options.WindowStyle = 'normal';
        prompt = {'Temperature [°C]','strain-rate [1/s]'};
        inputdlgtitle = 'Set T and strain-rate';
        defans = {'600','1e-14'};
        data_intensive = str2double(inputdlg(prompt,inputdlgtitle,1,defans,options));
        T = data_intensive(1);
        Eii = data_intensive(2);
        % change slider properties
        set(tempsld3,'ForeGroundColor',active_txt);
        set(tempsld2,'ForeGroundColor',active_txt);
        set(tempsld1,'Enable','on','BackGroundColor',[1 1 1])
        set(tempsld1,'Value',T);
        set(strainsld3,'ForeGroundColor',active_txt);
        set(strainsld2,'ForeGroundColor',active_txt);
        set(strainsld1,'Enable','on','BackGroundColor',[1 1 1]);
        set(strainsld1,'Value',log10(Eii));
        % change dynamic legend properties
        set(legend_intvar1,'String',['Log(',char(941),') = ',num2str(log10(Eii)),' [1/s]']);      
        set(legend_intvar1,'Visible','on');
        set(legend_intvar2,'String',['Temperature = ',num2str(T,'%.1f'),' [°C]']);      
        set(legend_intvar2,'Visible','on');
        
        % compute stress and stable deformation mechanism
        [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
        
        % compute deviations from mole-averaged median stress
        [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
        
        % update mean stress and peak stress values, deviations from mean stress and legend values
        set(legend_medianVal,'String',['mean stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
        set(legend_medianVal,'Visible','on');
        set(legend_peakVal,'String',['peak stress = ',num2str(max(max(Sii_xymap)),'%.1f'),' [MPa]']);
        set(legend_peakVal,'Visible','on');
        
        % plot the results (default shows the stress distribution resulting from
        Sii_map  = zeros(mapsize(1),mapsize(2));
        for i = 1:numel(mapgrain_D)
            Sii_map(mapgrain_ID{i}) = Sii_xymap(i);
        end
        
        % hide grain-boundaries
        for i = 1:mapsize(1)*mapsize(2)
            if mapbound(i)==0
                Sii_map(i) = nan;
            end
        end
        % plot phasemaps and deformation mechanism map
        SiiX = linspace(0,wmap,mapsize(2));
        SiiY = linspace(0,(wmap/mapsize(2))*mapsize(1),mapsize(1));
        [gridSiiX,gridSiiY] = meshgrid(SiiX,SiiY);
        % clear previous plots
        cla(h2,'reset');
        set(h2,'Visible','off')
        h2 = plotSii_xy(gridSiiX,gridSiiY,Sii_map,mapsize,siimap);
        set(legend_Siixy,'Visible','on');
        % enable data dynamic display 
        dcm.UpdateFcn = @cursorXYmap_voigt;
    end
end

function Eii_xy(~,~)
    % Colormap of strainrate pattern on 2D phase map. Strainrate values are
    % computed for each grain based on their composition and grainsize,
    % assuming that the stable deformation mechanism is the one that
    % locally maximixes strainrate (i.e., Reuss mechanical constraint
    % is required)
    % this projection requires: 1) specifying a temperature of deformation
    % and mean stress, 2) uploading a phase-map showing the boundaries between
    % different grains
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 6;
        SiiEiiMap = 2;
        % change menu properties
        set(gcbo, 'Checked', 'on');
        set(m1h(15:21),'Checked','off');  % cancel previous selection 
        set(m1h(22),'Enable','on');       % reset rheology
        set(m2h(7),'Enable','on');        % reset graph on
        set(radio1,'Enable','on','ForeGroundColor',active_txt);
        set(radio2,'Enable','on','ForeGroundColor',active_txt)
        set(radio3,'Enable','on','ForeGroundColor',active_txt);
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Strain rate field');
        % show selection message
        msg = msgbox(['Visualizing results in ',char(941),'/XY projection space!'],'results');
        uiwait(msg);
        % upload image to construct phasemap and analyze grain diameters
        [mapgrain_ID,grain_comp,mapgrain_D,mapsize,wmap,mapbound] = grainsize(model,listDB,currentFolder,phase_name,phase_ID);
        options.Interpreter = 'tex';
        options.Resize      = 'on';
        options.WindowStyle = 'normal';
        prompt = {'Temperature [°C]','stress [MPa]'};
        inputdlgtitle = 'Set T and stress';
        defans = {'600','100'};
        data_intensive = str2double(inputdlg(prompt,inputdlgtitle,1,defans,options));
        T = data_intensive(1);
        Sii = data_intensive(2);          
        % change slider properties
        set(tempsld3,'ForeGroundColor',active_txt);
        set(tempsld2,'ForeGroundColor',active_txt);
        set(tempsld1,'Enable','on','BackGroundColor',[1 1 1])
        set(tempsld1,'Value',T);
        set(stresssld3,'ForeGroundColor',active_txt);
        set(stresssld2,'ForeGroundColor',active_txt);
        set(stresssld1,'Enable','on','BackGroundColor',[1 1 1]);
        set(stresssld1,'Value',Sii);
        % change dynamic legend properties
        set(legend_intvar1,'String',[char(963),' = ',num2str(Sii,'%.1f'),' [MPa]']);      
        set(legend_intvar1,'Visible','on');
        set(legend_intvar2,'String',['Temperature = ',num2str(T,'%.1f'),' [°C]']);      
        set(legend_intvar2,'Visible','on');
        
        % compute log10 strainrate and stable deformation mechanism
        [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
        
        % compute deviations from mole-averaged median strain-rate
        [Eii_dev,Eii_median] = median_dev_Eii(Eii_xymap,mapgrain_D,grain_comp,listDB);
        
        % update mean strain rate value, deviations from mean strain rate and legend values
        set(legend_medianVal,'String',['median strainrate = 1e',num2str(Eii_median,'%.1f'),' [s-1]']);
        set(legend_medianVal,'Visible','on');
        set(legend_peakVal,'String',['max strainrate = 1e',num2str(max(max(Eii_xymap)),'%.1f'),' [s-1]']);
        set(legend_peakVal,'Visible','on');
        Eii_map  = zeros(mapsize(1),mapsize(2));
        for i = 1:numel(mapgrain_D)
            Eii_map(mapgrain_ID{i}) = Eii_xymap(i);
        end
        % hide grain-boundaries
        for i = 1:mapsize(1)*mapsize(2)
            if mapbound(i)==0
                Eii_map(i) = nan;
            end
        end
        % plot phasemaps and deformation mechanism map
        EiiX = linspace(0,wmap,mapsize(2));
        EiiY = linspace(0,(wmap/mapsize(2))*mapsize(1),mapsize(1));
        [gridEiiX,gridEiiY] = meshgrid(EiiX,EiiY);
        % plot the results
        cla(h2,'reset');
        set(h2,'Visible','off')
        [h2] = plotEii_xy(gridEiiX,gridEiiY,Eii_map,mapsize,siimap); 
        set(legend_Eiixy,'Visible','on');
        % enable data dynamic display 
        dcm.UpdateFcn = @cursorXYmap_reuss;
    end
end

function dmm(~,~)
    % plot deformation-mechanism maps: the field of various deformation mechanisms
    % (dislocation creep, diffusion(coble) creep and grain boundary
    % slidind) are plotted in (stress/strainrate)-grainsize space
    %----------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(legend_intvar1,'Visible','off');
        set(legend_intvar2,'Visible','off');
    else
        plottype = 7;
        % change menu properties
        set(gcbo, 'Checked', 'on');
        set(m1h(15:21),'Checked','off');  % cancel previous selection 
        set(m1h(22),'Enable','on');       % reset rheology
        set(m2h(7),'Enable','on');        % reset graph on
        % change datazone properties
        set(pmtxt1,'ForeGroundColor',active_txt);
        set(pmtxt1_sel,'String','Stress/grain size');
        if model==0
            % get rheological parameters of current phase
            flowp = set_flowp(phase_name);
        end
        % set temperature of deformation [K] and eval. stable deformation
        % mechanism
        T = str2double(inputdlg('Set temperature of deformation [°C]:'))+273.15;
        % change slider properties
        set(tempsld3,'ForeGroundColor',active_txt);
        set(tempsld2,'ForeGroundColor',active_txt);
        set(tempsld1,'Enable','on','BackGroundColor',[1 1 1])
        set(tempsld1,'Value',T);
        % change dynamic legend properties
        set(legend_intvar2,'String',['Temperature = ',num2str(T-273.15,'%.1f'),' [°C]']);      
        set(legend_intvar2,'Visible','on');
        if VR == 1
            % voigt condition, def.mech. are evaluated based on minimizing stress 
            % show selection message
            msg = msgbox(['Visualizing results in ',char(963),'/d projection space!'],'results');
            uiwait(msg);
            [mech] = dmm_voigt(flowp,nXY,T,Daxis,Eiiaxis);
            cla(h2,'reset');
            set(h2,'Visible','off')
            h2 = plot_dmm(Siiaxis,Daxis,mech,mechcmap,T,'voigt');
        else
            % reuss's condition, def.mech. are evaluated based on maximizing strainrate
            msg = msgbox(['Visualizing results in ',char(941),'/d projection space!'],'results');
            uiwait(msg);
            [mech] = dmm_reuss(flowp,nXY,T,Daxis,Siiaxis);
            cla(h2,'reset');
            set(h2,'Visible','off')
            h2 = plot_dmm(Siiaxis,Daxis,mech,mechcmap,T,'reuss');
        end    
    end
end

function plotselection(~,event)
    % callback controlling the behavior of the radiobuttons that allow
    % switching between deformation mechanism maps and stress/strain-rate
    % maps for 2D plots in the X,Y-V space
    %---------------------------------------------------------------------
    switch event.NewValue.String
        case 'stress/strain-rate map'
            set(legend_Siixy,'Visible','off');
            set(legend_Eiixy,'Visible','off');
            set(legend_deviation_Eii,'Visible','off');
            set(legend_deviation_Sii,'Visible','off');
            set(legend_dmm(1:9),'Visible','off')  
            cla(h2,'reset')
            set(h2,'Visible','off')
            if VR==1
                % voigt maps
                h2 = plotSii_xy(gridSiiX,gridSiiY,Sii_map,mapsize,siimap);
                plottype = 5;
                set(legend_Siixy,'Visible','on');
                % enable data dynamic display 
                dcm.UpdateFcn = @cursorXYmap_voigt;
            else 
                % reuss model
                h2 = plotEii_xy(gridEiiX,gridEiiY,Eii_map,mapsize,siimap);
                plottype = 6;
                set(legend_Eiixy,'Visible','on');
                % enable data dynamic display 
                dcm.UpdateFcn = @cursorXYmap_reuss;
            end
        case 'deformation mechanism map'
            set(legend_Siixy,'Visible','off');
            set(legend_Eiixy,'Visible','off');
            set(legend_deviation_Eii,'Visible','off');
            set(legend_deviation_Sii,'Visible','off');
            set(legend_dmm(1:9),'Visible','on')  
            plottype = 8;
            cla(h2,'reset')
            set(h2,'Visible','off')
            h2 = plotDMM_xy(wmap,mapsize,mapbound,mapgrain_D,mapgrain_ID,eff_mech_xy,mapsize,mechcmap); 
            % enable data dynamic display 
            dcm.UpdateFcn = @cursorXYmap_dmm;
        case 'standard deviation map'
            set(legend_Siixy,'Visible','off');
            set(legend_Eiixy,'Visible','off');
            set(legend_deviation_Eii,'Visible','off');
            set(legend_deviation_Sii,'Visible','off')
            set(legend_dmm(1:9),'Visible','off') 
            set(legend_medianVal,'Visible','on');
            plottype = 9;
            cla(h2,'reset')
            set(h2,'Visible','off')
            if SiiEiiMap == 1
                % switch is from Sii_xy maps 
                Sii_devmap  = zeros(mapsize(1),mapsize(2));
                for i = 1:numel(mapgrain_D)
                    Sii_devmap(mapgrain_ID{i}) = Sii_dev(i);
                end
                % hide grain-boundaries
                for i = 1:mapsize(1)*mapsize(2)
                    if mapbound(i)==0
                        Sii_devmap(i) = nan;
                    end
                end
                h2 = plot_deviationSii(gridSiiX,gridSiiY,Sii_devmap,mapsize,cmapdiff);
                set(legend_deviation_Sii,'Visible','on');
                set(legend_medianVal,'String',['mean stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
                % enable data dynamic display 
                dcm.UpdateFcn = @cursorXYmap_deviationSii;
            else
                % switch is from Eii_xy maps
                Eii_devmap  = zeros(mapsize(1),mapsize(2));
                for i = 1:numel(mapgrain_D)
                    Eii_devmap(mapgrain_ID{i}) = Eii_dev(i);
                end
                % hide grain-boundaries
                for i = 1:mapsize(1)*mapsize(2)
                    if mapbound(i)==0
                        Eii_devmap(i) = nan;
                    end
                end
                h2 = plot_deviationEii(gridEiiX,gridEiiY,Eii_devmap,mapsize,cmapdiff);
                set(legend_deviation_Eii,'Visible','on');
                set(legend_medianVal,'String',['median strainrate = 1e',num2str(Eii_median,'%.1f'),' [s-1]']);
                set(legend_peakVal,'String',['max strainrate = 1e',num2str(max(max(Eii_xymap)),'%.1f'),' [s-1]']);
                % enable data dynamic display 
                dcm.UpdateFcn = @cursorXYmap_deviationEii;
            end
    end
end

function slidEii(hObject,~)
    % the function returns the strain-rate value as set by the end user
    Eii = get(hObject,'Value');
    Eii = 10^Eii;
    % update dynamic legend properties
    set(legend_intvar1,'String',['Log(',char(941),') = ',num2str(log10(Eii)),' [1/s]']);
    set(m2h(1:5),'Checked','off');
    switch plottype
        case 1 % Sii/T projection space; variables Eii + d(single-phase models only)
            % update stress based on new strain-rate 
            if model==0
                [mech,Siimin] = Sii_T_voigt(model,phase_name,Taxis,nXY,cgs,Eii);
            else
                [mech,Siimin] = Sii_T_voigt(model,phase_name,Taxis,nXY,cgs,Eii,molp,cvol,mixr);
            end
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = [char(963),' [MPa]'];
            Yl = 'T [°C]';  
            lab1 = ['mean grainsize = ',num2str(mean(cgs),'%.3f'),'mm'];
            lab2 = ['Log_{10}(',char(941),')',num2str(log10(Eii),'%.1f'),'s^{-1}'];
            h2 = plot1d(lab1,lab2,Siimin,Taxis-273.15,Siiaxis,mech,Xl,Yl,model);      
        case 3 % Sii/d projection space; variables Eii,T
            [mech,Siimin] = Sii_d_voigt(phase_name,Daxis,nXY,T,Eii);
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = [char(963),' [MPa]'];
            Yl = 'Log_{10}d [mm]';
            lab1 = ['Temperature = ',num2str(T-273.15,'%.1f'),'°C'];
            lab2 = ['Log_{10}(',char(941),')',num2str(log10(Eii),'%.1f'),'s^{-1}'];
            h2 = plot1d(lab1,lab2,Siimin,Daxis,Siiaxis,mech,Xl,Yl,model,'logY');
        case 5 % Sii/xy map (Eii,T)
            % update stress and stable deformation mechanism based on new strain-rate
            [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
            Sii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Sii_map(mapgrain_ID{i}) = Sii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Sii_map(i) = nan;
                end
            end
            % clear previous plot and redraw the map
            cla(h2,'reset');
            colorbar off
            set(h2,'Visible','off')
            h2 = plotSii_xy(gridSiiX,gridSiiY,Sii_map,mapsize,siimap);
            % update mean stress value, deviations from mean stress and legend values
            [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
            set(legend_medianVal,'String',['mean stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
            set(legend_peakVal,'String',['peak stress = ',num2str(max(max(Sii_xymap)),'%.1f'),' [MPa]']);
        case 8 % XY deformation mechanism maps
            % update stress and stable deformation mechanism based on new strain-rate
            [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
            % update map of grain stresses
            Sii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Sii_map(mapgrain_ID{i}) = Sii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Sii_map(i) = nan;
                end
            end
            cla(h2,'reset')
            set(h2,'Visible','off')
            h2 = plotDMM_xy(wmap,mapsize,mapbound,mapgrain_D,mapgrain_ID,eff_mech_xy,mapsize,mechcmap);
            % update mean stress value, deviations from mean stress and legend values
            [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
            set(legend_medianVal,'String',['mean stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
        case 9 % stress deviation maps
            % compute deviations from mole-averaged median stress
            % update stress and stable deformation mechanism based on new strain-rate
            [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
            % update map of grain stresses
            Sii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Sii_map(mapgrain_ID{i}) = Sii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Sii_map(i) = nan;
                end
            end
            [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
            Sii_devmap  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Sii_devmap(mapgrain_ID{i}) = Sii_dev(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Sii_devmap(i) = nan;
                end
            end
            cla(h2,'reset');
            set(h2,'Visible','off')
            h2 = plot_deviationSii(gridSiiX,gridSiiY,Sii_devmap,mapsize,cmapdiff);
            % update median strain-rate value display
            set(legend_medianVal,'String',['median stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
            % enable data dynamic display 
            dcm.UpdateFcn = @cursorXYmap_deviationSii;
    end
end

function slidT(hObject,~)
    % the function returns the temperature value as set by the end user
    T = get(hObject,'Value');
    % update dynamic legend properties
    set(legend_intvar2,'String',['Temperature = ',num2str(T,'%.1f'),' [°C]']);  
    hold off
    set(m2h(1:5),'Checked','off');
    switch plottype
        case 3 % Sii/d projection space (Eii,T)
            [mech,Siimin] = Sii_d_voigt(phase_name,Daxis,nXY,T,Eii);
            % clear previous plot and redraw the map
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = [char(963),' [MPa]'];
            Yl = 'Log_{10}d [mm]';
            lab1 = ['Temperature = ',num2str(T-273.15,'%.1f'),'°C'];
            lab2 = ['Log_{10}(',char(941),')',num2str(log10(Eii),'%.1f'),'s^{-1}'];
            h2 = plot1d(lab1,lab2,Siimin,Daxis,Siiaxis,mech,Xl,Yl,model,'logY');
        case 4 % Eii/d projection space (Sii/T)
            [mech,Eiimin] = Eii_d_reuss(phase_name,Daxis,nXY,T,Sii);
            % clear previous plot and redraw the map
            Xl = ['Log_{10}(',char(941),') [s^{-1}]'];
            Yl = 'Log_{10}(d) [mm]';
            lab1 = ['Temperature = ',num2str(T-273.15,'%.1f'),'°C'];
            lab2 = [char(963),' = ',num2str(Sii,'%.1f'),'MPa'];
            h2 = plot1d(lab1,lab2,Eiimin,Daxis,Eiiaxis,mech,Xl,Yl,model,'logY','logX');
        case 5 % Sii/xy map (Eii,T)
            % update stress and stable deformation mechanism based on new T
            [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
            Sii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Sii_map(mapgrain_ID{i}) = Sii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Sii_map(i) = nan;
                end
            end
            % clear previous plot and redraw the map
            cla(h2,'reset');
            colorbar off
            set(h2,'Visible','off')
            h2 = plotSii_xy(gridSiiX,gridSiiY,Sii_map,mapsize,siimap);
            % update mean stress value, deviations from mean stress and legend values
            [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
            set(legend_medianVal,'String',['mean stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
            set(legend_peakVal,'String',['peak stress = ',num2str(max(max(Sii_xymap)),'%.1f'),' [MPa]']);
        case 6 % Eii/xy map (Sii,T)
            % update strainrate and stable deformation mechanism based on new T
            [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
            Eii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Eii_map(mapgrain_ID{i}) = Eii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Eii_map(i) = nan;
                end
            end
            % clear previous plot and redraw the map
            cla(h2,'reset');
            colorbar off
            set(h2,'Visible','off')
            h2 = plotEii_xy(gridEiiX,gridEiiY,Eii_map,mapsize,siimap);
            % update mean stress value, deviations from mean stress and legend values
            [Eii_dev,Eii_median] = median_dev_Sii(Eii_xymap,mapgrain_D,grain_comp,listDB);
            set(legend_medianVal,'String',['median strainrate = 1e',num2str(Eii_median,'%.1f'),' [s-1]']);
            set(legend_peakVal,'String',['max strainrate = 1e',num2str(max(max(Eii_xymap)),'%.1f'),' [s-1]']);
        case 7 % DM maps (T)
            % update deformation mechanisms based on new T
            if VR == 1 % voigt condition   
                [mech] = dmm_voigt(flowp,nXY,T,Daxis,Eiiaxis);
                cla(h2,'reset');
                set(h2,'Visible','off')
                h2 = plot_dmm(Siiaxis,Daxis,mech,mechcmap,T,'voigt');
            else % reuss's condition
                [mech] = dmm_reuss(flowp,nXY,T,Daxis,Siiaxis);
                cla(h2,'reset');
                set(h2,'Visible','off')
                h2 = plot_dmm(Siiaxis,Daxis,mech,mechcmap,T,'reuss');
            end  
        case 8 % XY deformation mechanism maps (switchable from Sii/xy and Eii/xy maps)
            % update deformation mechanism maps based on new T
            if VR == 1 
                % Voigt condition (isostrain-rate)
                % compute stress and stable deformation mechanism
                [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
                Sii_map  = zeros(mapsize(1),mapsize(2));
                for i = 1:numel(mapgrain_D)
                    Sii_map(mapgrain_ID{i}) = Sii_xymap(i);
                end
                % hide grain-boundaries
                for i = 1:mapsize(1)*mapsize(2)
                    if mapbound(i)==0
                        Sii_map(i) = nan;
                    end
                end
                % update mean stress value, deviations from mean stress and legend values
                [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
                set(legend_medianVal,'String',['mean stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
            else
                % Reuss condition (isostress)
                % compute strainrate and stable deformation mechanism
                [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
                Eii_map  = zeros(mapsize(1),mapsize(2));
                for i = 1:numel(mapgrain_D)
                    Eii_map(mapgrain_ID{i}) = Eii_xymap(i);
                end
                % hide grain-boundaries
                for i = 1:mapsize(1)*mapsize(2)
                    if mapbound(i)==0
                        Eii_map(i) = nan;
                    end
                end
                % update mean stress value, deviations from mean stress and legend values
                [Eii_dev,Eii_median] = median_dev_Eii(Eii_xymap,mapgrain_D,grain_comp,listDB);
                set(legend_medianVal,'String',['log10(',char(941),') = ',num2str(Eii_median,'%.1f'),' [s-1]']); 
            end
            cla(h2,'reset')
            set(h2,'Visible','off')
            h2 = plotDMM_xy(wmap,mapsize,mapbound,mapgrain_D,mapgrain_ID,eff_mech_xy,mapsize,mechcmap);
            set(legend_dmm(1:9),'Visible','on')      
        case 9 % SiiEii deviations maps
            if SiiEiiMap == 1
                % compute new stress and stress deviations
                [eff_mech_xy,Sii_xymap] = eff_mech_Sii_xy(mapgrain_D,grain_comp,listDB,T,Eii);
                % compute new Sii deviations
                [Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,mapgrain_D,grain_comp,listDB);
                % switch is from Sii_xy maps 
                Sii_devmap  = zeros(mapsize(1),mapsize(2));
                for i = 1:numel(mapgrain_D)
                    Sii_devmap(mapgrain_ID{i}) = Sii_dev(i);
                end
                % hide grain-boundaries
                for i = 1:mapsize(1)*mapsize(2)
                    if mapbound(i)==0
                        Sii_devmap(i) = nan;
                    end
                end
                cla(h2,'reset');
                set(h2,'Visible','off')
                h2 = plot_deviationSii(gridSiiX,gridSiiY,Sii_devmap,mapsize,cmapdiff);
                set(legend_medianVal,'String',['median stress = ',num2str(Sii_median,'%.1f'),' [MPa]']);
                % enable data dynamic display 
                dcm.UpdateFcn = @cursorXYmap_deviationSii;
            else
                % update strainrate and stable deformation mechanism based on new T
                [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
                % compute new strainrate deviations
                [Eii_dev,Eii_median] = median_dev_Eii(Eii_xymap,mapgrain_D,grain_comp,listDB);
                Eii_devmap  = zeros(mapsize(1),mapsize(2));
                for i = 1:numel(mapgrain_D)
                    Eii_devmap(mapgrain_ID{i}) = Eii_dev(i);
                end
                % hide grain-boundaries
                for i = 1:mapsize(1)*mapsize(2)
                    if mapbound(i)==0
                        Eii_devmap(i) = nan;
                    end
                end
                cla(h2,'reset');
                set(h2,'Visible','off')
                h2 = plot_deviationEii(gridEiiX,gridEiiY,Eii_devmap,mapsize,cmapdiff);
                % update median strainrate
                set(legend_medianVal,'String',['median log10(',char(941),') = ',num2str(Eii_median,'%.1f'),' [s-1]']);
                % enable data dynamic display 
                dcm.UpdateFcn = @cursorXYmap_deviationEii;
            end
    end
end

function slidSii(hObject,~)
    % the function returns the stress value as set by the end user
    Sii = get(hObject,'Value');
    % update dynamic legend properties
    set(legend_intvar1,'String',[char(963),' = ',num2str(Sii),' [MPa]']);  
    hold off
    set(m2h(1:5),'Checked','off');
    switch plottype
        case 2 % Eii/T projection space (Sii)
            if model==0
                [mech,Eiimax] = Eii_T_reuss(model,phase_name,Taxis,nXY,cgs,Sii);
            else
                [mech,Eiimax] = Eii_T_reuss(model,phase_name,Taxis,nXY,cgs,Sii,molp,cvol,mixr);
            end  
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = [char(941),' [s^{-1}]'];
            Yl = 'T [°C]';
            lab1 = ['grainsize = ',num2str(mean(cgs),'%.3f'),' mm'];
            lab2 = ['Mean stress ',char(963),num2str(Sii,'%.1f'),' MPa'];
            h2 = plot1d(lab1,lab2,Eiimax,Taxis-273.15,Eiiaxis,mech,Xl,Yl,model,'logX');
        case 4 % Eii/d projection space (Sii/T)
            % update stable deformation mechanism and strainrates
            [mech,Eiimin] = Eii_d_reuss(phase_name,Daxis,nXY,T,Sii);
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = ['Log_{10}(',char(941),') [s^{-1}]'];
            Yl = 'Log_{10}(d) [mm]';
            lab1 = ['Temperature = ',num2str(T-273.15,'%.1f'),'°C'];
            lab2 = [char(963),' = ',num2str(Sii,'%.1f'),'MPa'];
            h2 = plot1d(lab1,lab2,Eiimin,Daxis,Eiiaxis,mech,Xl,Yl,model,'logY','logX');  
        case 6 % Eii/xy map (Sii,T)
            % update strainrate and stable deformation mechanism based on new T
            [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
            Eii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Eii_map(mapgrain_ID{i}) = Eii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Eii_map(i) = nan;
                end
            end
            % clear previous plot and redraw the map
            cla(h2,'reset');
            colorbar off
            set(h2,'Visible','off')
            h2 = plotEii_xy(gridEiiX,gridEiiY,Eii_map,mapsize,siimap);
            % update mean stress value, deviations from mean stress and legend values
            [Eii_dev,Eii_median] = median_dev_Eii(Eii_xymap,mapgrain_D,grain_comp,listDB);
            set(legend_medianVal,'String',['median strainrate = 1e',num2str(Eii_median,'%.1f'),' [s-1]']);
            set(legend_peakVal,'String',['max strainrate = 1e',num2str(max(max(Eii_xymap)),'%.1f'),' [s-1]']);
        case 8 % XY deformation mechanism maps
            % update strainrate and stable deformation mechanism based on
            % new stress
            [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
            Eii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Eii_map(mapgrain_ID{i}) = Eii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Eii_map(i) = nan;
                end
            end
            cla(h2,'reset')
            set(h2,'Visible','off')
            h2 = plotDMM_xy(wmap,mapsize,mapbound,mapgrain_D,mapgrain_ID,eff_mech_xy,mapsize,mechcmap);
            set(legend_dmm(1:9),'Visible','on') 
            % update mean stress value, deviations from mean stress and legend values
            [Eii_dev,Eii_median] = median_dev_Eii(Eii_xymap,mapgrain_D,grain_comp,listDB);
            set(legend_medianVal,'String',['log10(',char(941),') = ',num2str(Eii_median,'%.1f'),' [s-1]']);   
        case 9
            % update strainrate and stable deformation mechanism based on new Sii
            [eff_mech_xy,Eii_xymap] = eff_mech_Eii_xy(mapgrain_D,grain_comp,listDB,T,Sii);
            % update map of grain stress
            Eii_map  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Eii_map(mapgrain_ID{i}) = Eii_xymap(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Eii_map(i) = nan;
                end
            end
            % compute new strainrate deviations
            [Eii_dev,Eii_median] = median_dev_Eii(Eii_xymap,mapgrain_D,grain_comp,listDB);
            Eii_devmap  = zeros(mapsize(1),mapsize(2));
            for i = 1:numel(mapgrain_D)
                Eii_devmap(mapgrain_ID{i}) = Eii_dev(i);
            end
            % hide grain-boundaries
            for i = 1:mapsize(1)*mapsize(2)
                if mapbound(i)==0
                    Eii_devmap(i) = nan;
                end
            end
            cla(h2,'reset');
            set(h2,'Visible','off')
            h2 = plot_deviationEii(gridEiiX,gridEiiY,Eii_devmap,mapsize,cmapdiff);
            % update median strain-rate
            set(legend_medianVal,'String',['median log10(',char(941),') = ',num2str(Eii_median,'%.1f'),' [s-1]']);
            % enable data dynamic display 
            dcm.UpdateFcn = @cursorXYmap_deviationEii;
    end
end

function slidd(hObject,~)
     % the function returns the grainsize value as set by the end user
     cgs = get(hObject,'Value');
     % change dynamic legend properties
     set(legend_intvar2,'String',['grain size = ',num2str(cgs),' [mm]']);      
     set(legend_intvar2,'Visible','on');
     hold off
     set(m2h(1:5),'Checked','off');
     switch plottype
         case 1  % Sii/T projection space; variables Eii(+d single-phase models only)
            % update stress based on new strain-rate 
            [mech,Siimin] = Sii_T_voigt(model,phase_name,Taxis,nXY,cgs,Eii);
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = [char(963),' [MPa]'];
            Yl = 'T [°C]';  
            lab1 = ['mean grainsize = ',num2str(mean(cgs),'%.3f'),'mm'];
            lab2 = ['Log_{10}(',char(941),')',num2str(log10(Eii),'%.1f'),'s^{-1}'];
            [h2] = plot1d(lab1,lab2,Siimin,Taxis,Siiaxis,mech,Xl,Yl,model);
         case 2 % Eii/T projection space; variables Sii (+d single-phase models only)
            [mech,Eiimax] = Eii_T_reuss(model,phase_name,Taxis,nXY,cgs,Sii);
            cla(h2,'reset');
            set(h2,'Visible','off')
            Xl = [char(941),' [s^{-1}]'];
            Yl = 'T [°C]';
            lab1 = ['grainsize = ',num2str(mean(cgs),'%.3f'),' mm'];
            lab2 = ['Mean stress ',char(941),')',num2str(Sii,'%.1f'),' MPa'];
            [h2] = plot1d(lab1,lab2,Eiimax,Taxis-273.15,Eiiaxis,mech,Xl,Yl,model,'logX');
     end
end

function reset_rheology(~,~)
    % reset all rheological parameters of current model (mechanical constraint,
    % projection space, mixture rule if composite model). The function 
    % also clear the graph and initialize the plot
    %------------------------------------------------------------------
    reset_rheology = questdlg(['Clicking OK will cancel the model and the data eventually displayed ',...
        'on the main plot! Do you want to reset the rheological model?'],'Reset rheological model','Yes','No','No');
    switch reset_rheology
        case 'Yes'
            model      = [];                  % composite (1), mono-phase model (0)
            phase_name = {};                  % name of phases in composite (this is = 1 for single phase materials!)
            phase_ID   = [];                  % ID of phase currently selected
            phase_pp   = [];                  % ID of phase currently selected for piezometric analysis
            flowp      = cell(3,1);           % flow law parameters
            mech       = [];                  % index of nodal effective deformation mechanism
            mixr       = '';
            cvol       = [];                  % phase proportions in composite
            cgs        = [];                  % median grain-size of phases in composite
            comp       = [];                  % name of composite
            molp       = [];                  % molar proportions of composite
            VR         = [];                  % initialize choice of Reuss/Voigt model
            plottype   = [];
            T          = [];
            Eii        = [];
            Sii        = [];
            SiiEiiMap  = [];
            % change menu properties
            set(m1h(6:7),'Checked','off','Enable','off');    % cancel & disable previous Voigt/Reuss selection
            set(m1h(1:4),'Checked','off','Enable','on');     % cancel previous model selection & enable new sel.
            set(m1h(9:13),'Checked','off','Enable','off');   % cancel & disable previous mixture rule 
            set(m1h(15:21),'Checked','off','Enable','off');  % cancel & disable previous projection space selection 
            % reset radio button group properties
            set(radio1,'Value',1);
            set(radio1,'Enable','off');
            set(radio2,'Enable','off')
            set(radio3,'Enable','off');
            set(radiotxt,'ForegroundColor',inact_txt);
            % change Data Zone properties
            set(modtxt1,'ForeGroundColor',inact_txt);
            set(mdtxt1_sel,'String',' ');
            set(pmtxt1,'ForeGroundColor',inact_txt);        
            set(pmtxt1_sel,'String',' ');
            set(modname,'ForeGroundColor',inact_txt);
            set(modname_sel,'String',' ');
            set(pmphase,'ForeGroundColor',inact_txt);
            set(pmphase_sel,'String',' ');
            set(modbehavior,'ForeGroundColor',inact_txt);
            set(modbehavior_sel,'String',' ');
            set(mixname,'ForeGroundColor',inact_txt);
            set(mixname_sel,'string',' '); 
            % reset slider properties
            set(tempsld3,'ForeGroundColor',inact_txt);
            set(tempsld2,'ForeGroundColor',inact_txt,'BackGroundColor',chkcol);
            set(tempsld1,'value',500,'BackGroundColor',chkcol);
            set(stresssld3,'ForeGroundColor',inact_txt);
            set(stresssld2,'ForeGroundColor',inact_txt,'BackGroundColor',chkcol);
            set(stresssld1,'value',100,'BackGroundColor',chkcol);
            set(strainsld3,'ForeGroundColor',inact_txt);
            set(strainsld2,'ForeGroundColor',inact_txt,'BackGroundColor',chkcol);
            set(strainsld1,'value',-12,'BackGroundColor',chkcol);
            set(dsld3,'ForeGroundColor',inact_txt);
            set(dsld2,'ForeGroundColor',inact_txt,'BackGroundColor',chkcol);
            set(dsld1,'value',0.1,'BackGroundColor',chkcol);
            cla(h2,'reset');
            set(h2,'Visible','off');
            set(legend_Siixy,'Visible','off');
            set(legend_Eiixy,'Visible','off');
            set(legend_dmm(1:9),'Visible','off') 
            set(legend_intvar1,'Visible','off');
            set(legend_intvar2,'Visible','off');
            set(legend_peakVal,'Visible','off');
            set(legend_medianVal,'Visible','off');
            set(legend_deviation_Sii,'Visible','off');
            set(legend_deviation_Eii,'Visible','off');
            h2 = initialize_plot;
        case 'No'
    end
end

%% (B) PLOTTING COMMANDS----------------------------------------------------
% this menu allow plotting some reference curve such as Byerlee envelope
% for various tensional states (normal, strike-slip and reverse), as well
% as paleopiezometric constraints. 
%---------------------------------------------------------------------------
function byen(~,~)
    % Plot Byerlee envelope for normal faults
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(bye_n,'visible','off');
    else
        set(gcbo, 'Checked', 'on');
        set(m2h(7),'Enable','on');         % enables resetting plot menu
        bye_n  = plot(h2,byerlee_norm,Taxis-273.15,'-','Color',[0.5 0.5 0.5]);
        % change axis labels
        Xl = [char(963),' [MPa]'];
        Yl = 'T [°C]';
        h2.YLabel.String = Yl;
        h2.XLabel.String = Xl;
        ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end        
end

function byess(~,~)
    % Plot Byerlee envelope for normal faults
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(bye_ss,'visible','off');
    else
        set(gcbo, 'Checked', 'on');
        set(m2h(7),'Enable','on');         % enables resetting plot menu
        bye_ss  = plot(h2,byerlee_strikeslip,Taxis-273.15,'--','Color',[0.5 0.5 0.5]);
        % change axis labels
        Xl = [char(963),' [MPa]'];
        Yl = 'T [°C]';
        h2.YLabel.String = Yl;
        h2.XLabel.String = Xl;
        ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end        
end

function byer(~,~)
    % Plot Byerlee envelope for normal faults
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(bye_r,'visible','off');
    else
        set(gcbo, 'Checked', 'on');
        set(m2h(7),'Enable','on');         % enables resetting plot menu
        bye_r  = plot(h2,byerlee_reverse,Taxis-273.15,':','Color',[0.5 0.5 0.5]);
        % change axis labels
        Xl = [char(963),' [MPa]'];
        Yl = 'T [°C]';
        h2.YLabel.String = Yl;
        h2.XLabel.String = Xl;
        ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end        
end

function qtzcreep(~,~)
    % Plot quartzite dislocation creep strength profiles for strain-rates
    % between 1e-8 to 1e-16
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(qtz_dc(:),'visible','off');
        set(t1,'Visible','off');
        set(t2,'Visible','off');
        set(t3,'Visible','off');
        set(t4,'Visible','off');
        set(t5,'Visible','off');
    else
        set(gcbo, 'Checked', 'on');
        set(m2h(7),'Enable','on');         % enables resetting plot menu
        if isempty(plottype)==1 || plottype==1
            eiiq = [1e-8 1e-10 1e-12 1e-14 1e-16];
            quartzite_creep    = zeros(nXY,numel(eiiq));
            for kk = 1:numel(eiiq)
                quartzite_creep(:,kk) = (eiiq(kk)./(1.198e-8*exp(-152e3./(8.314.*Taxis)))).^0.25;
            end
            % get coordinates of labels for quartz and olivine reference curves
            qtz_dc_label(1) = (-152e3/8.314)/(log(1e-8/(1.198e-8*250^4)));
            qtz_dc_label(2) = (-152e3/8.314)/(log(1e-10/(1.198e-8*500^4)));
            qtz_dc_label(3) = (-152e3/8.314)/(log(1e-12/(1.198e-8*300^4)));
            qtz_dc_label(4) = (-152e3/8.314)/(log(1e-14/(1.198e-8*500^4)));
            qtz_dc_label(5) = (-152e3/8.314)/(log(1e-16/(1.198e-8*350^4)));
            qtz_dc = zeros(1,numel(eiiq));
            for i = 1:numel(eiiq)
                qtz_dc(i)  = plot(h2,quartzite_creep(:,i),Taxis-273.15,'-b','LineWidth',1);
                hold on
            end
            t1 = text(250,qtz_dc_label(1)-273.15,'1e-8 s^{-1}','BackgroundColor','w');
            t2 = text(500,qtz_dc_label(2)-273.15,'1e-10 s^{-1}','Color','none');
            t3 = text(300,qtz_dc_label(3)-273.15,'1e-12 s^{-1}','BackgroundColor','w');
            t4 = text(500,qtz_dc_label(4)-273.15,'1e-14 s^{-1}','Color','none');
            t5 = text(350,qtz_dc_label(5)-273.15,'1e-16 s^{-1}','BackgroundColor','w');
            % change axis labels
            Xl = [char(963),' [MPa]'];
            Yl = 'T [°C]';
            h2.YLabel.String = Yl;
            h2.XLabel.String = Xl;
            ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        else
            siiq = [1 10 100 200 400];
            quartzite_creep    = zeros(nXY,numel(siiq));
            for kk = 1:numel(siiq)
                quartzite_creep(:,kk) = log10(1.198e-8*siiq(kk)^4*exp(-152e3./(8.314.*Taxis)));
            end
            % get coordinates of labels for quartz and olivine reference curves
            qtz_dc_label(1) = (-152e3/8.314)/log(1e-14/(1.198e-8*1^4));
            qtz_dc_label(2) = (-152e3/8.314)/log(1e-12/(1.198e-8*10^4));
            qtz_dc_label(3) = (-152e3/8.314)/log(1e-11/(1.198e-8*100^4));
            qtz_dc_label(4) = (-152e3/8.314)/log(1e-10/(1.198e-8*200^4));
            qtz_dc_label(5) = (-152e3/8.314)/log(1e-9/(1.198e-8*400^4));
            qtz_dc = zeros(1,numel(siiq));
            for i = 1:numel(siiq)
                qtz_dc(i)  = plot(h2,quartzite_creep(:,i),Taxis-273.15,'-b','LineWidth',1);
                hold on
                t1 = text(-14,qtz_dc_label(1)-273.15,'1 MPa','BackgroundColor','w');
                t2 = text(-12,qtz_dc_label(2)-273.15,'10 MPa','BackgroundColor','w');
                t3 = text(-11,qtz_dc_label(3)-273.15,'100 MPa','BackgroundColor','w');
                t4 = text(-10,qtz_dc_label(4)-273.15,'200 MPa','BackgroundColor','w');
                t5 = text(-9,qtz_dc_label(5)-273.15,'400 MPa','BackgroundColor','w');
            end
            % change axis labels
            Xl = ['Log_{10}(',char(941),') [s^{-1}]'];
            Yl = 'T [°C]';
            h2.YLabel.String = Yl;
            h2.XLabel.String = Xl;
            ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        end
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end        
end

function olcreep(~,~)
    % Plot olivine dislocation creep strength profiles for strain-rates
    % between 1e-8 to 1e-16
    %------------------------------------------------------------------
    if strcmp(get(gcbo, 'Checked'),'on')
        set(gcbo, 'Checked', 'off');
        set(ol_dc(:),'visible','off');
    else
        set(gcbo, 'Checked', 'on');
        set(m2h(7),'Enable','on');         % enables resetting plot menu
        if isempty(plottype)==1 || plottype==1
            eiiq = [1e-8 1e-10 1e-12 1e-14 1e-16];
            olivine_creep      = zeros(nXY,numel(eiiq));
            for kk = 1:numel(eiiq)
                olivine_creep(:,kk)   = (eiiq(kk)./(11e+4*exp(-530e3./(8.314.*Taxis)))).^0.25;
            end
            % get coordinates of labels for quartz and olivine reference curves
            ol_dc_label(1) = (-530e3/8.314)/(log(1e-8/(11e+4*300^4)));
            ol_dc_label(2) = (-530e3/8.314)/(log(1e-10/(11e+4*500^4)));
            ol_dc_label(3) = (-530e3/8.314)/(log(1e-12/(11e+4*400^4)));
            ol_dc_label(4) = (-530e3/8.314)/(log(1e-14/(11e+4*500^4)));
            ol_dc_label(5) = (-530e3/8.314)/(log(1e-16/(11e+4*500^4)));
            ol_dc = zeros(1,numel(eiiq));
            for i = 1:numel(eiiq)
                ol_dc(i)  = plot(h2,olivine_creep(:,i),Taxis-273.15,'-g','LineWidth',1);
                hold on
            end
            t1 = text(300,ol_dc_label(1)-273.15,'1e-8 s^{-1}','BackgroundColor','w');
            t2 = text(500,ol_dc_label(2)-273.15,'1e-10 s^{-1}','Color','none');
            t3 = text(400,ol_dc_label(3)-273.15,'1e-12 s^{-1}','BackgroundColor','w');
            t4 = text(500,ol_dc_label(4)-273.15,'1e-14 s^{-1}','Color','none');
            t5 = text(500,ol_dc_label(5)-273.15,'1e-16 s^{-1}','BackgroundColor','w');
            % change axis labels
            Xl = [char(963),' [MPa]'];
            Yl = 'T [°C]';
            h2.YLabel.String = Yl;
            h2.XLabel.String = Xl;
            ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        else
            siiq = [1 10 100 200 400];
            olivine_creep      = zeros(nXY,numel(siiq));
            for kk = 1:numel(siiq)
                olivine_creep(:,kk) = log10(11e4*siiq(kk)^4*exp(-530e3./(8.314.*Taxis)));
            end
            % get coordinates of labels for quartz and olivine reference curves
            ol_dc_label(1) = (-530e3/8.314)/log(1e-14/(1.198e-8*1^4))-273.15;
            ol_dc_label(2) = (-530e3/8.314)/log(1e-12/(1.198e-8*10^4))-273.15;
            ol_dc_label(3) = (-530e3/8.314)/log(1e-11/(1.198e-8*100^4))-273.15;
            ol_dc_label(4) = (-530e3/8.314)/log(1e-10/(1.198e-8*200^4))-273.15;
            ol_dc_label(5) = (-530e3/8.314)/log(1e-9/(1.198e-8*400^4))-273.15;
            ol_dc = zeros(1,numel(siiq));
            for i = 1:numel(siiq)
                ol_dc(i)  = plot(h2,olivine_creep(:,i),Taxis-273.15,'-b','LineWidth',1);
                hold on
                t1 = text(-14,ol_dc_label(1),'1 MPa','BackgroundColor','w');
                t2 = text(-12,ol_dc_label(2),'10 MPa','BackgroundColor','w');
                t3 = text(-11,ol_dc_label(3),'100 MPa','BackgroundColor','w');
                t4 = text(-10,ol_dc_label(4),'200 MPa','BackgroundColor','w');
                t5 = text(-9,ol_dc_label(5),'400 MPa','BackgroundColor','w');
            end
            % change axis labels
            Xl = ['Log_{10}(',char(941),') [s^{-1}]'];
            Yl = 'T [°C]';
            h2.YLabel.String = Yl;
            h2.XLabel.String = Xl;
            ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        end
        % activate object properties manager 
        hlines = findall(h2,'Type','line');
        % Attach the context menu to each line
        for line = 1:length(hlines)
            set(hlines(line),'uicontextmenu',hcmenu)
        end
    end        
end

function cancel_plot(~,~)
    % clear the main plot without deleting variables and rheological model 
    %------------------------------------------------------------------
    % change menu properties 
    set(m2h(7),'Enable','off');        % reset graph off
    set(m2h(1:5),'Checked','off');
    cla(h2,'reset');
    set(h2,'Visible','off');
    plottype = [];
    h2 = initialize_plot;
end

%% (C) PALEOPIEZOMETRY FUNCTIONS------------------------------------------
function new_cal(~,~)
    % the function enables loading a new piezometric calibration for either one
    % of the minerals already included in the DB, or for a new mineral specie
    %----------------------------------------------------------------------
    uiwait(msgbox(['Piezometric calibration of the form:              ';...
        '             sigma = B*exp(Q/T)*d^-m              ';...
        'Please, set Q = 0  for calibration not dependent  ';...
        'on  T.                                            '],'set DB','modal'));
    namecal = inputdlg('Set name of new piezometer','New calibration'); 
    options.Interpreter = 'tex';
    options.Resize      = 'on';
    options.WindowStyle = 'normal';
    prompt = {'mineral: ','Material constant [MPa mm_p]:  ','Temperature-dependent term [K]','stress exponent'};
    inpdlgtitle = cell2mat(namecal);
    defans = {'quartz','500e-3','0','-0.9'};
    par = inputdlg(prompt,inpdlgtitle,1,defans,options);                                
    save(['DB_calibrations\',char(namecal)],'par')

    % update piezometric calibration database as new calibration has been
    % added
    listCAL = dir('DB_calibrations\*.mat');
    [r2,~]  = size(listCAL);
    mincal = cell(1,r2);
    indV = 1:r2;
    for IDcal = 1:r2
        cnamewithext    = listCAL(IDcal).name;
        cnamewithoutext = strsplit(cnamewithext,'.');
        ccurrname       = cnamewithoutext(1);
        listcal(IDcal)   = ccurrname;
        calpar = load(['DB_calibrations\',cnamewithext]);
        calpar.par(1);
        mincal(IDcal) = calpar.par(1);
    end
end

function mod_cal(~,~)
    % the function allow modifying one of the existent piezometric
    % calibrations from the database
    %----------------------------------------------------------------------
    [selcal,~] = listdlg('ListString',listcal,'Name','Select calibration to be modified: ','SelectionMode','single');
    if isempty(selcal)==0
        options.Interpreter = 'tex';
        caledit     = load(['MYflow\DB_calibrations\',listcal{selcal},'.mat']);
        valedit = caledit.par;                   % extract parameters of selected calibration
        prompt = {'mineral: ','Material constant [MPa mm_p]:  ','Temperature-dependent term [K]','stress exponent'};
        titledit = cell2mat(valedit(1));
        uiwait(msgbox(['Modify piezometric calibration of the form:       ';...
        '             sigma = B*exp(Q/T)*d^-m              ';...
        'Please, set Q = 0  for calibration not dependent  ';...
        'on  T.                                            '],'set DB','modal'));
        par = inputdlg(prompt,titledit,1,valedit,options);
        save(['DB_calibrations\',listcal{selcal}],'par');
    end
end

function d_sing(~,~)
    % the function allow to add single grainsize/T data for paleostress 
    % analysis. Paleopiezometric analysis could run by itself without 
    % previously setting the compositional model (i.e., single
    % phase/composite).
    %------------------------------------------------------------------
    [phase_ID,~] = listdlg('ListString',listDB,'Name','Select phase: ','SelectionMode','single');
    phase_name = listDB(phase_ID);
    % show selection message
    figw=warndlg({'Single phase model authomatically selected! '},'Model setup');
    uiwait(figw);
    % change datazone properties
    set(modname,'ForeGroundColor',active_txt);
    set(modname_sel,'String',cell2mat(phase_name),'ForeGroundColor',active_txt);
    phase_pp = phase_name;
    prompt_addsingle = {['Size [',char(956),'m]'],'d err. ','T [°C]','T err. [°C]'};
    name_addsingle = 'Set grainsize/T data';
    numlines = [1 47];
    defaultanswer = {'20','0','700','0'};
    SINGLE_DATA = inputdlg(prompt_addsingle,name_addsingle,numlines,defaultanswer);
    if isempty(SINGLE_DATA)==0
        newrow = str2double(SINGLE_DATA);
        for SDC = 1:4
            if isnan(newrow(SDC))
                newrow(SDC)= 0;
            end
        end
        datin = datin+1;                       % update number of data in current dataset
        dat(datin,1) = {false};                % update checkbox value for new data entry
        dat(datin,2) = phase_pp;               % update name of current mineral
        dat(datin,4) = {round(newrow(1))};     % update mean grain size
        dat(datin,5) = {round(newrow(2))};     % update grainsize error
        dat(datin,6) = {round(newrow(3))};     % update mean T  
        dat(datin,7) = {round(newrow(4))};     % update T error
        set(tab1,'Data',dat);   
    else
    end
    % change menu properties
    set(m3h(3:4),'Enable','on');
end

function d_imp(~,~)
    % the function allow to import grainsize/T data from .xlsx
    % spreadsheets. d_imp has the same behavior of function  d_sing, so
    % the user is asked to select a mineral phase from the data base 
    % and the model type is authomatically set to 'single phase', as
    % long as model=[]
    %------------------------------------------------------------------
    [phase_ID,~] = listdlg('ListString',listDB,'Name','Select phase: ','SelectionMode','single');
    phase_name = listDB(phase_ID);
    % change datazone properties
    set(modname,'ForeGroundColor',active_txt);
    set(modname_sel,'String',cell2mat(phase_name),'ForeGroundColor',active_txt);
    phase_pp = phase_name;
    [FileName,PathName,~] = uigetfile([currentFolder,'\paleopiezometry\example.xlsx']);
    if FileName~=0
        [data] = xlsread([PathName,FileName]);
        [rxls,cxls] = size(data);
        datin = datin+rxls;    % update number of data in current dataset
        hold on
        if cxls==2
            % data as single points without errors
            for j = datin-rxls+1:datin
                dat(j,1) = {false};
                dat(j,2) = phase_pp;
                dat(j,4) = {round(data(j-(datin-rxls),1))};
                dat(j,5) = {0};
                dat(j,6) = {round(data(j-(datin-rxls),2))};
                dat(j,7) = {0};
            end
        elseif cxls==4
            % d,T data + errors (2 sigma)
            diammean = data(:,1);
            diamSD   = data(:,2);
            Tmean    = data(:,3);
            TSD      = data(:,4);
            for b = datin-rxls+1:datin
                dat(b,1) = {false};
                dat(b,2) = phase_pp;
                dat(b,4) = {round(diammean(b-(datin-rxls)))};
                dat(b,5) = {round(diamSD(b-(datin-rxls)))};
                dat(b,6) = {round(Tmean(b-(datin-rxls)))};
                dat(b,7) = {round(TSD(b-(datin-rxls)))};
            end
        else
            % makes warning message if data are not properly formatted
            warndlg({'WARNING: data should be organized as: ';...
                ['1. grain diameter [',char(956),'m]/T[°C]'];...
                ['2. grain diameter [',char(956),'m]/min T[°C]/max T[°C].'];...
                'Check data format'}); 
        end 
        set(tab1,'Data',dat);
        % change menu properties
        set(m3h(3:4),'Enable','on');
    end
end

function calcPP(~,~)
% Specify mineral phase to calculate differential stress from
% 1. show available piezometric relationships.
% 2. select calibration
% 3. specify marker color for current (sub)-set of
% paleopiezometric data
%------------------------------------------------------------------     
    pdat = cell2mat(dat(datcalALL+1:datin,4:7));
    % calculate and plot paleostress values based on latest uploaded minerals
    % NOTE that adding one grain-size at once require to specify a
    % calibration each time
    if isempty(pdat)==0
        % some grainsize data has been uploaded from external file or
        % manually set up through the interactive data table
        switch cell2mat(phase_pp)
            case {'anorthite dry','anorthite wet','anorthite','albite'}
                % load plg piezometric relationshipsy
                indAct = indV(strcmp(mincal,'plagioclase')); % index of available piezometers
                piezoname = 'plagioclase';
            case {'feldspar','Kfeldspar'}
                % load plg piezometric relationshipsy
                indAct = indV(strcmp(mincal,'Kfeldspar')); % index of available piezometers
                piezoname = 'Kfeldspar';
            case {'quartz','a-quartz','B-quartz'}
                % load qtz piezometers
                indAct = indV(strcmp(mincal,'quartz')); % index of available piezometers
                piezoname = 'quartz';
            case {'calcite','aragonite'}
                % load Cal piezometers
                indAct = indV(strcmp(mincal,'calcite')); % index of available piezometers
                piezoname = 'calcite';
            case {'olivine','fayalite','forsterite','olivine_wet','olivine_dry'}
                % load ol piezometers
                indAct = indV(strcmp(mincal,'olivine')); % index of available piezometers
                piezoname = 'olivine';
            case {'diopside','pyroxene','jadeite','orthopyroxene','clinopyroxene'}
                % load opx piezometers
                indAct = indV(strcmp(mincal,'pyroxene')); % index of available piezometers
                piezoname = 'pyroxene';
            otherwise
                % no available piezometric relationships for this mineral -
                % update the calibration DB, if possible!
                indAct = [];
                piezoname = [];
                warndlg(['Sorry, there''s no available piezometric calibration for ',phase_pp,'. Update requested'])
        end
        if isempty(indAct)==0
            listAct = listcal(indAct);                   % list of available piezometers 
            [calibration,~] = listdlg('ListString',listAct,'Name','Select calibration: ','SelectionMode','single');
            set(piezo_name_sel,'String',piezoname); 
            calAct = load(['DB_calibrations\',cell2mat(listAct(calibration)),'.mat']);
            B   = str2double(calAct.par(2));
            nud = str2double(calAct.par(4));
            Qp  = str2double(calAct.par(3));
            
            % change datazone properties
            set(piezo_name,'ForegroundColor',active_txt);
            set(cal_name,'ForegroundColor',active_txt);
            set(piezo_name_sel,'ForegroundColor',active_txt);
            set(cal_name_sel,'ForegroundColor',active_txt);
            set(cal_name_sel,'String',listAct(calibration)); 

            % get number of data stress should be calculated for
            datcal = datin-datcalALL;
            datcalALL = datcalALL+datcal;

            % ask for marker color and update color in data table
            color_string = {'white','red','blue','green','yellow','black','cyan','gray'};
            markcolor = listdlg('ListString',color_string,'Name','Set marker''s color','SelectionMode','single');
            % compute average paleostress (all calibrations but the ones based on

            % Twiss's model (and Mehl & Hirth 2008) which use grain size in mm!)
            switch cell2mat(listAct(calibration))
                case {'Twiss (1977) - quartz','Twiss (1977) - olivine','Twiss (1977) - plagioclase','Twiss (1977) - calcite','Mehl & Hirth (2008) - plagioclase'}
                    pdat(:,1) = pdat(:,1)./1000;         % convert gran-size in mm
                    pdat(:,2) = pdat(:,2)./1000;
            end

            pstress_mean = B.*exp(Qp./pdat(:,3)).*pdat(:,1).^nud; 
            pstress_max = B.*exp(Qp./pdat(:,3)).*(pdat(:,1)-pdat(:,2)).^nud;
            pstress_min = B.*exp(Qp./pdat(:,3)).*(pdat(:,1)+pdat(:,2)).^nud; 
            pTmin       = pdat(:,3)-pdat(:,4);
            pTmax       = pdat(:,3)+pdat(:,4); 
            % add calculated mean paleostress values to the table
            indPstress = 0;
            for i = datin-datcal+1:datin
                indPstress = indPstress+1;
                dat(i,8) = {pstress_mean(indPstress)};
                dat(i,3) = color_string(markcolor);
            end
            % change table data properties
            set(tab1,'Data',dat);
            
            % display authomatically piezometric data as they are calculated.
            % Plot paleostress data as error ellipses in T-sigma space (the
            % shape and dimension of ellipses are proportional to the errors
            % relative to grainsize and T
            indST = 0;
            for g = datin-datcal+1:datin
                indST = indST+1;
                hold on
                if pstress_min(indST)==pstress_mean(indST) || pTmin(indST)==pTmax(indST)
                    ellcoo = [pstress_min(indST) pTmin(indST) 7 13];
                else
                    ellcoo = [pstress_min(indST) pTmin(indST) pstress_max(indST)-pstress_min(indST) pTmax(indST)-pTmin(indST)];             
                end
                hmdata(datin-datcal+indST) = rectangle(h2,'Position',ellcoo,'Curvature',[1 1],'FaceColor',color_string{markcolor},'EdgeColor','k');
            end 
        else
        end
        % change table properties
        coleditable = [true,false,false,false,false,false,false,false];
        set(tab1,'ColumnEditable',coleditable);
        for i = datin-datcal+1:datin
            dat(i,1) = {true};
        end
        set(tab1,'Data',dat);
        % change axis labels
        Xl = [char(963),' [MPa]'];
        Yl = 'T [°C]';
        h2.YLabel.String = Yl;
        h2.XLabel.String = Xl;
        set(m2h(7),'Enable','on');         % enables resetting plot menu
        ylim(h2,[Taxis(1)-273.15,Taxis(end)-273.15]);
        % activate object properties manager 
        hmark = findall(h2,'Type','line');
        % Attach the context menu to each marker
        for mark = 1:length(hmark)
            set(hmark(mark),'uicontextmenu',hcmenu)
        end
        % change menu properties
        set(m3h(5),'Enable','on');
    else
        % there are no grain-size data 
        warndlg('MYflow requires uploading new grain size data to calculate paleostress!');
    end
end

function hidePP(~,~)
    % hide paleopiezometric data from the plot
    % the function doesn't reset the piezometric dataset.
    %------------------------------------------------------------------
    set(hmdata,'Visible','off');
    % change menu properties
    set(m3h(5),'Checked','on');
end

function cancel_pp(~,~)
    % reset paleopiezometric data (calculations, results, selected
    % calibration);
    %------------------------------------------------------------------   
    set(hmdata,'Visible','off');           % clear plot
    hmdata               = [];             % handles to paleostress data added to current plot
    pstress_mean         = [];             % mean paleostress values computed from grainsize dataset (REQUIRED!)
    pstress_min          = [];             % minimum paleostress values (if available)
    pstress_max          = [];             % maximum paleostress values (if available)
    pTmin                = [];             % minimum temperature estimate associated to grainsize data (if available)
    pTmax                = [];             % maximum temperature associated to grainsize data (if available)
    datin                = 0;              % counter of uploaded grainsize/T data
    datcal               = 0;              % counter of calculated piezometric data (datin <= datcal!!)
    datcalALL            = 0; 
    dat                  = cell(999,8);
    % change table properties
    set(tab1,'Data',dat);
    % change menu properties
    set(m3h(3:5),'Enable','off');
    set(m3h(3:5),'Checked','off');
    % reset Data Zone 
    set(cal_name,'ForeGroundColor',inact_txt); % make inactive calibration name (datazone)
    set(cal_name_sel,'String',' ');            % hide selected calibration name
    set(piezo_name_sel,'String',' ');
    set(piezo_name,'ForegroundColor',inact_txt);
end

function piezo_data(~,EventData)
    % manage single piezometric data
    val = EventData.NewData;
    pos = EventData.Indices;
    if isempty(hmdata)==0 % check a plot has been produced
        if val==0         % check user's choice to display/not display the result
            set(hmdata(pos(1)),'Visible','off');
        elseif val==1
            set(hmdata(pos(1)),'Visible','on');
        end
    end
end

%% GUIDE function
function help(~,~)
    open guide.pdf
end

end
