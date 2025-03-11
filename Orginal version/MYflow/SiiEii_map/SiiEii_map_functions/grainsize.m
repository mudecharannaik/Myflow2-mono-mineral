function[grain_ID,grain_comp,grain_D,isize,w,bound] = grainsize(model,listDB,currentFolder,phase_name,phase_ID) 
% the function allow to load a phase map of either monophase or composite
% aggregates composed of N phases. Unknown phases are returned as not
% indexed pixels.
% RETURNS:
% grain_ID   = [1-by-N of grains] cell array. Each cell corresponds to a detected
%              grain and contains the GLOBAL indices of pixels belonging to
%              that grain
% grain_comp = [1-by-N of grains]composition vector. Each element of the 
%              vector provides the index of the mineral contituting that 
%              grain
% grain_D    = 1-by-N vector of grain equivalent-circle diameters [um]
% grain_xy   = 2-by-N matrix of grain centroid coordinates
% isize      = [a b] vector, size of indexed image
% w          = scalar, width of image [mm]
%--------------------------------------------------------------------------        
if model==0
    % single phase models
    N = 1;
else
    % composite models
    N = numel(phase_name);
end
msg1 = msgbox({[char(963),' and ',char(941),' colormaps require uploading '];...
    'color-coded phase map composed of at least 1 mineral specie.';...
    'The selected phase(s) must be already included in the mineral database.';...
    'Make phase boundaries BLACK to optimize grain-size analysis!'},'Phase Map Setup');
uiwait(msg1);
[I,path] = uigetfile([currentFolder,'\SiiEii_map\microstructural_map\*.png'],'load figure');
prompt1 = {'Set width of phasemap [mm]:'};
options.Interpreter = 'tex';
options.Resize      = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
inpdlgtitle = ['Set phasemap ',I,': '];
defans = {'5'};
w = str2double(inputdlg(prompt1,inpdlgtitle,1,defans,options));

% load RGB image of some microstructure; could be either monophase
% aggregate or composites made of up to 5 mineral phases
im = imread([path,I]);                           % RGB image
imb = rgb2ind(im,N+1);                           % indexed image
imb = -int8(imb);                                % negative indexes to avoid unwanted matching with minerals ID
isize = size(imb);
bound = ones(isize(1),isize(2));
for i = 1:isize(1)
    for j = 1:isize(2)
        if imb(i,j)==0
            % extract grain boundaries (black)
            bound(i,j)=0;
        end
    end
end
% find grains in phase map
cc = bwconncomp(bound,4);
    
% extract number of grains
ngrains = cc.NumObjects;

if model==0
    %% monophase models do not require phase selection - directly show the number of grains detected
    phaseID = phase_ID;
    phind = -1;
    f1 = figure('Name','Phase Map','NumberTitle','off');
    imshow(im);
    msg = msgbox(['MYflow detected a total of ',num2str(ngrains),' grains in phasemap ',I]);
    uiwait(msg);
else
%% composites require phase selection
    % composites require specificating phase composition
    % select one reference point within each phase to setup colormap
    phaseID = zeros(N,1);              % initialize array of ID to selected mineral phases
    msg2 = msgbox({['Set the composition of phase map ',I];...
        '1. left click at the center of one grain on the phasemap';...
        '2. set composition of selected phase by loading the appropriate mineral';...
        '   from the mineral database.';...
        '3. repeat click/selection until all phases are classified';...
        ['The current phasemap requires ',num2str(N),' selections!']},'Phase Map setup');
    uiwait(msg2);
    f1 = figure('Name','Phase Map','NumberTitle','off');
    imshow(im);
    % get coordinate of reference point of current mineral phase
    phind = zeros(N,1);
    for j = 1:N 
        % select mineral phase from the phasemap
        [x,y] = ginput(1);
        % assign composition  to selected mineral phase
        [local_ind,~] = listdlg('ListString',phase_name,'Name','Set mineral ','SelectionMode','single');
         % get global index of currently selected mineral phase
        ind = 1:numel(listDB);                         
        global_ind = ind(strcmp(phase_name(local_ind),listDB));   
        if sum(phaseID~=global_ind)==N
            phaseID(j) = global_ind;
            % detect pixel ID (in the indexed image reference frame, values from 0 to N) 
            % corresponding to the selected point of coordinates x,y. 
            % Use bisection method
            xmin = 1;
            xmax = isize(2);
            while (xmax-xmin)>1
                xi = double(int16((xmax+xmin)./2-0.5));
                if x<xi
                    xmax = xi;
                else
                    xmin = xi;
                end
            end
            xi = xmin;
            if xi<1
                xi = 1;
            end
            if xi>isize(2)-1
                xi = isize(2)-1;
            end
            X = xi;
            ymin = 1;
            ymax = isize(1);
            while (ymax-ymin)>1
                yi = double(int16((ymax+ymin)./2-0.5));           
                if y<yi
                    ymax = yi;
                else
                    ymin = yi;
                end
            end
            yi = ymin;
            if yi<1
                yi = 1;
            end
            if yi > isize(1)-1
                yi = isize(1)-1;
            end
            Y = yi;
            phind(j) = imb(Y,X);
        else
            % menage duplicate phase selection
            msg = msgbox(['WARNING: the phase',listDB(local_ind),' has been already selected.','Chose another phase ','or hit cancel']);
            uiwait(msg);
        end
    end
end

if sum(phaseID==0)~=0
    % The user hit 'cancel' button before selecting the required number of
    % phases
    msg = msgbox(['WARNING: at least one phase ','is required for colormaps']);
    uiwait(msg);
    grain_ID = [];
    grain_comp = [];
    grain_D = []; 
else
    close(f1);
    % analyze grainsize 
    % change ID of indexed image to make them consistent with codes in the
    % mineral database
    for i = 1:N  
        imb(imb==phind(i)) = phaseID(i);
    end
    % detect composition of each grain and assign the corresponding ID 
    grain_comp = zeros(1,ngrains);
    for i = 1:ngrains
        grain_comp(i) = mode(imb(cc.PixelIdxList{i})); 
    end
    grain_ID = cc.PixelIdxList; 
    
    % get pixel resolution [mm^2]
    pixres = ((w)/isize(2))^2;

    % compute area-based statistics and get grainsize data 
    graindata = regionprops(cc,'basic');
    grainArea = [graindata.Area].*pixres;
    grain_D = 2*sqrt(grainArea./pi);
    
    % visualize grain-size statistics
    f2 = figure('Name','Grain size statistics','NumberTitle','off');
    colorpalette = {'k','b','m','g','r','c','y'};
    phase_area = zeros(1,N);
    valid_phase = 0;
    for i = 1:N
        if strcmp(listDB(phaseID(i)),'not defined')==0
            valid_phase = valid_phase+1;
            phase_area(i) = sum(grainArea(grain_comp==phaseID(i))); % area fraction of each phase (excluding not defined!)
        end
    end
    phase_total = sum(phase_area);
    for i = 1:N
        if strcmp(listDB(phaseID(i)),'not defined')==0
            phase_diameter = grain_D(grain_comp==phaseID(i));
            area_fraction = phase_area(i)/phase_total; % area fraction of each phase 
            subplot(ceil(valid_phase/3),3,i)
            nbins = round(numel(phase_diameter)/5);
            h1 = histogram(phase_diameter/area_fraction,nbins,'Normalization','Probability');
            
            title([listDB(phaseID(i)),', N = ',num2str(numel(phase_diameter))],'FontSize',7);
            hold on
            xav = median(phase_diameter);
            yav = max(h1.Values);
            h1.FaceColor = colorpalette{i};
            %xmedian = round(numel(phase_diameter)/2)+1;
            %ymedian = max(phase_diameter);
            %plot([xmedian,xmedian],[0,ymedian],'-r','LineWidth',2);
            text(xav,yav,['d = ',num2str(median(phase_diameter).*1000,'%.1f'),' [um]'],...
                'BackgroundColor','w','Color',colorpalette{i},'FontSize',7);
            %text(xmedian,ymedian/2.75,['\phi = ',num2str(area_fraction,'%.2f'),' [A_{phase}/A_{tot}]'],...
            %    'BackgroundColor','w','Color',colorpalette{i},'FontSize',7);
        end
    end 
    msg = msgbox('Phase Map successfully created!');
    uiwait(msg);
    if isvalid(f2)==1
        close(f2)
    end
end
end