function[mech,Siimin] = Sii_T_voigt(model,phase_name,Taxis,nXY,gs,Eii,varargin)
% the function Sii_d_voigt.m calculates the effective deformation
% mechanism at nodal points of the 2D grid nXYZ-by-nXYZ.
% inputs:
% phase_name = could be either the name of phase in single-phase materials,
%              or the name of a composite. In this latter case, the
%              function requires a total of 10 inputs, the last 5 required
%              only for composites
% Taxis      = temperature vector
% nXY        = number of gridpoints in X,Y directions
% gs         = grain size values (this is just one scalar for monophase
%              models, and a vector for composites
% Eii        = reference strain rate
% vararg{1}  = molar proportion in composites
% vararg{2}  = volume proportion of phases in composites
% vararg{3}  = ID of selected mixture rule (compatible with Sii/T
%              projection: 1-mean, 2-arithmetic, 3-harmonic, 4-huet, 5-hobbs) 
% returns:
% mech   = 1-by-nXYZ vector of index to active 
%          deformation mechanism, specified this way: 1. dislocation creep, 
%          2. diffusion creep, 3. GBS. 
% Siimin = 1-by-nXYZ vector of stress values for stable deformation
%          mechanism (at any given point)
%--------------------------------------------------------------------------
if nargin~=6 && model==0
    msgbox('Error! Function Sii_T_voigt.m requires 6 input to be used for monophase');
end
if nargin~=9 && model==1
    msgbox('Error! Function Sii_T_voigt.m requires 9 inputs to be used for composites - including compositionally homogeneous rocks with variable grain size');
end
R = 8.314;                                       % gas constant
mech = zeros(1,nXY);                             % initialize arrays 
Siimin = zeros(1,nXY);

if model==0 || numel(phase_name)==1 
    % monophase models
    dc = zeros(4,1);                                 % dislocation creep parameters
    diff = zeros(4,1);                               % diffusion creep parameters
    gbs = zeros(4,1);                                % GBS parameters
    mat = load(['DB_mineral_parameters\',char(phase_name),'.mat']);
    parameters  = mat.par;                           % extract material parameters
    % dislocation creep, array dc
    dc(1) = parameters(2,1);                         % pre-exponential constant 
    dc(2) = 4;                                       % stress exp. 
    dc(3) = 0;                                       % grain-size exp. 
    dc(4) = parameters(2,2);                         % activation energy 
    % diffusion creep
    diff(1) = (82*parameters(1,3)^3*...              % pre-exponential constant
        parameters(1,2)^3/3)*parameters(2,1);   
    diff(2) = 1;                                     % stress exp. 
    diff(3) = 3;                                     % grain-size exp. 
    diff(4) = parameters(2,3);                       % activation energy 
    % GBS
    gbs(1)  = (10*parameters(1,3)^2*...              % pre-exponential 
        parameters(1,2)^2/3)*parameters(2,1);   
    gbs(2) = 2;                                      % stress exp. 
    gbs(3) = 2;                                      % grain-size exp. 
    gbs(4) = parameters(2,3);                        % activation energy
    % monophase models with uniform phase grain-size
    ddc   = gs.^-dc(3);
    ddiff = gs.^-diff(3);
    dgbs  = gs.^-gbs(3);
    % compute T-dependent exponential term 
    exp_dc = exp(-(dc(4)*1000)./(R.*Taxis));   
    exp_diff = exp(-(diff(4)*1000)./(R.*Taxis));   
    exp_gbs = exp(-(gbs(4)*1000)./(R.*Taxis));                    
    for w = 1:nXY  
        % compute stress for each deformation mechanism
        Sii_dc   = (Eii./((dc(1)*ddc.*exp_dc(w)))).^(1/dc(2));
        Sii_diff = (Eii./((diff(1).*ddiff.*exp_diff(w)))).^(1/diff(2));
        Sii_gbs  = (Eii./((gbs(1).*dgbs.*exp_gbs(w)))).^(1/gbs(2));
        % find the mechanism that minimizes stress
        [Siimin(w),mech(w)] = min([Sii_dc,Sii_diff,Sii_gbs]);
    end
else
    % composite models
    molp = varargin{1};
    cvol = varargin{2};
    mixr = varargin{3};
    ivalid = 0;
    num = numel(phase_name) - sum(strcmp(phase_name,'not defined')); % get the number of valid phases
    mech = zeros(num,nXY);
    flowp = cell(1,num);
    valid_phase = strcmp(phase_name,'not defined');
    GS = gs(valid_phase~=1);
    
    % compute composite creep parameters (for each deformation mechanism) from selected mixing rule
    for i = 1:numel(phase_name)
        if strcmp(phase_name(i),'not defined')==0
            ivalid = ivalid+1;
            % loop through valid phases
            mat = load(['DB_mineral_parameters\',char(phase_name(i)),'.mat']);
            parameters  = mat.par;                       % extract material parameters
            % dislocation creep, array dc
            dc(1) = parameters(2,1);                     % pre-exponential constant 
            dc(2) = 4;                                   % stress exp. 
            dc(3) = 0;                                   % grain-size exp. 
            dc(4) = parameters(2,2);                     % activation energy 
            % diffusion creep
            diff(1) = (82*parameters(1,3)^3*...          % pre-exponential constant
                parameters(1,2)^3/3)*parameters(2,1);   
            diff(2) = 1;                                 % stress exp. 
            diff(3) = 3;                                 % grain-size exp. 
            diff(4) = parameters(2,3);                   % activation energy 
            % GBS
            gbs(1)  = (10*parameters(1,3)^2*...          % pre-exponential 
                parameters(1,2)^2/3)*parameters(2,1);   
            gbs(2) = 2;                                  % stress exp. 
            gbs(3) = 2;                                  % grain-size exp. 
            gbs(4) = parameters(2,3);                    % activation energy
            flowp(ivalid) = {[dc',diff',gbs']};          % parameters by column (1st column is dc, 2nd diff and 3rd gbs), each cell belongs to a mineral phase
        end
    end
    
    % evaluate stable deformation mechanism at various temperatures phase
    % by phase
    mech_phase = zeros(num,nXY);
    
    for i = 1:num
        % loop through phase to get the stable deformation mechanism at
        % each temperature as defined by the nXY axis
        for w = 1:nXY
            % loop through T
            % compute exponential term 
            exp_dc = exp(-(flowp{i}(4,1)*1000)./(R.*Taxis(w)));   
            exp_diff = exp(-(flowp{i}(4,2)*1000)./(R.*Taxis(w)));   
            exp_gbs = exp(-(flowp{i}(4,3)*1000)./(R.*Taxis(w))); 
            % evaluate stable deformation mechanism (for all T) of current (i) phase in composite
            Sii_dc   = (Eii./((flowp{i}(1,1)*exp_dc))).^(1/4);
            Sii_diff = (Eii./((flowp{i}(1,2).*GS(i)^-flowp{i}(3,2).*exp_diff))).^(1/1);
            Sii_gbs  = (Eii./((flowp{i}(1,3).*GS(i)^-flowp{i}(3,3).*exp_gbs))).^(1/2);
            % find the mechanism that minimizes stress
            [~,mech_phase(i,w)] = min([Sii_dc,Sii_diff,Sii_gbs]);
        end
    end
    % calculate flow law parameters of the composite at various T steps as defined by 
    % the stable deformation mechanism specified by mech_phase, and T axis. 
    % Different mixing rules defined by 'mixr' input argument (varargin{3})
    
    phases = 1:num;   % local indices of phases in composite
    
    for i = 1:nXY
        % loop through temperatures
        switch mixr
            case 'hobbs'
                % thermodynamic mixing rule
                A = 1;
                Q = 0;
                D = 1;
                numer = 1;
                den = 0;
                for j = 1:num
                    numer = numer*(flowp{j}(2,mech_phase(j,i)));
                    if num>1
                        den = den+molp(j)*prod(flowp{j}(2,mech_phase(phases~=j,i)));
                    else
                        den = 1;
                    end
                end
                n = numer/den;
                for j = 1:num
                    A = A*flowp{j}(1,mech_phase(j,i))^(molp(j)*n/flowp{j}(2,mech_phase(j,i)));
                    Q = Q+molp(j)*flowp{j}(4,mech_phase(j,i))/flowp{j}(2,mech_phase(j,i));
                    D = D*GS(j)^(-molp(j)*flowp{j}(3,mech_phase(j,i))*n/flowp{j}(2,mech_phase(j,i)));
                end
                Q = Q*n;
            case 'huet'
                % Minimized Power Geometric Mean 
                % loop through deformation mechanisms
                dmean = mean(GS);
                ai = zeros(1,num);        % ai = Prod(j~i){nj+1}
                aj = 1;                   % aj = Prod(i){nj+1}
                a1 = 1;
                a2 = 0;
                a3 = 1;
                n1 = 0;
                n2 = 0;
                n3 = 0;
                den1 = 0;
                den2 = 0;
                for j = 1:num
                    aj = aj*(flowp{j}(2,mech_phase(j,i))+1);
                    kk = phases(phases~=j); % get all phase indices but the current one
                    AI = 1;
                    for x = 1:num-1
                        AI = AI*prod(flowp{j}(2,mech_phase(kk(x),i))+1);
                    end
                    ai(j) = AI;
                    n1 = n1+cvol(j)*AI*flowp{j}(2,mech_phase(j,i));
                    n2 = n2+cvol(j)*AI*flowp{j}(4,mech_phase(j,i));
                    n3 = n3+cvol(j)*AI*flowp{j}(3,mech_phase(j,i));
                    den1 = den1+cvol(j)*AI;
                    den2 = den2+cvol(j)*aj;
                end
                n = n1/den1;                % mean stress exponent
                Q = n2/den1;                % mean activation energy of composite
                m = n3/den1;
                D = dmean^-m;
                for j = 1:num                  
                    a1 = a1*flowp{j}(1,mech_phase(j,i))^(ai(j)*cvol(j)/den2);
                    a2 = (a2+(cvol(j)*flowp{j}(2,mech_phase(j,i)))/(flowp{j}(2,mech_phase(j,i))+1));
                    a3exp = cvol(j)*flowp{j}(2,mech_phase(j,i))*ai(j)/den2;
                    a3 = a3*(flowp{j}(2,mech_phase(j,i))/(flowp{j}(2,mech_phase(j,i))+1))^a3exp;  
                end               
                A = a1*(a2^-n)*a3;
            case 'median'
                % median (central) value
                N = zeros(1,num);
                q = zeros(1,num);
                d = zeros(1,num);
                a = zeros(1,num);
                for j = 1:num
                    N(j) = flowp{j}(2,mech_phase(j,i));
                    q(j) = flowp{j}(4,mech_phase(j,i));
                    d(j) = GS(j)^-flowp{j}(3,mech_phase(j,i));
                    a(j) = flowp{j}(1,mech_phase(j,i));
                end
                n = median(N);
                Q = median(q);
                D = median(d);
                A = median(a);
            case 'mean'
                % arithmetic mean mixing rule
                N = zeros(1,num);
                q = zeros(1,num);
                d = zeros(1,num);
                a = zeros(1,num);
                for j = 1:num
                    N(j) = flowp{j}(2,mech_phase(j,i));
                    q(j) = flowp{j}(4,mech_phase(j,i));
                    d(j) = GS(j)^-flowp{j}(3,mech_phase(j,i));
                    a(j) = flowp{j}(1,mech_phase(j,i));
                end
                n = mean(N);
                Q = mean(q);
                D = mean(d);
                A = mean(a);
            case 'harmonic'
                % armonic mean mixing rule
                N = zeros(1,num);
                q = zeros(1,num);
                d = zeros(1,num);
                a = zeros(1,num);
                for j = 1:num
                    N(j) = flowp{j}(2,mech_phase(j,i));
                    q(j) = flowp{j}(4,mech_phase(j,i));
                    d(j) = GS(j)^-flowp{j}(3,mech_phase(j,i));
                    a(j) = flowp{j}(1,mech_phase(j,i));
                end
                n = num/sum(1./N);
                A = num/sum(1./a);
                Q = num/sum(1./q);
                D = num/sum(1./d);
        end
        % compute maximum stress at various T for the composite (NB: the active
        % grain-scale deformation mechanisms minimize stress at the local
        % scale)
        Siimin(i) = (Eii/(A*exp(-(Q*1000)/(R*Taxis(i)))*D))^(1/n);
        % compute dominant deformation mechanism (based on molar
        % proportions between phases
        mech(i) = mech_phase(phases(molp==max(molp)),i);      
    end
end
end
