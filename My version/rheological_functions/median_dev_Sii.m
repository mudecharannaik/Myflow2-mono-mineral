function[Sii_dev,Sii_median] = median_dev_Sii(Sii_xymap,grain_D,grain_comp,listDB)
% the function eff_mech_SiiEii_xy.m calculates the effective deformation
% mechanism at nodal points of a 2D indexed image. The image should
% represents either a monophase rock or composite made of up to 5
% different phases
% returns:
% Sii_dev  = M-by-N array of percent deviations from a central stress value
%            calculated as median stress averaged by molar volume of
%            the contributing phases
% Sii_median = median stress value calculated based on molar proportions
%              bentween grains
%--------------------------------------------------------------------------
N = numel(grain_D);                              % total number of grains in phasemap
ip = nan(1,N);                                   % LOCAL index (1 to n_of_phases ~nan) of phases used to extract flow parameters
phaseID = unique(grain_comp);                    % ID of existent phases in phasemap;
n = numel(phaseID);                              % total number of phases in phasemap (max=5)
Vm = nan(1,n);                                   % molar volume of each phase
Sp = nan(1,n);                                   % total surface (equivalent 2D volume) of each phase
local_ID = 0;
for i = 1:n
    % get molar volume of ith phase
    mat = load(['DB_mineral_parameters\',char(listDB(phaseID(i))),'.mat']);
    parameters  = mat.par; 
    if ~isnan(parameters(1,1))
        local_ID = local_ID+1;
        % change global index of phase (from DB_mineral_parameters) to local
        % index ranging from 1 to n_local!
        ip(grain_comp==phaseID(i)) = local_ID;
        Vm(i) = parameters(1,1);
        % get partial surface of each phase (sum = total surface area)
        Sp(i)  = sum(pi.*(grain_D(ip==local_ID)./2).^2);
    end
end
% compute denominator in the expression of molar volumetric proportions
Sp = Sp(~isnan(Sp));                   % remove nan
Vm = Vm(~isnan(Vm));                   % remove nan
den = 0;
for i = 1:numel(Sp)
    den = den+Sp(i)*prod(Vm(Vm~=Vm(i)));
end
% compute molar volumetric proportions of n phases
alpha = nan(1,numel(Sp));                              
for i = 1:numel(Sp)
    alpha(i) = Sp(i)*prod(Vm(Vm~=Vm(i)))/den;
end
S = sum(Sp);          % total surface

% calculate mean grain stress (mgs) based on the molar volume proportion of each phase
% and the surface of each grain
mgs = zeros(1,N);
for i = 1:N
    if strcmp(listDB(grain_comp(i)),'not defined') ~= 1
        % calculate grain molar proportion as molar volumetric proportion*area_grain/total_area
        mgs(i) = alpha(ip(i))*(pi*(grain_D(i)/2)^2)/S;   
    end
end

% compute grain stresses averaged on molar proportions
Siim = Sii_xymap.*mgs;
% compute 'central' stress
Sii_median = sum(Siim,'omitnan');
% calculate deviations from central value
Sii_dev = (Sii_xymap-Sii_median)/Sii_median*100;
end