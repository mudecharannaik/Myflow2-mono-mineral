function[mech,Siimin] = Sii_d_voigt(phase_name,Daxis,nXY,T,Eii)
% the function Sii_d_voigt.m calculates the effective deformation
% mechanism at nodal points of the 2D grid nXYZ-by-nXYZ.
% returns:
% mech   = 1-by-nXYZ vector of index to active 
%          deformation mechanism, specified this way: 1. dislocation creep, 
%          2. diffusion creep, 3. GBS. 
% Siimin = 1-by-nXYZ vector of stress values for stable deformation
%          mechanism (at any given point)
%--------------------------------------------------------------------------
T = T+273.15;
R = 8.314;                                       % gas constant
mech = zeros(1,nXY);                             % initialize arrays 
Siimin = zeros(1,nXY);
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

% compute grainsize-dependent term
ddc   = Daxis.^-dc(3);
ddiff = Daxis.^-diff(3);
dgbs  = Daxis.^-gbs(3);
% compute T-dependent exponential term 
exp_dc = exp(-(dc(4)*1000)./(R*T));   
exp_diff = exp(-(diff(4)*1000)./(R*T));   
exp_gbs = exp(-(gbs(4)*1000)./(R*T));                    
for w = 1:nXY    
    % compute stress for each deformation mechanism
    Sii_dc   = (Eii./((dc(1)*ddc(w).*exp_dc))).^(1/dc(2));
    Sii_diff = (Eii./((diff(1).*ddiff(w).*exp_diff))).^(1/diff(2));
    Sii_gbs  = (Eii./((gbs(1).*dgbs(w).*exp_gbs))).^(1/gbs(2));
    % find the mechanism that minimizes stress
    [Siimin(w),mech(w)] = min([Sii_dc,Sii_diff,Sii_gbs]);
end
end