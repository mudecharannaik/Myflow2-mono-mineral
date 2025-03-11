function[eff_mech_xy,Sii_xy] = eff_mech_Sii_xy(grain_D,grain_comp,listDB,Txy,meanEii)
% the function eff_mech_SiiEii_xy.m calculates the effective deformation
% mechanism at nodal points of a 2D indexed image. The image should
% represents either a monophase rock or composite made of up to 5
% different phases
% returns:
% eff_mech_xy  = 1-by-N array of index to active deformation mechanism, where 
%                N corresponds to the number of grains detected in the indexed 
%                image. The code of deformation
%                mechanism is specified this way: 1. dislocation creep, 
%                2. diffusion creep, 3. GBS. 
% Sii_xy = 1-by-N array of grain stress values 
%--------------------------------------------------------------------------
Txy = Txy+273.15;
R = 8.314;                                       % gas constant
N = numel(grain_D);                              % total number of grains in phasemap
indf = zeros(1,N);                               % LOCAL index (1 to n) of phases used to extract flow parameters
eff_mech_xy = nan(1,N);                          % initialize array of active deformation mechanisms in grains
Sii_xy = nan(1,N);                               % initialize array of stress value in grains
phaseID = unique(grain_comp);                    % ID of existent phases in phasemap;
n = numel(phaseID);                              % total number of phases in phasemap (max=5)
for i = 1:n
    indf(grain_comp==phaseID(i))=i;
end

% evaluate stable deformation mechanism as the one that locally minimizes
% stress (VR==1, Voigt) 
flowp = cell(1,n);
for i = 1:n                                      % loop through phases
    % first get phases flow parameters (based on composition)
    FP = zeros(4,3);
    mat = load(['DB_mineral_parameters\',char(listDB(phaseID(i))),'.mat']);
    parameters  = mat.par;                       % extract mineral parameters
    % dislocation creep, array dc
    FP(1,1) = parameters(2,1);                   % pre-exponential constant 
    FP(2,1) = 4;                                 % stress exp. 
    FP(3,1) = 0;                                 % grain-size exp. 
    FP(4,1) = parameters(2,2);                   % activation energy 
    % diffusion creep
    FP(1,2) = (82*parameters(1,3)^3*...          % pre-exponential constant
        parameters(1,2)^3)/3*parameters(2,1);   
    FP(2,2) = 1;                                 % stress exp. 
    FP(3,2) = 3;                                 % grain-size exp. 
    FP(4,2) = parameters(2,3);                   % activation energy 
    % GBS
    FP(1,3)  = (10*parameters(1,3)^2*...         % pre-exponential 
        parameters(1,2)^2)/3*parameters(2,1);   
    FP(2,3) = 2;                                 % stress exp. 
    FP(3,3) = 2;                                 % grain-size exp. 
    FP(4,3) = parameters(2,3);                   % activation energy
    flowp(i) = {FP};
end
for i = 1:N
    % loop through grains
    if ~isnan(flowp{indf(i)}(1,1))
        % loop through grains of indexed material
        % compute grainsize-dependent term 
        ddc = 1;
        ddiff = grain_D(i)^-flowp{indf(i)}(3,2);
        dgbs  = grain_D(i)^-flowp{indf(i)}(3,3);  
        % compute T-dependent exponential term 
        % Q is expressed in [J*mol-1], R = [J*K-1*mol-1], T = K so that
        % exp(-Q/RT) is balanced (dimensionless)
        exp_dc = exp(-(flowp{indf(i)}(4,1)*1000)/(R*Txy));    
        exp_diff = exp(-(flowp{indf(i)}(4,2)*1000)/(R*Txy)); 
        exp_gbs = exp(-(flowp{indf(i)}(4,3)*1000)/(R*Txy));    
        % Voigt (isostrain-rate)                 
        % compute stress for each deformation mechanism
        Sii_dc     = (meanEii./((flowp{indf(i)}(1,1)*ddc.*exp_dc))).^(1/flowp{indf(i)}(2,1));
        Sii_diff   = (meanEii./((flowp{indf(i)}(1,2)*ddiff.*exp_diff))).^(1/flowp{indf(i)}(2,2));
        Sii_gbs    = (meanEii./((flowp{indf(i)}(1,3)*dgbs.*exp_gbs))).^(1/flowp{indf(i)}(2,3));
        % find the mechanism that minimizes stress
        [val,pind] = min([Sii_dc,Sii_diff,Sii_gbs]);
        % save index of effective deformation mechanism in current grain
        eff_mech_xy(i) = pind; % save index of active deformation mechanism
        Sii_xy(i) = val; % save effective stress in MPa
    end
end
end