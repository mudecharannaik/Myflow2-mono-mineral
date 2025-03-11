function[eff_mech_xy,Eii_xy] = eff_mech_Eii_xy(grain_D,grain_comp,listDB,Txy,meanSii)
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
% flowp_xy = 1-by-3 cell array containing the flow parameters of dislocation
%            creep (flowp{1}), diffusion creep (flowp{2}) and GBS
%            (flowp{3}), of all phases present within the indexed
%            image
%--------------------------------------------------------------------------
Txy = Txy+273.15;
R = 8.314;                                       % gas constant
N = numel(grain_D);                              % total number of grains in phasemap
indf = zeros(1,N);                               % LOCAL index (1 to n) of phases used to extract flow parameters
eff_mech_xy = nan(1,N);                          % initialize array of active deformation mechanisms in grains
Eii_xy = nan(1,N);                               % initialize array of strainrate value in grains
phaseID = unique(grain_comp);                    % ID of existent phases in phasemap;
n = numel(phaseID);                              % total number of phases in phasemap 
for i = 1:n
    indf(grain_comp==phaseID(i))=i;
end

% evaluate stable deformation mechanism as the one that locally minimizes
% strainrate (VR==0, Reuss)
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
        parameters(1,2)^3/3)*parameters(2,1);   
    FP(2,2) = 1;                                 % stress exp. 
    FP(3,2) = 3;                                 % grain-size exp. 
    FP(4,2) = parameters(2,3);                   % activation energy 
    % GBS
    FP(1,3)  = (10*parameters(1,3)^2*...         % pre-exponential 
        parameters(1,2)^2/3)*parameters(2,1);   
    FP(2,3) = 2;                                 % stress exp. 
    FP(3,3) = 2;                                 % grain-size exp. 
    FP(4,3) = parameters(2,3);                   % activation energy
    flowp(i) = {FP};
end

for i = 1:N                                      % loop through grains
    if ~isnan(flowp{indf(i)}(1,1))
        % compute grainsize-dependent term 
        ddc = grain_D(i)^-flowp{indf(i)}(3,1);
        ddiff = grain_D(i)^-flowp{indf(i)}(3,2);
        dgbs  = grain_D(i)^-flowp{indf(i)}(3,3);  
        % compute T-dependent exponential term 
        exp_dc = exp(-(flowp{indf(i)}(4,1)*1000)/(R*Txy));
        exp_diff = exp(-(flowp{indf(i)}(4,2)*1000)/(R*Txy));
        exp_gbs = exp(-(flowp{indf(i)}(4,3)*1000)/(R*Txy));    
        % Reuss (isostress)     
        % compute strain-rate for each deformation mechanism
        Eii_dc   = flowp{indf(i)}(1,1).*meanSii.^flowp{indf(i)}(2,1).*ddc.*exp_dc;
        Eii_diff = flowp{indf(i)}(1,2).*meanSii.^flowp{indf(i)}(2,2).*ddiff.*exp_diff;
        Eii_gbs  = flowp{indf(i)}(1,3).*meanSii.^flowp{indf(i)}(2,3).*dgbs.*exp_gbs;
        % find the mechanism that minimizes strain-rate
        [val,pind] = min([Eii_dc,Eii_diff,Eii_gbs]);
        % save index of effective deformation mechanism in current grain
        eff_mech_xy(i) = pind; 
        Eii_xy(i) = log10(val);
    end
end
end