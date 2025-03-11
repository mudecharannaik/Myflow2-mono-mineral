function[flowp] = set_flowp(phase_name)
% the function extract the flow parameters from the mineral database based 
% on the phase specified by phase_name
% returns:
% flowp = 1-by-3 cell array containing the flow parameters of dislocation
%                 creep (flowp{1}), diffusion creep (flowp{2}) and GBS
%                 (flowp{3}), of the only phase in monophase materials
%--------------------------------------------------------------------------
dc = zeros(4,1);                                 % dislocation creep parameters
diff = zeros(4,1);                               % diffusion creep parameters
gbs = zeros(4,1);                                % GBS parameters
mat = load(['DB_mineral_parameters\',char(phase_name),'.mat']);
parameters  = mat.par;                       % extract material parameters
% dislocation creep, array dc
dc(1) = parameters(2,1);                   % pre-exponential constant 
dc(2) = 4;                                 % stress exp. 
dc(3) = 0;                                 % grain-size exp. 
dc(4) = parameters(2,2);                   % activation energy 
% diffusion creep
diff(1) = (82*parameters(1,3)^3*...        % pre-exponential constant
    parameters(1,2)^3/3)*parameters(2,1);   
diff(2) = 1;                               % stress exp. 
diff(3) = 3;                               % grain-size exp. 
diff(4) = parameters(2,3);                 % activation energy 
% GBS
gbs(1)  = (10*parameters(1,3)^2*...        % pre-exponential 
    parameters(1,2)^2/3)*parameters(2,1);   
gbs(2) = 2;                                % stress exp. 
gbs(3) = 2;                                % grain-size exp. 
gbs(4) = parameters(2,3);                  % activation energy 
flowp    = cell(1,3);
flowp(1) = {dc};
flowp(2) = {diff};
flowp(3) = {gbs};
end