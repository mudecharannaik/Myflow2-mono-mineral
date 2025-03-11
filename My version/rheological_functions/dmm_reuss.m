function[mech] = dmm_reuss(flowp,nXY,T,Daxis,Siiaxis)
% the function eff_mech_dmm.m calculates the effective deformation
% mechanism at nodal points of the 2D grid nXYZ-by-nXYZ.
% returns:
% mech = nXYZ-by-nXYZ matrix of index to active 
%        deformation mechanism, specified this way: 1. dislocation creep, 
%        2. diffusion creep, 3. GBS. Only for single-phase materials
%--------------------------------------------------------------------------
T = T+273.15;
R = 8.314;                                       % gas constant
mech = zeros(nXY,nXY);                           % initialize arrays 

% compute T-dependent exponential term 
exp_dc = exp(-(flowp{1}(4)*1000)./(R.*T));   
exp_diff = exp(-(flowp{2}(4)*1000)./(R.*T));   
exp_gbs = exp(-(flowp{3}(4)*1000)./(R.*T));
% compute grainsize-dependent term
ddc   = Daxis.^-flowp{1}(3);
ddiff = Daxis.^-flowp{2}(3);
dgbs  = Daxis.^-flowp{3}(3);
for k = 1:nXY                        
    for w = 1:nXY    
        % compute strain-rate for each deformation mechanism
        Eii_dc   = flowp{1}(1).*(Siiaxis(w).^flowp{1}(2)).*ddc(k).*exp_dc;
        Eii_diff = flowp{2}(1).*(Siiaxis(w).^flowp{2}(2)).*ddiff(k).*exp_diff;
        Eii_gbs  = flowp{3}(1).*(Siiaxis(w).^flowp{3}(2)).*dgbs(k).*exp_gbs;
        % find the mechanism that minimizes strain-rate
        [~,pind] = min([Eii_dc,Eii_diff,Eii_gbs]);
        % save index of effective deformation mechanism
        mech(k,w) = pind;             
    end
end
end