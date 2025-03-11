function[mech] = dmm_voigt(flowp,nXY,T,Daxis,Eiiaxis)
% the function eff_mech_dmm.m calculates the effective deformation
% mechanism at nodal points of the 2D grid nXYZ-by-nXYZ.
% returns:
% mech = nXYZ-by-nXYZ matrix of index to active 
%        deformation mechanism, specified this way: 1. dislocation creep, 
%        2. diffusion creep, 3. GBS. Only for single-phase materials
%--------------------------------------------------------------------------
T = T+273.15;
R = 8.314;                                       % gas constant
mech = zeros(nXY,nXY);                         % initialize arrays 
% compute T-dependent exponential term 
exp_dc = exp(-(flowp{1}(4)*1000)./(R.*T));   
exp_diff = exp(-(flowp{2}(4)*1000)./(R.*T));   
exp_gbs = exp(-(flowp{3}(4)*1000)./(R.*T));
% compute grainsize-dependent term
ddc   = Daxis.^-flowp{1}(3);
ddiff = Daxis.^-flowp{2}(3);
dgbs  = Daxis.^-flowp{3}(3);
for k = 1:nXY       % grainsize loop                     
    for w = 1:nXY   % strainrate loop
        % compute strain-rate for each deformation mechanism
        Sii_dc   = (Eiiaxis(w)./((flowp{1}(1)*ddc(k).*exp_dc))).^(1/flowp{1}(2));
        Sii_diff = (Eiiaxis(w)./((flowp{2}(1)*ddiff(k).*exp_diff))).^(1/flowp{2}(2));
        Sii_gbs  = (Eiiaxis(w)./((flowp{3}(1)*dgbs(k).*exp_gbs))).^(1/flowp{3}(2));
        % find the mechanism that minimizes stress
        [~,pind] = min([Sii_dc,Sii_diff,Sii_gbs]);
        % save index of effective deformation mechanism
        mech(k,w) = pind;             
    end
end
end