function[MOLP] = molar_proportion(comp)
% the function calculates the molar proportion of each phase in composite
% materials
%--------------------------------------------------------------------------
phase_name = comp.phase_name;       % name of phases in composite
cvol = comp.phase_vol;              % volume proportions of phases
if sum(strcmp(phase_name,'not defined'))~=0
    % the composite contains some not defined mineral
    % recalculating volume proportions of mineral phases
    newvol = sum(cvol(strcmp(phase_name,'not defined')~=1));
    for i = 1:numel(phase_name)
        if strcmp(phase_name(i),'not defined')~=1
            cvol(i) = cvol(i)+(cvol(i)*(1.0-newvol))/1.0;
        end
    end
end 
cvol = nonzeros(cvol);
% initialize array of results
molp = zeros(1,numel(phase_name));
% array of phase names
ca_plg    = {'anorthite dry','anorthite wet','plagioclase'};
na_plg    = {'albite dry','albite wet'};
feldspar  = {'Kfeldspar'};
amph      = {'amphibole','hornblende'};
calcite   = {'calcite','aragonite'};
quartz    = {'quartz'};
coesite   = {'coesite'};
garnet    = {'garnet','almandine'};
pyroxene  = {'diopside','pyroxene'};
jadeite   = {'jadeite'};
muscovite = {'muscovite','phengite'};
biotite   = {'biotite','siderophyllite','phlogopite'};
olivine   = {'olivine','olivine dry','olivine wet','fayalite','forsterite'};
magnetite = {'magnetite'};

% calculate molar volumes o
for i = 1:numel(phase_name)
    if sum(strcmp(phase_name(i),ca_plg))>0 
        molp(i) = 10.07;
    elseif sum(strcmp(phase_name(i),na_plg))>0
        molp(i) = 10.02;
    elseif sum(strcmp(phase_name(i),feldspar))>0
        molp(i) = 10.86;
    elseif sum(strcmp(phase_name(i),amph))>0
        molp(i) = 27.85;
    elseif sum(strcmp(phase_name(i),calcite))>0
        molp(i) = 3.69;
    elseif sum(strcmp(phase_name(i),quartz))>0
        molp(i) = 2.27;
    elseif sum(strcmp(phase_name(i),garnet))>0
        molp(i) = 11.52;
    elseif sum(strcmp(phase_name(i),pyroxene))>0
        molp(i) = 6.61;
    elseif sum(strcmp(phase_name(i),jadeite))>0
        molp(i) = 6.09;
    elseif sum(strcmp(phase_name(i),muscovite))>0
        molp(i) = 14.05;
    elseif sum(strcmp(phase_name(i),biotite))>0
        molp(i) = 14.96;
    elseif sum(strcmp(phase_name(i),olivine))>0
        molp(i) = 62.299;
    elseif sum(strcmp(phase_name(i),magnetite))>0
        molp(i) = 1.2;
    elseif sum(strcmp(phase_name(i),coesite))>0
        molp(i) = 2.06;
    end   
end
molp = nonzeros(molp);
% calculate molar proportions
% as (i.e., three phase of volume v1, v2, and v3 and molar volumes
% vm1, vm2 and vm3): 
% alpha1 = v1*vm2*vm3/(v1*vm2*vm3+v2*vm1*vm3+v3*vm1*vm2)
MOLP = zeros(1,numel(molp));
den = 0;
local_ID = 1:numel(molp);
for i = 1:numel(molp)
    den = den+cvol(i)*prod(molp(local_ID~=i));
end
for i = 1:numel(molp)
    MOLP(i) = cvol(i)*prod(molp(local_ID~=i))/den;
end
end