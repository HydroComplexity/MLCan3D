function [znc, zhc, dzc, zlai] = Canopy_Grid(PARAMS)
%
%   Construct the canopy grid by evenly spacing 'nl_can' nodes between the
%   soil surface and canopy top at 'hcan'
%
%   INPUTS: 
%       nl_can = # of canopy layers
%       hcan = height of the canopy         [m]
%       LAInorm = normalized LAI profile
%
%   OUTPUTS:
%       znc = canopy node heights           [m]
%       zhc = canopy interface heights      [m]
%
%   Written By: Darren Drewry
%

% De-reference Structure Values
    nl_can = PARAMS.CanStruc.nl_can;
    hcan = PARAMS.CanStruc.hcan;

% CALCULATE NODE AND LAYER SPACING
    zhc = [1:nl_can] * (hcan / nl_can);
    
    dzc = zhc(1) / 2;
    
    znc = zhc - dzc;
    
% CALCULATE NORMALIZED LEAF AREA PROFILE
    zb = find(znc>(floor(hcan/4)), 1, 'first');
    zlai = zeros(size(znc));
    zlai(zb:end) = 1;
    zlai = zlai / sum(zlai);
