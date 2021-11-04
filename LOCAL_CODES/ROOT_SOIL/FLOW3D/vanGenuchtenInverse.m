%% inverse vG that calculates H from theta.
function [H] = vanGenuchtenInverse(theta, alpha, theta_S, theta_R, n, m) 
%=========================================================================
% This function uses soil moisture and vanGenutchen parameters to calculate
% matric head of the soil
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       theta         % [-] soil moisture
%       alpha         % [1/cm] inverse air entry suction parameter
%       theta_S       % [-] Saturated water content
%       theta_R       % [-] Residual water content
%       n             % [-] Pore-size distributions
%       m             % (n-1)/n
%------------------------- Output Variables ------------------------------
%       H             % [m] soil matric potential
%-------------------------------------------------------------------------  

    H = -(1./alpha) .* (((theta_S-theta_R)./(theta-theta_R)).^(1./m)-1.0).^(1./n) * 0.01; % [cm] to [m]
    H(theta >= theta_S) = 0;
end

