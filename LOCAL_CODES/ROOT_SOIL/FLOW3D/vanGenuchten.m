function [C,K,theta] = vanGenuchten(h, alpha, theta_S, theta_R, n, m, Ksat)    
%=========================================================================
% This function uses vanGenuchten parameters to calculate soil moisture and
% conductivity from matric pressure
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       h             % [m] soil matric potential
%       alpha         % [1/cm] inverse air entry suction parameter
%       theta_S       % [-] Saturated water content
%       theta_R       % [-] Residual water content
%       n             % [-] Pore-size distributions
%       m             % (n-1)/n
%       Ksat          % [m/hr] Saturated hydraulic conductivity 
%------------------------- Output Variables ------------------------------
%       C             % [1/m] Specific moisture storage
%       K             % [m/hr] hydraulic conductivity given 
%       theta         % [-] soil moisture
%-------------------------------------------------------------------------  
    % Convert unit of h from [m] to [cm] to match with alpha
    h = h*100;

    % Compute the volumetric moisture content [eqn 21]
    theta = (theta_S - theta_R)./(1 + (alpha.*abs(h)).^n).^m + theta_R;
    if (length(theta_S) > 1)  % input is array
        theta(h>=0.0) = theta_S(h>=0.0);
    else
        theta(h>=0.0) = theta_S;    % input is single value
    end
    
    % Compute the effective saturation [eqn 2]
    Se = ((theta - theta_R)./(theta_S - theta_R));
 
    % Compute the hydraulic conductivity [eqn 8]
    K = Ksat.*Se.^(1/2).*(1 - (1 - Se.^(1./m)).^m).^2;
 
    % Compute the specific moisture storage (derivative of eqn 21: 
    % C = d(theta)/dh
    C = -alpha.*n.*sign(h).*(1./n - 1).*(alpha.*abs(h)).^(n - 1).*(theta_R - theta_S).*((alpha.*abs(h)).^n + 1).^(1./n - 2)*100;
    C(h>=0.0) = 0;
end