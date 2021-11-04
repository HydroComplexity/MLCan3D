function [ Conc_can ] = SIMILARITY_ADJ(FORCING, PARAMS, CONSTANTS, H, E, Conc_top, case_var)

%=========================================================================
% Use similarity to calculate water, heat, windspeed at the top of the
% canopy using water, heat, windspeed at the tower height
%
% Written By: Esther Lee
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING         % FORCING structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%------------------------- Output Variables ------------------------------
%       U               % [m/s] Wind speed 
%       km              % Km momentum diffusivity 
%       tau             % parameter momentum model
%       Conc_can        % [vary] Concentration of tracer considered 
%                         at the top of the canopy
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    Utop = FORCING.U;                                                      %[m/s] Wind Speed at the top (Forcing from atmospheric measurements)  
    Ustar = FORCING.Ustar;                                                 %[m/s] Friction Velocity
    hcan = PARAMS.CanStruc.hcan;                                           %[m] Height canopy  
    hobs = PARAMS.CanStruc.hobs;                                           %[m] Height flux tower %added by Esther to use for Utop calculation
    vonk = CONSTANTS.vonk;                                                 % [] von Karman constant 
    %case_var;                                                              %1. for U  2. for T (assume it is the same as potential T)   3. for specific humidity
    z0 = PARAMS.Soil.z0;                                                   % [m] surface roughness length
    rho = CONSTANTS.rho_dry_air;                                           %density of dry air [kg/m3] check if is is the right rho!!!
    cp = CONSTANTS.cp_JkgK;                                                   %specific heat of air at constant pressure [J/kg/K]

%*************************************************************************


d0 = hcan*2/3;                                                                %Esther: zero plane displacement
% Ustar = Utop*vonk/(log((hobs-d0)/z0));                                        %Esther: friction velocity



if case_var==1                  %case 1. for U      added by Esther
    Conc_can = Conc_top + Ustar/vonk * log((hcan-d0) / (hobs-d0));
    
elseif case_var==2              %case 2. for T (assume it is the same as potential T)       added by Esther  
%     H=H/60/30;
    if Ustar<.01 && abs(H)>10
        Ustar;
    end
    Conc_can = Conc_top - H/vonk/Ustar/rho/cp * log((hcan-d0) / (hobs-d0));
    if abs(Conc_can - Conc_top) > 15
%         Conc_can = Conc_top;
    end
elseif case_var==3              %case 3. for specific humidity      added by Esther
    Conc_can = Conc_top - E/vonk/Ustar/rho * log((hcan-d0) / (hobs-d0));
end
