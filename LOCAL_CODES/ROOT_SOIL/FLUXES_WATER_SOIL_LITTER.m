
function [VARIABLES] = FLUXES_WATER_SOIL_LITTER (PARAMS, VARIABLES, CONSTANTS,...
    FORCING, SWITCHES)

%=========================================================================
% Solve surface energy balance with a snow-litter pack
%
% Written by Juan Quijano, UIUC, 2013
% Modified by Kunxuan Wang, UIUC, 2020
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       PARAMS          % PARAMS structure
%       VARIABLES       % VARIABLES structure
%       CONSTANTS       % CONSTANTS structure
%       FORCING         % CONSTANTS structure
%       SWITCHES        % SWITCHES structure
%------------------------- Output Variables ------------------------------
%       VARIABLES       % VARIABLES structure
% 
%========================================================================
    
    % De reference block
    %*********************************************************************
    % PARAMS
    thetals = PARAMS.Soil.thetals;                                          % Litter Soil Moisture at saturation
    thetafc = PARAMS.Soil.thetafc;                                          % Litter moisture at field capacity
    km = PARAMS.Soil.km;                                                    % parameter to compute drainage from litter
    bm = PARAMS.Soil.bm;                                                    % parameter to compute drainage from litter
    Tcr_atm = PARAMS.Tcr_atm;                                               % [K]Critical Atmospheric Temperature around 273 K
    Tf = PARAMS.Soil.Tf;                                                    % [K] Freezing Temperature of Fresh Water
    rho_liq = PARAMS.Soil.rho_liq;                                          % [kg / m^3] Density of Liquid Water
    rho_ice = PARAMS.Soil.rho_ice;                                          % [km / m^2] Density of ice
    HC_snow = PARAMS.Soil.HC_snow;                                          % [] Maximum fraction water capacity that snow can hold
    rhosn_min = PARAMS.Soil.rhosn_min;                                      % [kg / m^3] Minimum Snow Density  
    fcomp = PARAMS.Soil.fcomp;                                              % Compactation parameter           
    slfactor = PARAMS.Soil.slfactor;                                        % Percentage above which is considered as snow for energy balance
            
    % VARIABLES
    wliqsl = VARIABLES.SOIL.wliqsl;                                         % [kg / m^2] Liquid Water Density per unit area
    wicesl = VARIABLES.SOIL.wicesl;                                         % [kg / m^2] Ice Water Density per unit area
    wsn = VARIABLES.SOIL.wsn;                                               % [kg / m^2] Total Snow Density per unit area
    zliqsl = VARIABLES.SOIL.zliqsl;                                         % [mm] Liquid water depth
    zicesl = VARIABLES.SOIL.zicesl;                                         % [mm] Ice water depth
    zsn = VARIABLES.SOIL.zsn;                                               % [mm] Snow water depth
    volliqli = VARIABLES.SOIL.volliqli;                                     % [] Volumetric Water Content Snow-Litter Pack
    voliceli = VARIABLES.SOIL.voliceli;                                     % [] Volumetric Ice Content Snow-Litter Pack
    deltawice = VARIABLES.SOIL.deltawice;                                   % [] Change in volumetric volume of ice
    rhosn = VARIABLES.SOIL.rhosn;                                           % [kg / m^3] Density of snow
    Tli = VARIABLES.SOIL.Tli;                                               % [C]Temperature of Snow-Litter pack
    ppt_ground = VARIABLES.SOIL.ppt_ground;                                 % [mm] precipitation reaching ground
    snow_tcount = VARIABLES.SOIL.snow_tcount;                               % [] Number of time steps to account snow age
    if SWITCHES.soil3D
        pond_depth = VARIABLES.OVERLAND2D.pond_depth;                           % [m] 2D matrix of ponding depth from 3D flow model
    else
        pond_depth = 0;
    end
    
    case1sl = VARIABLES.SOIL.case1sl;                                       % Index that defines whether is treated as snow or as a litter
    
    dzlit_mm = VARIABLES.SOIL.dzlit_m*1000;                                 % Thickness of litter [mm]

    grav = CONSTANTS.grav;                                                  % [m / s^2] Gravity Acceleration 
    dtime = CONSTANTS.dtime;                                                % [s] Time Step
    Lv_kg = CONSTANTS.Lv_kg;                                                % [J/kg] Latent Heat of Evaporation 
    Lf_kg = CONSTANTS.Lv_kg;                                                % [J/kg] Latent Heat of Fusion 
    Ls_kg = Lv_kg + Lf_kg;                                                  % [J/kg] Latent Heat of Sublimation
    
    LEsl = VARIABLES.SOIL.LEl;                                              % [W/m^2] Latent Heat from snow pack
    LEs = VARIABLES.SOIL.LEs;                                               % [W/m^2] Latent Heat from the soil

    Ta = FORCING.Ta;                                                        % [C] Atmospheric temperature
    %*********************************************************************
    %i. Compute LE fluxes 
    Esoil = LEs / Lv_kg;                                                   % Evaporation From Soil only [mm liq /s]
    if case1sl == 1        
        Esl_ice = LEsl / Ls_kg * (1000/rho_ice);                           %[mm ice / s]Sublimation From Sow-Litter Pack 
        Esl = Esl_ice;
        Esl_liq = 0;
    else
        Esl_ice = 0;                  
        Esl_liq = LEsl / (Lv_kg) * (1000/rho_liq);                         %[mm liq / s] Evaporation From Sow-Litter Pack 
        Esl = Esl_liq;
    end 
    
    % ii. snow density of snow reaching ground [kg/m3]
    if Ta < -15 
        rhosn_new = 50;
    elseif (Ta >= -15 && Ta <= 2)
        rhosn_new = 50 + 1.7*(Ta+15)^(1.5); 
    elseif Ta > 2
        rhosn_new = 50 + 1.7*(17)^(1.5);
    end    
    % iii. compute rainfall and ponded water
    pond_avg = mean(pond_depth, 'all')*1000;   %[mm]   
    pptrate_ground = (ppt_ground+pond_avg) / dtime;  % add ponded water to ppt to calc together
    if (Ta < (Tcr_atm - Tf) && pptrate_ground > 0)
        zice_gro = pptrate_ground*dtime *rho_liq/rho_ice;                  %[mm ice] 
        zliq_gro = 0;
        snow_tcount = 0;
    else        
        zice_gro = 0;
        zliq_gro = pptrate_ground*dtime;                                   %[mm liq] 
        rhosn_new = rho_liq;
        snow_tcount = snow_tcount + 1;
    end
        
    % iv. Computation of Drainage
    Clitter = volliqli*dzlit_mm;   % Water retained n the litter    
    if volliqli < thetafc
         drainlitter = 0;
    else
         drainlitter = km*exp(bm*Clitter)*dtime;                        % Water that drains into soil from litter [mm]            
    end   
    if case1sl == 1
       drainlitter = max((wliqsl/wsn - HC_snow)*wsn/(rho_liq)*1000,0);    % Water that drains into soil from litter [mm]
       %drainlitter = min(drainlitter_snow,drainlitter); 
    end
    
    % v. Compute infiltration into snow-litter pack and into soil    
    if SWITCHES.soilevap   
        qinfl = drainlitter / dtime - Esoil;                                %[mm/s] Infiltration into soil with extraction of soil evaporation 
    else
        qinfl = drainlitter / dtime;                                        %[mm/s] Infiltration into soil with no extraction of soil evaporation
    end
    qinflL = pptrate_ground - Esl;                           % Water that infiltrates into litter    
    net_qinflL = pptrate_ground - Esl - drainlitter/dtime;   % Net Water that infiltrates into litter [mm] (with no evaporation no drainage)
    
    % vi. Compute mass balance
    %remains 
    remainsl = wliqsl + wicesl - (drainlitter/1000).*rho_liq ...
               - (Esl_liq/1000)*dtime*rho_liq - (Esl_ice/1000)*dtime*rho_ice; %[kg/m2]
    % liq and ice
    dwliq_mm = -deltawice/rho_liq*1000;                                           %[mm liq]
    zliqsl_new = zliqsl + (zliq_gro - drainlitter - Esl_liq*dtime + dwliq_mm);    %[mm liq]    
    if zliqsl_new < 0 
        qinfl = qinfl - (0 - zliqsl_new)/dtime;
        drainlitter = max(drainlitter - (0 - zliqsl_new),0);                        %[mm liq]
        zliqsl_new = 0;                                                             %[mm liq]
        % correct qinfl
    end    
%    zliqsl_new = max(0,zliqsl_new);                                               %[mm liq]
    wliqsl_new = (zliqsl_new/1000)*rho_liq;                                       %[kg H2O /m2]
     
    dwice_mm = deltawice/rho_ice*1000;                                            %[mm ice]
    zicesl_new = zicesl + (zice_gro - Esl_ice*dtime + dwice_mm);                  %[mm ice]    
    if (zicesl_new < 0 && SWITCHES.soilevap)
        qinfl = qinfl + zicesl_new/dtime;                                         %[mm ice]
    end    
    zicesl_new = max(0,zicesl_new);                                               %[mm ice]
    wicesl_new = (zicesl_new/1000)*rho_ice;                                       %[kg H2O /m2] Is Ice but Equivalent of water
    
    
    
    
    %vii. Compute snow density
    % ***** compute CRm ***** 
    if rhosn <= 150 
        c3 = 1;
    else
        c3 = exp(-0.046*(rhosn-150));
    end
    if zliqsl == 0
        c4 = 1;
    else
        c4 = 2;
    end
    CRm = 2.778*10^(-6)*c3*c4*exp(-0.04*(-Tli)); 
    % ***** compute CRo ***** 
    c5 = 0.08;
    c6 = 0.021;
    Ps = 0.5*grav*rho_liq*(zice_gro + zliq_gro + fcomp*zsn);
    CRo = (Ps/(3.6*10^6))*exp(-c5*(-Tli))*exp(-c6*rhosn);
    % ***** compute change *****
    drhosn = (CRm + CRo)*rhosn*dtime; 
    rhosn_updt = rhosn + drhosn;
    if (zicesl < 0.01*dzlit_mm*thetals || rhosn_updt > rho_liq)                         
        rhosn_updt = rho_liq;                               % rho_updt is equal to rho_liq if ice is very small
    end
    if rhosn_updt < rhosn_min
        rhosn_updt = rhosn_min;
    end    
    %compact_net = rhosn_updt/rhosn;
    % ponderate rho of water in snow-litter pack.    
    
    % The final snow density is computed by a mass wighted average of the
    % different densities involved
    rhosn = ((rhosn_updt*wsn) - ((Esl_liq/1000)*dtime*rho_liq)*rho_liq ... 
               - ((Esl_ice/1000)*dtime*rho_ice)*rho_ice - ((drainlitter/1000).*rho_liq)*rho_liq ...
               + ((zice_gro/1000)*rho_liq)*rhosn_new + ((zliq_gro/1000)*rho_liq)*rhosn_new) / ...
                 (remainsl + (zice_gro/1000)*rho_liq + (zliq_gro/1000)*rho_liq);
    rhosn = max(50,rhosn);                                  % Snow density can not be lower than 50

               
    % vii. Save new variables.
    % w
    VARIABLES.SOIL.wliqsl = wliqsl_new;                                     % [kg / m^2] Liquid Water Density per unit area  
    VARIABLES.SOIL.wicesl = wicesl_new;                                     % [kg / m^2] Ice Water Density per unit area 
    % z
    VARIABLES.SOIL.zicesl = zicesl_new;                                     % [mm] Ice water depth  
    VARIABLES.SOIL.zliqsl = zliqsl_new;                                     % [mm] Liquid water depth 
    VARIABLES.SOIL.zicesl_prev = zicesl;                                    % [mm] Ice water depth in previous time step
    VARIABLES.SOIL.zliqsl_prev = zliqsl;                                    % [mm] Liquid water depth in previous time step

    if zicesl_new < slfactor*dzlit_mm*thetals    
        zsn = 0;
        wsn = 0;
        rhosn = rho_liq;                                                    % rho is equal to rho_liq if rho snow is very small        
    else
        wsn = wicesl_new + wliqsl_new;                                      % [kg / m^2]
        zsn = (wsn/rhosn)*1000;                                             % [mm] Snow depth
    end
    % rhosn
    VARIABLES.SOIL.rhosn = rhosn;                                           % [kg / m^3] Density of snow
    
    VARIABLES.SOIL.wsn = wsn;                                               % [kg / m^2] Snow Water Density per unit area     
    VARIABLES.SOIL.zsn = zsn;                                               % [mm]
    % vols
    VARIABLES.SOIL.volliqli = zliqsl_new/dzlit_mm;                          % [] Volumetric Water Litter Content     
    VARIABLES.SOIL.voliceli = zicesl_new/dzlit_mm;                          % [] Volumtetric Ice Litter Content     

    VARIABLES.SOIL.volliqsn = zliqsl_new/zsn;                               % [] Volumetric Liquid Water Content in Snow      
    VARIABLES.SOIL.volicesn = zicesl_new/zsn;                               % [] Volumetric Ice Water Content in Snow         
    VARIABLES.SOIL.voltotsn = (zliqsl_new + zicesl_new)/zsn;                % [] Total Volumetric Ice and Water Snow Content         
    VARIABLES.SOIL.voltotli = (zliqsl_new + zicesl_new)/dzlit_mm;           % [] Total Volumetric Ice and Water Litter Content 
    
    %Other variables 
    VARIABLES.SOIL.CRm = CRm;                                               % Snow Compaction Factor
    VARIABLES.SOIL.CRo = CRo;                                               % Snow Compaction Factor
    VARIABLES.SOIL.drhosn = drhosn;                                         % [kg / m^3] Change in snow density 
    
    VARIABLES.SOIL.qinflL = qinflL;                                         %[mm/s] Water that infiltrates into snow
    VARIABLES.SOIL.net_qinflL = net_qinflL;                                 %[mm/s] Net water that infiltrates into snow excluding drainage in soil  
    VARIABLES.SOIL.qinfl = qinfl;                                           %[mm/s] Net Infiltration into the soil
    VARIABLES.SOIL.drainlitter = drainlitter;                               %[mm] Drainage from snow into soil
    VARIABLES.SOIL.pptrate_ground =  pptrate_ground;                        %[mm/s] Precipitation reaching the ground
    VARIABLES.SOIL.Esoil = Esoil;                                           %[mm/s] Evaporation soil 
    VARIABLES.SOIL.Esl = Esl;                                               %[mm/s] Sublimation from snow 
    VARIABLES.SOIL.zliq_gro = zliq_gro;                                     %[mm] Depth of liquid water reaching snow
    VARIABLES.SOIL.zice_gro = zice_gro;                                     %[mm] Depth of ice water reaching snow
    %
    VARIABLES.SOIL.Esl_ice = Esl_ice;                                       %[mm/s] Sublimation from ice in the snow
    VARIABLES.SOIL.Esl_liq = Esl_liq;                                       %[mm/s] Evaporation from snow (assumed zero)
    % Count time since last snow
    VARIABLES.SOIL.snow_tcount = snow_tcount;                               %[] Number of time steps with snow  
    