function [An_can, Ph_can, LE_can, H_can, dHcan, Rnrad_can, TR_can, ...
          Fc_soil, LE_soil, H_soil, Rnrad_soil, G, Tsurf, remainsoil, remaincan,remaineco, ...
          Rnrad_sun, Rnrad_shade, Rnrad_eco, ...
          An_sun, An_shade, LE_sun, LE_shade, H_sun, H_shade, TR_sun, TR_shade, ...
          Tl_sun, Tl_shade, psil_sun, psil_shade, gsv_sun, gsv_shade, fsvg_sun,fsvm_sun, ...
          fsvg_shade,fsvm_shade,Ci_sun, Ci_shade, CAz, TAz, EAz, Uz, gbv_sun, gbh_sun, gbv_shade, gbh_shade, ...
          LAIsun, LAIshade, fsun, fshade, ...
          Ph_limit_sun, Jc_C3_sun, Jj_C3_sun, Js_C3_sun, Jc_C4_sun, Jj_C4_sun, Js_C4_sun, ...
          Ph_limit_shade, Jc_C3_shade, Jj_C3_shade, Js_C3_shade, Jc_C4_shade, Jj_C4_shade, Js_C4_shade, ...
          PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, SWout, ...
          LWabs_can, LWabs_soil, LWemit_soil, LWemit_can, LWemit_sun, LWemit_shade, LWout, LWoutM, RH_soil, fdiff, ...
          Sh2o_prof, Sh2o_can, ppt_ground, Ch2o_prof, Ch2o_can, Evap_prof, Evap_can, ...
          dryfrac, wetfrac, Vz, VARIABLES, FORCING, ...
          SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out, SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
          LWabs_canM, LWabs_soilM, LSshaCON, LSsunCON, countingmaxedout, diverging, largeremainder, outlier, redo_noturb] = ...    
        CANOPY_MODEL(SWITCHES, VERTSTRUC, FORCING, PARAMS, VARIABLES, CONSTANTS, tt, VARIABLES_last_tt, countingmaxedout, diverging, largeremainder, outlier, redo_noturb)

      % Last updated by Meredith Richardson, October 29, 2018
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % SWITCHES
    turb_on = SWITCHES.turb_on;                                            % Switch for simulation of turbulence in atmosphere [1] Yes [0] No     
    Sh2o_prof = VARIABLES.CANOPY.Sh2o_prof;                                % [mm] Intercepted water in canopy, for all layers     
    Tsoil = VARIABLES.SOIL.Ts(1);                                          % [C] Temperature in the Soil top layer 
    canstorheat = SWITCHES.canstorheat;                                    %[] Switch for considering canopy storage of heat [1] yes [0] no
    
    % Mere SWITCHES
    removing_turb_outliers = 1;                                            % Removing outliers in turbulence calculations; [1] yes [0] no
    replace_remainder_errors = 0;                                          % Replacing timesteps with large energy balance errors; [1] yes [0] no; Not recommended if LAI goes to zero during simulation period!
    similarity_adjustment = 1;                                             % Use similarity adjustment to scale concentration from tower to canopy top; [1] yes [0] no
    if (redo_noturb) 
        turb_on = 0;
        redo_noturb = 0;
    end
    
    % VERTSTRUC
    vinds = VERTSTRUC.vinds;                                               % [] Indicator of Layers where LAI is zero of lower  
    nvinds = VERTSTRUC.nvinds;                                             % [] Indicator of Layers where LAI is zero of lower for all species         
    LAIz = VERTSTRUC.LAIz;                                                 % [m^2 leaf/ m^2 ground] LAI for all layers 
    LAIzall = VERTSTRUC.LAIzall;                                           % [m^2 leaf/ m^2 ground] LAI for all layers and species  
    fLAIz = VERTSTRUC.fLAIz;                                               % [] Fraction of total LAI for each species at each canopy layer 
    
    % PARAMS
    Ro = PARAMS.Resp.Ro;                                                   % [umol / m^2 / s] Ecosystem Respiration Parameter Q10 model
    Q10 = PARAMS.Resp.Q10;                                                 % [] Ecosystem respiration parameter Q10 model 
    nl_can = PARAMS.CanStruc.nl_can;                                       % [] Number of layers in canopy 
    nspecies = PARAMS.CanStruc.nspecies;                                   % [] Number of species 
    
    % CONSTANTS
    dtime = CONSTANTS.dtime;                                               % [s] Time step
    
    % Variables (Dongkook: For the output variables defined here)
    CAz = VARIABLES.CANOPY.CAz;                                            % [umol/mol] Atmospheric Concentration of CO2, for all layers
    TAz = VARIABLES.CANOPY.TAz;                                            % [C] Air Temperature, for all layers 
    EAz = VARIABLES.CANOPY.EAz;                                            % [kPa] Vapor Pressure Air, for all layers 
    
    Utop = FORCING.U;                                                      %[m/s] Wind Speed at the top (Forcing from atmospheric measurements)  
    
    
%*************************************************************************
%                          ALLOCATE MATRICES

%    ALLOCATE_STORAGE_CAN();
%percdiffprof = logical(zeros(nl_can,1));
%*************************************************************************
%*************************************************************************
             
% Vertical distribution of photosynthetic capacity
    [Vz] = PH_Dist(VERTSTRUC, PARAMS);
        VARIABLES.CANOPY.Vz = Vz;

    [Uz, Km] = ORDER_1_CLOSURE_U(VERTSTRUC, PARAMS, Utop);  
%       uinds = find(Uz < 0.05);
%       Uz(uinds) = 0.05;
        VARIABLES.CANOPY.Uz = Uz;
        VARIABLES.CANOPY.Km = Km;
    
        
% INITIALIZE CANOPY STATES
    VARIABLES.CANOPY.gsv_sun = 0.01*ones(nl_can,nspecies);
    VARIABLES.CANOPY.gsv_shade = 0.01*ones(nl_can,nspecies);
           
% CANOPY PRECIPITATION INTERCEPTION
%   Smax [mm / LAI]
%   Sh2o_prof
%   ppt [mm]
    [Sh2o_prof, Smaxz, ppt_ground, wetfrac, dryfrac] = PRECIP_INTERCEPTION(FORCING, VARIABLES, VERTSTRUC, PARAMS);   
        VARIABLES.CANOPY.Sh2o_prof = Sh2o_prof;
        VARIABLES.CANOPY.Smaxz = Smaxz;
        VARIABLES.CANOPY.wetfrac = wetfrac;
        VARIABLES.CANOPY.dryfrac = dryfrac;
        VARIABLES.SOIL.ppt_ground = ppt_ground;
        
% PREVIOUS
% TAz_prev = VARIABLES.CANOPY.TAz;
% VARIABLES_prev = VARIABLES;      
% diver_control = false;           
    %======================================================================
    %                   SHORTWAVE RADIATION PROFILES
    %======================================================================
    [fsun, fshade, LAIsun, LAIshade, ...
     SWabs_sun, SWabs_shade, PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, ...
     PARabs_sun_lai, PARabs_shade_lai, ...
     SWabs_soil, PARabs_soil, NIRabs_soil, ...
     SWout, PARout, NIRout, PARtop, NIRtop, fdiff,...
     SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
     SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
     Kbm, taud, refl_soil] = ...
     SWRAD (FORCING, VERTSTRUC, PARAMS, CONSTANTS, VARIABLES); 
    
        % STORE VARIABLES
        VARIABLES.CANOPY.PARabs_sun = PARabs_sun_lai;
        VARIABLES.CANOPY.PARabs_shade = PARabs_shade_lai;
        VARIABLES.CANOPY.fsun = fsun;
        VARIABLES.CANOPY.fshade = fshade;
        VARIABLES.CANOPY.LAIsun = LAIsun;
        VARIABLES.CANOPY.LAIshade = LAIshade;
        VARIABLES.CANOPY.taud = taud;
        
        %Mere Add Variables for next step if replaced
        VARIABLES.CANOPY.PARabs_sun_rep = PARabs_sun;
        VARIABLES.CANOPY.PARabs_shade_rep = PARabs_shade;
        VARIABLES.CANOPY.NIRabs_sun = NIRabs_sun;
        VARIABLES.CANOPY.NIRabs_shade = NIRabs_shade;
        VARIABLES.CANOPY.SWout = SWout;
        VARIABLES.CANOPY.fdiff = fdiff;
        VARIABLES.CANOPY.SWcandir_in = SWcandir_in;
        VARIABLES.CANOPY.SWcandir_out = SWcandir_out;
        VARIABLES.CANOPY.SWcandif_in = SWcandif_in;
        VARIABLES.CANOPY.SWcandif_out = SWcandif_out;
        VARIABLES.CANOPY.SWsoildir_in = SWsoildir_in;
        VARIABLES.CANOPY.SWsoildir_out = SWsoildir_out;
        VARIABLES.CANOPY.SWsoildif_in = SWsoildif_in;
        VARIABLES.CANOPY.SWsoildif_out = SWsoildif_out;
    
        
    % LONGWAVE CONVERGENCE LOOP  
    converged_LW = 0; cnt_LW = 0; maxiters = 20; percdiff = 0.1;
    replaced = 0; %Mere
    while (~converged_LW)     
     
        % LONGWAVE RADIATION ABSORPTION    
        [LWabs_can, LWabs_canM, LWabs_sun, LWabs_shade, LWabs_soil, LWabs_soilM, LWin, LWin2, LWout, LWoutM, ...
         LWemit_can, LWemit_sun, LWemit_shade, LWemit_soil] = ...
            LWRAD (FORCING, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);          
       
        % TOTAL ABSORBED RADIATION [W/m^2 ground]
        Totabs_sun = LWabs_sun + SWabs_sun;         
        Totabs_shade = LWabs_shade + SWabs_shade;   

        % TOTAL ABSORBED RADIATION PER UNIT LEAF AREA [W/m^2 leaf area]
        Rabs_sun_lai = Totabs_sun./LAIsun;
        Rabs_sun_lai(LAIsun==0)=0;
        Rabs_shade_lai = Totabs_shade./LAIshade;
        Rabs_shade_lai(LAIshade==0)=0;
        
        % SOIL ABSORBED ENERGY
        Totabs_soil = SWabs_soil + LWabs_soil;
        Rnrad_soil = Totabs_soil - LWemit_soil;

            % ASSIGN
            VARIABLES.CANOPY.Rnrad_soil = Rnrad_soil;
            VARIABLES.CANOPY.Rabs_sun = Rabs_sun_lai;
            VARIABLES.CANOPY.Rabs_shade = Rabs_shade_lai;
            VARIABLES.SOIL.Totabs_soil = Totabs_soil;
            FORCING.LWdn = LWin;
            FORCING.LWdn2 = LWin2;
            
            %Mere Add Variables for next step if replaced
            VARIABLES.CANOPY.LWabs_canM = LWabs_canM;
            VARIABLES.CANOPY.LWabs_soilM = LWabs_soilM;
            VARIABLES.CANOPY.LWabs_can = LWabs_can;
            VARIABLES.CANOPY.LWabs_soil = LWabs_soil;
            VARIABLES.CANOPY.LWemit_soil = LWemit_soil;
            VARIABLES.CANOPY.LWemit_can = LWemit_can;
            VARIABLES.CANOPY.LWemit_sun = LWemit_sun;
            VARIABLES.CANOPY.LWemit_shade = LWemit_shade;
            VARIABLES.CANOPY.LWout = LWout; 
            VARIABLES.CANOPY.LWoutM = LWoutM; 
            
           
        % ALLOCATE THE MEMORY FOR VARIABLES.    
        ALLOCATE_LEAF_VARS;    
        %==================================================================
        %                   SHADED CANOPY SOLUTION
        %   Calculations performed per [m^2 ground area], and canopy fluxes
        %   are calculated by integrating over the shaded leaf area
        %==================================================================     
        sunlit = 0;
        
        for ii=1:nspecies
         [Ph_shade(:,ii), An_shade(:,ii), Ci_shade(:,ii), gsv_shade(:,ii), Tl_shade(:,ii), LE_shade(:,ii),...
          dHcan_shade(:,ii), TR_shade(:,ii), Evap_shade(:,ii), H_shade(:,ii), ...
          psil_shade(:,ii), fsvg_shade(:,ii),fsvm_shade(:,ii), Ch2o_shade(:,ii), gbv_shade(:,ii), gbh_shade(:,ii), ...
          Ph_limit_shade(:,ii), Jc_C3_shade(:,ii), Jj_C3_shade(:,ii), Js_C3_shade(:,ii), Jc_C4_shade(:,ii),...
          Jj_C4_shade(:,ii), Js_C4_shade(:,ii), ...
          VARIABLES,LSshaCON(ii)] = ...
                     LEAF_SOLUTION(FORCING, VARIABLES, PARAMS, CONSTANTS, VERTSTRUC, SWITCHES, sunlit, ii, tt);
        end 
                         
        %==================================================================
        %                       SUNLIT CANOPY SOLUTION
        %   Calculations performed per [m^2 leaf area], and canopy fluxes
        %       are calculated by integrating vertically over the sunlit 
        %       leaf area
        %==================================================================
        if (sum(fsun)==0)   % under nocturnal conditions all leaf area is 
                            % considered to be shaded
                            
            Ph_sun = zeros(nl_can,nspecies);
            An_sun = zeros(nl_can,nspecies);
            LE_sun = zeros(nl_can,nspecies);
            dHcan_sun = zeros(nl_can,nspecies);
            H_sun = zeros(nl_can,nspecies);
            Phtype_sun = NaN(nl_can,nspecies);
            gsv_sun = gsv_shade;
            gbv_sun = gbv_shade;
            gbh_sun = gbh_shade;
            Ci_sun = Ci_shade;
            VARIABLES.CANOPY.Cs_sun = zeros(nl_can,nspecies);
%             Tl_sun = NaN(nl_can,nspecies);
            Tl_sun = Tl_shade;   
            psil_sun = psil_shade;
            fsvg_sun = fsvg_shade;
            fsvm_sun = fsvm_shade;
            Evap_sun = zeros(nl_can,nspecies);
            TR_sun = zeros(nl_can,nspecies);
            Ch2o_sun = zeros(nl_can,nspecies);           
            
            Ph_limit_sun = zeros(nl_can,nspecies); 
            Jc_C3_sun = zeros(nl_can,nspecies); 
            Jj_C3_sun = zeros(nl_can,nspecies); 
            Js_C3_sun = zeros(nl_can,nspecies); 
            Jc_C4_sun = zeros(nl_can,nspecies); 
            Jj_C4_sun = zeros(nl_can,nspecies); 
            Js_C4_sun = zeros(nl_can,nspecies);
            LSsunCON = ones(1,nspecies);
            VARIABLES.CANOPY.Tl_sun = Tl_sun;
        else
            sunlit = 1;
            
            for ii=1:nspecies
            [Ph_sun(:,ii), An_sun(:,ii), Ci_sun(:,ii), gsv_sun(:,ii), Tl_sun(:,ii), LE_sun(:,ii),...
             dHcan_sun(:,ii), TR_sun(:,ii), Evap_sun(:,ii), H_sun(:,ii), ...
             psil_sun(:,ii), fsvg_sun(:,ii),fsvm_sun(:,ii), Ch2o_sun(:,ii), gbv_sun(:,ii), gbh_sun(:,ii), ...
             Ph_limit_sun(:,ii), Jc_C3_sun(:,ii), Jj_C3_sun(:,ii), Js_C3_sun(:,ii), Jc_C4_sun(:,ii),...
             Jj_C4_sun(:,ii), Js_C4_sun(:,ii),...
             VARIABLES,LSsunCON(ii)] = ...
                    LEAF_SOLUTION(FORCING, VARIABLES, PARAMS, CONSTANTS, VERTSTRUC, SWITCHES, sunlit, ii, tt);                                
            end
        end
        % COMPUTE TEMPERATURE IN ALL CANOPY (AVERAGING OVER THE SPECIES)
            Tl_can_shade = nansum(((Tl_shade.*LAIzall)./repmat(LAIz,1,nspecies)),2);
            Tl_can_sun = nansum(((Tl_sun.*LAIzall)./repmat(LAIz,1,nspecies)),2);
        % CHECK THOSE LAYERS WHEN THERE IS NOT LAI AT ALL AND SET Tl_can = nan    
            Tl_can_shade(nvinds) = nan;
            Tl_can_sun(nvinds) = nan;
        % CHECK AGAIN THAT ALL THE SUNLIT FLUXES ARE nan IF fsun = 0
            inddark = fsun == 0;
            Tl_can_sun(inddark) = nan;
        % SAVE CANOPY TEMPERATURES
            VARIABLES.CANOPY.Tl_can_sun = Tl_can_sun;
            VARIABLES.CANOPY.Tl_can_shade = Tl_can_shade;                            
        % ASSIGN
        
            VARIABLES.CANOPY.Ph_sun = Ph_sun;
            VARIABLES.CANOPY.Ph_shade = Ph_shade;        
            VARIABLES.CANOPY.An_sun = An_sun;
            VARIABLES.CANOPY.An_shade = An_shade;
            VARIABLES.CANOPY.LE_sun = LE_sun;
            VARIABLES.CANOPY.LE_shade = LE_shade;
            VARIABLES.CANOPY.H_sun = H_sun;
            VARIABLES.CANOPY.H_shade = H_shade;
            
            %Mere Add Variables for next step if replaced
            VARIABLES.CANOPY.TR_sun = TR_sun;
            VARIABLES.CANOPY.TR_shade = TR_shade;
            VARIABLES.CANOPY.Tl_sun = Tl_sun;
            VARIABLES.CANOPY.Tl_shade = Tl_shade;
            VARIABLES.CANOPY.psil_sun = psil_sun;
            VARIABLES.CANOPY.psil_shade = psil_shade;
            VARIABLES.CANOPY.gsv_sun = gsv_sun;
            VARIABLES.CANOPY.gsv_shade = gsv_shade;
            VARIABLES.CANOPY.gbv_sun = gbv_sun;
            VARIABLES.CANOPY.gbv_shade = gbv_shade;
            VARIABLES.CANOPY.gbh_sun = gbh_sun;
            VARIABLES.CANOPY.gbh_shade = gbh_shade;
            VARIABLES.CANOPY.fsvg_sun = fsvg_sun;
            VARIABLES.CANOPY.fsvg_shade = fsvg_shade;
            VARIABLES.CANOPY.fsvm_sun = fsvm_sun;
            VARIABLES.CANOPY.fsvm_shade = fsvm_shade;
            VARIABLES.CANOPY.Ci_sun = Ci_sun;
            VARIABLES.CANOPY.Ci_shade = Ci_shade;
            VARIABLES.CANOPY.Ph_limit_sun = Ph_limit_sun;
            VARIABLES.CANOPY.Ph_limit_shade = Ph_limit_shade;
            VARIABLES.CANOPY.Jc_C3_sun = Jc_C3_sun;
            VARIABLES.CANOPY.Jc_C3_shade = Jc_C3_shade;
            VARIABLES.CANOPY.Jj_C3_sun = Jj_C3_sun;
            VARIABLES.CANOPY.Jj_C3_shade = Jj_C3_shade;
            VARIABLES.CANOPY.Js_C3_sun = Js_C3_sun;
            VARIABLES.CANOPY.Js_C3_shade = Js_C3_shade;
            VARIABLES.CANOPY.Jc_C4_sun = Jc_C4_sun;
            VARIABLES.CANOPY.Jc_C4_shade = Jc_C4_shade;
            VARIABLES.CANOPY.Jj_C4_sun = Jj_C4_sun;
            VARIABLES.CANOPY.Jj_C4_shade = Jj_C4_shade;
            VARIABLES.CANOPY.Js_C4_sun = Js_C4_sun;
            VARIABLES.CANOPY.Js_C4_shade = Js_C4_shade;
            VARIABLES.CANOPY.LSshaCON = LSshaCON;
            VARIABLES.CANOPY.LSsunCON = LSsunCON;
            
          
        % SOIL RESPIRATION [umol CO2/ m^2 ground / s] 
            Fc_soil = Ro .* Q10.^((Tsoil - 10)/10); 
            
       % SOIL ENERGY FLUXES
       % determine the case
           zicesl = VARIABLES.SOIL.zicesl;     
           
         
        if SWITCHES.litter;         
           [H_soil, LE_soil, G, Gsl, RH_soil, Tsurf, remainsoil, dH, dS, VARIABLES] =...
           SOIL_SURFACE_FLUXES_LITTER(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES,FORCING);
        else        
           if zicesl > 0;                                  % Assumed ice to solve energy balance                               
                   [H_soil, LE_soil, G, Gsl, RH_soil, Tsurf, remainsoil, dH, dS, VARIABLES] =...
                   SOIL_SURFACE_FLUXES_SNOW(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES,FORCING);
           else                                            % Assumed no ice to solve energy balance
                  [H_soil, LE_soil, G, RH_soil, Tsurf, remainsoil, VARIABLES] = ... 
                  SOIL_SURFACE_FLUXES(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS,SWITCHES,FORCING);
                  dH = 0;
                  dS = 0;
           end
        end  
        
                % ASSIGN VARIABLES
                VARIABLES.SOIL.G = G;                
                VARIABLES.SOIL.Fc_soil = Fc_soil;
                VARIABLES.SOIL.LE_soil = LE_soil;
                VARIABLES.SOIL.H_soil = H_soil;
                VARIABLES.SOIL.Tsurf=Tsurf;
                VARIABLES.SOIL.remainsoil=remainsoil;
                
                %Mere Add Variables for next step if replaced
                VARIABLES.SOIL.RH_soil = RH_soil;
                
        out = 0;%Mere determine if outliers are replaced in this iteration        
        if (turb_on)
            if cnt_LW ==0
                Sh=zeros(length(TAz));
            end
            [CAz, EAz, TAz, tt, out, VARIABLES, Sh] = MICROENVIRONMENT(FORCING, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, tt, out, removing_turb_outliers, similarity_adjustment, cnt_LW, Sh);  %Mere
            % ASSIGN
                VARIABLES.CANOPY.CAz = CAz;
                VARIABLES.CANOPY.TAz = TAz;
                VARIABLES.CANOPY.EAz = EAz;
        end  
                           
        % TEST LONGWAVE CONVERGENCE
        cnt_LW = cnt_LW + 1; 
        if (cnt_LW>1)
            diffprof = LWabs_can-LWabs_prev;
            percdiffprof(vinds) = diffprof(vinds)./LWabs_prev(vinds);
            
            Tl_new = [Tl_can_shade Tl_can_sun];
            difftl = max((max(abs(Tl_new-Tl_prev))));

            if  isempty(percdiffprof)
                converged_LW = 1;
            end
            if (max(abs(percdiffprof(vinds))) < percdiff)
                if (difftl <= 0.5)
                    converged_LW = 1;
                end
            end
%             if (max(abs(percdiffprof(vinds))) < percdiff && difftl <= 0.5)
%                 converged_LW = 1;
%             end
            if (difftl > 10 && cnt_LW > 2)
                diverging = diverging + 1; %Mere Count the number of diverging timesteps
                disp(['*** DIVERGENCE IN CANOPY MODEL!!! --> Timestep:', num2str(tt)]);
%                 redo_noturb = 1; 
                replaced = 1; %Mere Indicates we will replace this timestep with the previous one
                converged_LW = 1; %Indicates that we are leaving the loop
            end
        end
        LWabs_prev = LWabs_can;
        Tl_prev = [Tl_can_shade Tl_can_sun];
        
        if (cnt_LW>maxiters && converged_LW==0)
            countingmaxedout = countingmaxedout + 1; %Mere Count the number of timesteps exceeding max # of iterations
            disp(['*** TOO MANY ITERATIONS IN CANOPY MODEL!!! --> Timestep:', num2str(tt)]);%VARIABLES.niters_driver)]); %Mere uncommented
%             redo_noturb = 1; 
            replaced = 1; %Mere Indicates we will replace this timestep with the previous one
            converged_LW = 1; %Indicates that we are leaving the loop
        end      
        
    end %LONGWAVE ITERATION
    
    %Notes that outliers were removed during turbulence in this time step
    if out > 0
        outlier = outlier + 1;
    end

%%  %%%%%%%%%%%%%%%%% ALL VARIABLES CALCULATED REGULARLY %%%%%%%%%%%%%%%%%%     

if (~replaced) %Mere  
    
   
    % NET RADIATION
    if (~isnan(FORCING.LWdn))
           Rnrad_eco = (FORCING.Rg - SWout) + (FORCING.LWdn - LWout);          %[W/m^2] Net Radiation in the Ecosystem 
    else
           Rnrad_eco = (FORCING.Rg - SWout) + (LWin - LWout);                  %[W/m^2] Net Radiation in the Ecosystem                     
    end
    Rnrad_sun = SWabs_sun + LWabs_sun - LWemit_sun;                            %[W/m^2] Net Radiation Canopy Sunlit Fraction  
    Rnrad_shade = SWabs_shade + LWabs_shade - LWemit_shade;                    %[W/m^2] Net Radiation Canopy Shade Fraction

        % H2O storage on foliage --> Precipitation and Condensation
    Evap_prof = (sum(Evap_sun.*fLAIz,2).*LAIsun + sum(Evap_shade.*fLAIz,2).*LAIshade) * dtime;  % [mm] Evaporation at each layer
    Evap_can = sum(Evap_prof);                                                 % [mm] Total canopy evaporation

    Ch2o_prof = -(sum(Ch2o_sun.*fLAIz,2).*LAIsun + sum(Ch2o_shade.*fLAIz,2).*LAIshade) * dtime; % [mm] Condensation at each layer
    Ch2o_can = sum(Ch2o_prof);                                                 % [mm] Total canopy condensation

       % ASSIGN
    VARIABLES.CANOPY.Ch2o_prof = Ch2o_prof;                                    % [mm] Condensation Canopy at every Layer 
    VARIABLES.CANOPY.Evap_prof = Evap_prof;                                    % [mm] Evaporation from Canopy at every Layer 
    VARIABLES.CANOPY.Ch2o_can = Ch2o_can;                                      % [mm] Total Condensation Canopy 
    VARIABLES.CANOPY.Evap_can = Evap_can;                                      % [mm] Total Evaporation from Canopy 

        % Adjust Canopy Water Storage
    [Sh2o_prof, Sh2o_can, VARIABLES, Evap_can, Evap_prof, ppt_ground] = EVAP_CONDENSATION_ADJUST(VARIABLES, VERTSTRUC, PARAMS);

        % COMPUTE CANOPY TOTAL FLUXES

    Ph_can = sum(sum(Ph_sun.*fLAIz,2).*LAIsun) + sum(sum(Ph_shade.*fLAIz,2).*LAIshade);      %[umol/m^2/s] Total Photosynthesis canopy  
    Ph_can_all = sum(Ph_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
                + sum(Ph_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[umol/m^2/s] Total Photosynthesis canopy for all  species

    An_can = sum(sum(An_sun.*fLAIz,2).*LAIsun) + sum(sum(An_shade.*fLAIz,2).*LAIshade);      %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy 
    An_can_all = sum(An_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
                + sum(An_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy for species

    LE_can = sum(sum(LE_sun.*fLAIz,2).*LAIsun) + sum(sum(LE_shade.*fLAIz,2).*LAIshade);      %[W /m^2] Total Latent Heat canopy     
    LE_can_all = sum(LE_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
                + sum(LE_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[W /m^2] Total Latent Heat canopy for all species  

    LE_can_all_lay_shade = LE_shade.*fLAIz.*(repmat(LAIshade,1,nspecies));                   %[W /m^2] Total Latent Heat canopy for all layers and species, shade
    LE_can_all_lay_sun = LE_sun.*fLAIz.*(repmat(LAIsun,1,nspecies));                         %[W /m^2] Total Latent Heat canopy for all layers and species, sunlit

    H_can = sum(sum(H_sun.*fLAIz,2).*LAIsun) + sum(sum(H_shade.*fLAIz,2).*LAIshade);         %[W /m^2] Total Sensible Heat canopy   
    H_can_all = sum(H_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
                + sum(H_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                        %[W /m^2] Total Sensible Heat canopy for all species     

    H_can_all_lay_shade = H_shade.*fLAIz.*(repmat(LAIshade,1,nspecies));                     %[W /m^2] Total Sensible Heat canopy for all layers and species, shade
    H_can_all_lay_sun = H_sun.*fLAIz.*(repmat(LAIsun,1,nspecies));                           %[W /m^2] Total Sensible Heat canopy for all layers and species, sunlit


    dHcan = sum(sum(dHcan_sun.*fLAIz,2).*LAIsun) + sum(sum(dHcan_shade.*fLAIz,2).*LAIshade); %[W /m^2] Change in heat storage in canopy
    dHcan_all = sum(dHcan_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...                       
                + sum(dHcan_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                    %[W /m^2] Change in heat storage in canopy for all species 

    dHcan_all_lay_shade = dHcan_shade.*fLAIz.*(repmat(LAIshade,1,nspecies));                 %[W /m^2] Change in heat storage in canopy for all layers and species, shade fraction   
    dHcan_all_lay_sun = dHcan_sun.*fLAIz.*(repmat(LAIsun,1,nspecies));                       %[W /m^2] Change in heat storage in canopy for all layers and species, sunlit fraction


    TR_can = sum(sum(TR_sun.*fLAIz,2).*LAIsun) + sum(sum(TR_shade.*fLAIz,2).*LAIshade);      %[mm/s] = [g/m^2/s] Total Transpiration from Canopy
    TR_can_all = sum(TR_sun.*fLAIz.*(repmat(LAIsun,1,nspecies)))...
                + sum(TR_shade.*fLAIz.*(repmat(LAIshade,1,nspecies)));                       %[mm/s] = [g/m^2/s] Total Transpiration from Canopy for all species

    Rnrad_can = sum(Rnrad_sun(vinds)) + sum(Rnrad_shade(vinds));                             %[mm/s] Net Radiation Canopy 

        
    % COMPUTE REMAINDER
   
    %Mere added in the loop; dHcan needs to be included in remainder
    if (canstorheat)
        remaincan = Rnrad_can - H_can - LE_can - An_can*0.506 - dHcan;                       % [W / m^2] Energy Balance Error in Canopy 
        remaineco = Rnrad_eco - (H_can + LE_can + An_can*0.506 + dHcan) ...
                            - (H_soil  + LE_soil + G + dH + dS);                             % [W / m^2] Energy Balance Error in all Ecosystem
    else
        remaincan = Rnrad_can - H_can - LE_can - An_can*0.506;
        remaineco = Rnrad_eco - (H_can + LE_can + An_can*0.506) ...
                            - (H_soil  + LE_soil + G + dH + dS); 
    end

    
    %ASSIGN   

    VARIABLES.CANOPY.Ph_can = Ph_can;                                          %[umol/m^2/s] Total Photosynthesis canopy
    VARIABLES.CANOPY.Ph_can_all = Ph_can_all;                                  %[umol/m^2/s] Total Photosynthesis canopy for all  species   

    VARIABLES.CANOPY.An_can = An_can;                                          %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy
    VARIABLES.CANOPY.An_can_all = An_can_all;                                  %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy for species                          

    VARIABLES.CANOPY.LE_can = LE_can;                                          %[W /m^2] Total Latent Heat canopy 
    VARIABLES.CANOPY.LE_can_all = LE_can_all;                                  %[W /m^2] Total Latent Heat canopy for all species 
    VARIABLES.CANOPY.LE_can_all_lay_shade = LE_can_all_lay_shade;              %[W /m^2] Total Latent Heat canopy for all layers and species, shade 
    VARIABLES.CANOPY.LE_can_all_lay_sun = LE_can_all_lay_sun;                  %[W /m^2] Total Latent Heat canopy for all layers and species, sunlit 


    VARIABLES.CANOPY.H_can = H_can;                                            %[W /m^2] Total Sensible Heat canopy 
    VARIABLES.CANOPY.H_can_all = H_can_all;                                    %[W /m^2] Total Sensible Heat canopy for all species    
    VARIABLES.CANOPY.H_can_all_lay_shade = H_can_all_lay_shade;                %[W /m^2] Total Sensible Heat canopy for all layers and species, shade 
    VARIABLES.CANOPY.H_can_all_lay_sun = H_can_all_lay_sun;                    %[W /m^2] Total Sensible Heat canopy for all layers and species, sunlit 

    VARIABLES.CANOPY.Tl_prev_dt = (nansum(Tl_sun.*fLAIz,2)).*fsun + ...
                                    (nansum(Tl_shade.*fLAIz,2)).*fshade;       %[C] Temperature canopy for all layers  

    VARIABLES.CANOPY.dHcan = dHcan;                                            %[W /m^2] Change in heat storage in canopy 
    VARIABLES.CANOPY.dHcan_all = dHcan_all;                                    %[W /m^2] Change in heat storage in canopy for all species   
    VARIABLES.CANOPY.dHcan_all_lay_shade = dHcan_all_lay_shade;                %[W /m^2] Total Sensible Heat canopy for all layers and species, shade 
    VARIABLES.CANOPY.dHcan_all_lay_sun = dHcan_all_lay_sun;                    %[W /m^2] Total Sensible Heat canopy for all layers and species, sunlit 

    VARIABLES.CANOPY.TR_can = TR_can;                                          %[mm/s] = [g/m^2/s] Total Transpiration from Canopy  
    VARIABLES.CANOPY.TR_can_all = TR_can_all;                                  %[mm/s] = [g/m^2/s] Total Transpiration from Canopy for all species
    VARIABLES.CANOPY.Evap_can = Evap_can;                                      %[mm] Total canopy evaporation

    VARIABLES.CANOPY.Sh2o_prof = Sh2o_prof;                                    %[mm] Canopy moisture storage for all layers                                      
    VARIABLES.CANOPY.Sh2o_can = Sh2o_can;                                      %[mm] Canopy canopy moisture storage  

    VARIABLES.SOIL.Tlprev = VARIABLES.SOIL.Tli;                                % [C] Temperature Snow-Litter pack   
    VARIABLES.refl_soil = refl_soil;                                           % [] SW Reflection from Soil
    
    %Mere Add Variables for next step if replaced
    VARIABLES.CANOPY.Rnrad_eco = Rnrad_eco;
    VARIABLES.CANOPY.Rnrad_sun = Rnrad_sun;
    VARIABLES.CANOPY.Rnrad_shade = Rnrad_shade;
    VARIABLES.CANOPY.Rnrad_can = Rnrad_can;
    VARIABLES.CANOPY.remaincan = remaincan;
    VARIABLES.CANOPY.remaineco = remaineco;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Check for large remainder %%%%%%%%%%%%%%%%%%%%%%
    %Mere If remainder is large, replace with previous timestep

    if (replace_remainder_errors)
        if tt>1
            xbar_remaincan = VARIABLES.CANOPY.xbar_remaincan; 
            sigma2_remaincan = VARIABLES.CANOPY.sigma2_remaincan;           
            remaincan_h = VARIABLES.CANOPY.remaincan_history;
        end


        %Create Timestep n
        n = tt;
        starter_threshold = 50; %Mere Timesteps exceeding this will be replaced in first calc_set timesteps.
        calc_set = 100; %Mere size of set used to calculate mean/variance for removal 
        set_i = calc_set - 1;
        num_sd = 4; %Mere number of standard deviations from mean to remove outliers for energy balance issues

        % Initialize variables
        if n==1
            sigma2_remaincan = 0;
            xbar_remaincan = remaincan;
        elseif n==2
            sigma2_remaincan = 1/(n-1).*(remaincan-xbar_remaincan).^2; %variance
            xbar_remaincan = 1/(n).*(remaincan + remaincan_h(n-1)); %mean           
        elseif n>2 && n<=calc_set
            %check for outliers before variance calculated
            if abs(remaincan) > starter_threshold            
                largeremainder = largeremainder + 1; %Mere Count the number of timesteps with large remainder
                disp(['*** LARGE REMAINDER!!! --> Timestep:', num2str(tt)]);
                replaced = 1; %Mere Indicates we will replace this timestep with the previous one
            end

            %set variance and mean for next time step
            sigma2_remaincan = (n-2)./(n-1).*sigma2_remaincan + (1./(n-1).*(remaincan - ((n-1)/n).*(xbar_remaincan + 1/(n-1).*remaincan)).^2); %variance function of previous
            xbar_remaincan = (n-1)/n.*(xbar_remaincan + 1/(n-1).*remaincan); %mean as a function of previous mean   
        end

        % Remove Outliers
        if n>calc_set            
           sigma_remaincan = sqrt(sigma2_remaincan); %standard deviation

           % Outlier if > or < 99.9 percentile
           if remaincan > xbar_remaincan + num_sd*sigma_remaincan
              largeremainder = largeremainder + 1;
              disp(['*** LARGE REMAINDER!!! --> Timestep:', num2str(tt)]);
              replaced = 1; 
           end
           if remaincan < xbar_remaincan - num_sd*sigma_remaincan
              largeremainder = largeremainder + 1;
              disp(['*** LARGE REMAINDER!!! --> Timestep:', num2str(tt)]);
              replaced = 1;
           end

           % Set mean and variance based on previous calc_set time steps
            sumr=0;
            difr=0;
            for h=n-set_i:n-1
                sumr=remaincan_h(h)+sumr;
            end
            xbar_remaincan = (sumr+remaincan)/calc_set;
            for h=n-set_i:n-1
                difr=difr+(remaincan_h(h)-xbar_remaincan)^2;
            end
            difr=difr+(remaincan-xbar_remaincan)^2;
            sigma2_remaincan = difr/set_i;
        end            

    %   %uncomment if replacing all values based on a single threshold 
    %   %(then must comment out previous outlier calcs)
        %     if (abs(remaincan) > starter_threshold)
        %         largeremainder = largeremainder + 1;
        %         disp(['*** LARGE REMAINDER!!! --> Timestep:', num2str(tt)]);
        %         replaced = 1;
        %     end

        %Assign Variables
            VARIABLES.CANOPY.xbar_remaincan = xbar_remaincan; 
            VARIABLES.CANOPY.sigma2_remaincan = sigma2_remaincan;
            VARIABLES.CANOPY.remaincan_history(tt) = remaincan; 

    end
end


%% %%%%%%%%%% VARIABLES CONVERTED TO PREVIOUS TIMESTEP %%%%%%%%%%%%%%%%%

if (replaced == 1)
    
    
    VARIABLES = VARIABLES_last_tt;
    if (turb_on) && (removing_turb_outliers)
        VARIABLES.CANOPY.TAz_history(:,tt) = VARIABLES.CANOPY.TAz_history(:,tt-1); %Mere need for outliers calc for next step
        VARIABLES.CANOPY.CAz_history(:,tt) = VARIABLES.CANOPY.CAz_history(:,tt-1);
        VARIABLES.CANOPY.EAz_history(:,tt) = VARIABLES.CANOPY.EAz_history(:,tt-1);             
    end
    VARIABLES.timestep = tt; %Kunxuan: reset timestep to current, not last_tt.
    
    %ASSIGN
    
    %Mere These variables are not called again after leaving the function
    %before further asssignment and calculations. Need to be reassigned
    %here.
    
    dHcan = VARIABLES.CANOPY.dHcan; 
    Ph_can = VARIABLES.CANOPY.Ph_can;                                          %[umol/m^2/s] Total Photosynthesis canopy
    An_can = VARIABLES.CANOPY.An_can;                                          %[umol/m^2/s] Total Photosynthesis minus leaf respiration canopy
    LE_can = VARIABLES.CANOPY.LE_can;                                          %[W /m^2] Total Latent Heat canopy 
    H_can = VARIABLES.CANOPY.H_can;                                            %[W /m^2] Total Sensible Heat canopy    
    TR_can = VARIABLES.CANOPY.TR_can;                                          %[mm/s] = [g/m^2/s] Total Transpiration from Canopy  
    Evap_can = VARIABLES.CANOPY.Evap_can;                                      %[mm] Total canopy evaporation
    Ch2o_prof = VARIABLES.CANOPY.Ch2o_prof;                                    % [mm] Condensation Canopy at every Layer 
    Evap_prof = VARIABLES.CANOPY.Evap_prof;                                    % [mm] Evaporation from Canopy at every Layer 
    Ch2o_can = VARIABLES.CANOPY.Ch2o_can;                                      % [mm] Total Condensation Canopy
    Sh2o_prof = VARIABLES.CANOPY.Sh2o_prof;                                    %[mm] Canopy moisture storage for all layers                                      
    Sh2o_can = VARIABLES.CANOPY.Sh2o_can;                                      %[mm] Canopy canopy moisture storage     
    CAz = VARIABLES.CANOPY.CAz;                                            % [umol/mol] Atmospheric Concentration of CO2, for all layers
    TAz = VARIABLES.CANOPY.TAz;                                            % [C] Air Temperature, for all layers 
    EAz = VARIABLES.CANOPY.EAz;                                            % [kPa] Vapor Pressure Air, for all layers 
    An_sun = VARIABLES.CANOPY.An_sun;
    An_shade = VARIABLES.CANOPY.An_shade;
    LE_sun = VARIABLES.CANOPY.LE_sun;
    LE_shade = VARIABLES.CANOPY.LE_shade;
    H_sun = VARIABLES.CANOPY.H_sun;
    H_shade = VARIABLES.CANOPY.H_shade;
    Tsurf = VARIABLES.SOIL.Tsurf;
    remainsoil = VARIABLES.SOIL.remainsoil;
    G = VARIABLES.SOIL.G;     
    fsun = VARIABLES.CANOPY.fsun;
    fshade = VARIABLES.CANOPY.fshade;
    LAIsun = VARIABLES.CANOPY.LAIsun;
    LAIshade = VARIABLES.CANOPY.LAIshade;   
    Fc_soil = VARIABLES.SOIL.Fc_soil;
    LE_soil = VARIABLES.SOIL.LE_soil;
    H_soil = VARIABLES.SOIL.H_soil;
    Rnrad_soil = VARIABLES.CANOPY.Rnrad_soil;
    Vz = VARIABLES.CANOPY.Vz;
    Uz = VARIABLES.CANOPY.Uz;
    wetfrac = VARIABLES.CANOPY.wetfrac;
    dryfrac = VARIABLES.CANOPY.dryfrac;
    ppt_ground = VARIABLES.SOIL.ppt_ground;
    
    Rnrad_eco = VARIABLES.CANOPY.Rnrad_eco;
    Rnrad_sun = VARIABLES.CANOPY.Rnrad_sun;
    Rnrad_shade = VARIABLES.CANOPY.Rnrad_shade;
    Rnrad_can = VARIABLES.CANOPY.Rnrad_can;
    remaincan = VARIABLES.CANOPY.remaincan;
    remaineco = VARIABLES.CANOPY.remaineco;
    
    if (replace_remainder_errors)
        %Mere Calculate Remainder for next step
        xbar_remaincan = VARIABLES.CANOPY.xbar_remaincan; 
        sigma2_remaincan = VARIABLES.CANOPY.sigma2_remaincan;           
        remaincan_h = VARIABLES.CANOPY.remaincan_history;
        n=tt;
        calc_set = 100; %Mere size of set used to calculate mean/variance for removal 
        set_i = calc_set - 1;

        if n > calc_set
            % Set mean and variance based on previous calc_set time steps
            sumr=0;
            difr=0;
            for h=n-set_i:n-1
                sumr=remaincan_h(h)+sumr;
            end
            xbar_remaincan = (sumr+remaincan)/calc_set;
            for h=n-set_i:n-1
                difr=difr+(remaincan_h(h)-xbar_remaincan)^2;
            end
            difr=difr+(remaincan-xbar_remaincan)^2;
            sigma2_remaincan = difr/set_i;
        else
            sigma2_remaincan = (n-2)./(n-1).*sigma2_remaincan + (1./(n-1).*(remaincan - ((n-1)/n).*(xbar_remaincan + 1/(n-1).*remaincan)).^2); %variance function of previous
            xbar_remaincan = (n-1)/n.*(xbar_remaincan + 1/(n-1).*remaincan); %mean as a function of previous mean   
        end

        %Assign Variables
            VARIABLES.CANOPY.xbar_remaincan = xbar_remaincan; 
            VARIABLES.CANOPY.sigma2_remaincan = sigma2_remaincan;
            VARIABLES.CANOPY.remaincan_history(tt) = remaincan;
    end
    
    %Mere Add Variables for next step if replaced
        PARabs_sun = VARIABLES.CANOPY.PARabs_sun_rep;
        PARabs_shade = VARIABLES.CANOPY.PARabs_shade_rep ;
        NIRabs_sun = VARIABLES.CANOPY.NIRabs_sun;
        NIRabs_shade = VARIABLES.CANOPY.NIRabs_shade;
        SWout = VARIABLES.CANOPY.SWout;
        fdiff = VARIABLES.CANOPY.fdiff;
        SWcandir_in = VARIABLES.CANOPY.SWcandir_in;
        SWcandir_out = VARIABLES.CANOPY.SWcandir_out;
        SWcandif_in = VARIABLES.CANOPY.SWcandif_in;
        SWcandif_out = VARIABLES.CANOPY.SWcandif_out;
        SWsoildir_in = VARIABLES.CANOPY.SWsoildir_in;
        SWsoildir_out = VARIABLES.CANOPY.SWsoildir_out;
        SWsoildif_in = VARIABLES.CANOPY.SWsoildif_in;
        SWsoildif_out = VARIABLES.CANOPY.SWsoildif_out;
        
        LWabs_canM = VARIABLES.CANOPY.LWabs_canM;
        LWabs_soilM = VARIABLES.CANOPY.LWabs_soilM;
        LWabs_can = VARIABLES.CANOPY.LWabs_can;
        LWabs_soil = VARIABLES.CANOPY.LWabs_soil;
        LWemit_soil = VARIABLES.CANOPY.LWemit_soil;
        LWemit_can = VARIABLES.CANOPY.LWemit_can;
        LWemit_sun = VARIABLES.CANOPY.LWemit_sun;
        LWemit_shade = VARIABLES.CANOPY.LWemit_shade;
        LWout = VARIABLES.CANOPY.LWout; 
        LWoutM = VARIABLES.CANOPY.LWoutM; 
        
        TR_sun = VARIABLES.CANOPY.TR_sun;
        TR_shade = VARIABLES.CANOPY.TR_shade;
        Tl_sun = VARIABLES.CANOPY.Tl_sun;
        Tl_shade = VARIABLES.CANOPY.Tl_shade;
        psil_sun = VARIABLES.CANOPY.psil_sun;
        psil_shade = VARIABLES.CANOPY.psil_shade;
        gsv_sun = VARIABLES.CANOPY.gsv_sun;
        gsv_shade = VARIABLES.CANOPY.gsv_shade;
        gbv_sun = VARIABLES.CANOPY.gbv_sun;
        gbv_shade = VARIABLES.CANOPY.gbv_shade;
        gbh_sun = VARIABLES.CANOPY.gbh_sun;
        gbh_shade = VARIABLES.CANOPY.gbh_shade;
        fsvg_sun = VARIABLES.CANOPY.fsvg_sun;
        fsvg_shade = VARIABLES.CANOPY.fsvg_shade;
        fsvm_sun = VARIABLES.CANOPY.fsvm_sun;
        fsvm_shade = VARIABLES.CANOPY.fsvm_shade;
        Ci_sun = VARIABLES.CANOPY.Ci_sun;
        Ci_shade = VARIABLES.CANOPY.Ci_shade;
        Ph_limit_sun = VARIABLES.CANOPY.Ph_limit_sun;
        Ph_limit_shade = VARIABLES.CANOPY.Ph_limit_shade;
        Jc_C3_sun = VARIABLES.CANOPY.Jc_C3_sun;
        Jc_C3_shade = VARIABLES.CANOPY.Jc_C3_shade;
        Jj_C3_sun = VARIABLES.CANOPY.Jj_C3_sun;
        Jj_C3_shade = VARIABLES.CANOPY.Jj_C3_shade;
        Js_C3_sun = VARIABLES.CANOPY.Js_C3_sun;
        Js_C3_shade = VARIABLES.CANOPY.Js_C3_shade;
        Jc_C4_sun = VARIABLES.CANOPY.Jc_C4_sun;
        Jc_C4_shade = VARIABLES.CANOPY.Jc_C4_shade;
        Jj_C4_sun = VARIABLES.CANOPY.Jj_C4_sun;
        Jj_C4_shade = VARIABLES.CANOPY.Jj_C4_shade;
        Js_C4_sun = VARIABLES.CANOPY.Js_C4_sun;
        Js_C4_shade = VARIABLES.CANOPY.Js_C4_shade;
        LSshaCON = VARIABLES.CANOPY.LSshaCON;
        LSsunCON = VARIABLES.CANOPY.LSsunCON;
        RH_soil = VARIABLES.SOIL.RH_soil;     
    
end

           
                
                 