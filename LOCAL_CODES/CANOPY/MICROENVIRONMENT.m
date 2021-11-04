% function [CAz, EAz, TAz, znc] = MICROENVIRONMENT(VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, tt)
function [CAz, EAz, TAz, tt, outlier, VARIABLES, Sh] = MICROENVIRONMENT(FORCING, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, tt, outlier, removing_turb_outliers, similarity_adjustment, cnt_LW, Sh)
%=========================================================================
% Calculate the canopy microenvironment (Ca, Ta, ea) using a first-order
% canopy closure model
%
% Written By: Darren Drewry, Modified by Meredith Richardson
% Last updated by Meredith Richardson, May 26, 2018
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING         % FORCING structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%       CONSTANTS       % CONSTANTS structure
%------------------------- Output Variables ------------------------------
%       CAz             %[umol/mol] Atmospheric Concentration of CO2, for all layers
%       EAz             %[kPa] Vapor Pressure Air, for all layers
%       TAz             %[C] Air Temperature, for all layers
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
% Calculate the canopy microenvironment (Ca, Ta, ea) using a first-order
% canopy closure model
%
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES
    CAz = VARIABLES.CANOPY.CAz;                                            % [umol/mol] Atmospheric Concentration of CO2, for all layers 
    TAz = VARIABLES.CANOPY.TAz;                                            % [C] Air Temperature, for all layers 
    EAz = VARIABLES.CANOPY.EAz;                                            % [kPa] Vapor Pressure Air, for all layers 
    PAz = VARIABLES.CANOPY.PAz;                                            % [kPa] Pressure Air, for all layers  
    Km = VARIABLES.CANOPY.Km;                                              % Momentum Diffusivity Constant      
    An_sun = VARIABLES.CANOPY.An_sun;                                      % [umol CO2/ m^2 leaf / s] Photosynthetic minus Leaf Respiration flux from canopy sunlit  
    An_shade = VARIABLES.CANOPY.An_shade;                                  % [umol CO2/ m^2 leaf / s] Photosynthetic minus Leaf Respiration flux from canopy shade
    LE_sun = VARIABLES.CANOPY.LE_sun;                                      % [W / m^2 leaf] LE flux from canopy sunlit     
    LE_shade = VARIABLES.CANOPY.LE_shade;                                  % [W / m^2 leaf] LE flux from canopy shade  
    H_sun = VARIABLES.CANOPY.H_sun;                                        % [W / m^2 leaf] H flux from canopy shade     
    H_shade = VARIABLES.CANOPY.H_shade;                                    % [W / m^2 leaf] H flux from canopy shade 
    LAIsun = VARIABLES.CANOPY.LAIsun;                                      % [m^2 leaf / m^2 ground] LAI sunlit
    LAIshade = VARIABLES.CANOPY.LAIshade;                                  % [m^2 leaf / m^2 ground] LAI shade
    Fc_soil = VARIABLES.SOIL.Fc_soil;                                      % [m^2 leaf / m^2 ground] CO2 flux from soil
    LE_soil = VARIABLES.SOIL.LE_soil;                                      % [m^2 leaf / m^2 ground] LE flux from soil 
    H_soil = VARIABLES.SOIL.H_soil;                                        % [m^2 leaf / m^2 ground] H flux from soil 
    
    if tt>1 && (removing_turb_outliers)
        xbar_TAz = VARIABLES.CANOPY.xbar_TAz;                              %Mere Mean of previous calc_set of TAz
        xbar_CAz = VARIABLES.CANOPY.xbar_CAz;                              %Mere Mean of previous calc_set of CAz   
        xbar_EAz = VARIABLES.CANOPY.xbar_EAz;                              %Mere Mean of previous calc_set of EAz   
        sigma2_TAz = VARIABLES.CANOPY.sigma2_TAz;                          %Mere Variance of previous calc_set of TAz
        sigma2_CAz = VARIABLES.CANOPY.sigma2_CAz;                          %Mere Variance of previous calc_set of CAz   
        sigma2_EAz = VARIABLES.CANOPY.sigma2_EAz;                          %Mere Variance of previous calc_set of EAz   
        TAz_h = VARIABLES.CANOPY.TAz_history;                              %Mere History of TAz
        CAz_h = VARIABLES.CANOPY.CAz_history;                              %Mere History of CAz
        EAz_h = VARIABLES.CANOPY.EAz_history;                              %Mere History of EAz
    elseif tt==1 && (removing_turb_outliers)
        TAz_h(:,tt) = TAz;                                                  %Mere Intialize History of TAz
        CAz_h(:,tt) = CAz;                                                  %Mere Initialize History of CAz
        EAz_h(:,tt) = EAz;                                                  %Mere Initialize History of EAz
    end
    if cnt_LW > 0 
        Sh_prev = Sh;
    end
    
    % VERTSTRUC
    znc = VERTSTRUC.znc;                                                   % [m] Height of canopy levels 
    dzc = VERTSTRUC.dzc;                                                   % [m] Thickness of layers in canopy  
    fLAIz = VERTSTRUC.fLAIz;                                               % [] Fraction of total LAI for each species at each canopy layer    
            
    % PARAMS
    hcan = PARAMS.CanStruc.hcan;                                           % [h] Canopy Height  
    
    % CONSTANTS
    Lv = CONSTANTS.Lv;                                                     % Latent heat of vaporization [J/mol]
    Lv_kg = CONSTANTS.Lv_kg;                                               % Latent heat of vaporization [J/kg]
    cp_mol = CONSTANTS.cp_mol;                                             %[J/mol/K] specific heat of air at constant pressure      
    
%*************************************************************************
%*************************************************************************

    molar_density = 44.6 * PAz * 273.15 ./ (101.3 * (TAz + 273.15));
    %psy = 6.66 * 10^-4; % [1/C]
    if tt == 97
        tt;
    end

     % CO2
        %Ca = CAz .* molar_density;  % [umol / m^3]
        Sc = ((sum(An_sun.*fLAIz,2).*LAIsun) + (nansum(An_shade.*fLAIz,2).*LAIshade))./dzc;
        Sc = -Sc./molar_density;        % [umol/mol / s]
        [CAz] = ORDER_1_CLOSURE_ALL ( CAz, znc, dzc, Km, Sc, Fc_soil, hcan ); 

     % VAPOR
        q = (EAz./PAz).*molar_density;                 
        Sv = ((sum(LE_sun.*fLAIz,2).*LAIsun) + (sum(LE_shade.*fLAIz,2).*LAIshade))./dzc; % [W / m^3]
        Sv = (Sv./Lv);
        Sv_soil =  LE_soil/Lv;
        if (similarity_adjustment) && cnt_LW==0
            case_var=3;
            E = sum((sum(LE_sun.*fLAIz,2).*LAIsun+sum(LE_shade.*fLAIz,2).*LAIshade))/Lv_kg;          %E = [kg/m2/s]
%             E_layers=sum(LE_sun.*fLAIz,2)+sum(LE_shade.*fLAIz,2);             %E = [kg/m2/s]
%             E_top=E_layers(size(E_layers,1));
            [q] = SIMILARITY_ADJ(FORCING, PARAMS, CONSTANTS, -9999, E, q, case_var);    % Esther: 10/26/18 
        end
        [q] = ORDER_1_CLOSURE_ALL ( q, znc, dzc, Km, Sv, Sv_soil, hcan );
        EAz = (q ./ molar_density) .* PAz;
        
     % HEAT
        %heat = TAz .* molar_density .* psy;
        Sh = ((sum(H_sun.*fLAIz,2).*LAIsun) + (sum(H_shade.*fLAIz,2).*LAIshade))./dzc; % [W / m^3]        
        Sh = Sh./cp_mol./molar_density;
        Sh_soil = H_soil./cp_mol./molar_density(1); 
        if (similarity_adjustment) && cnt_LW==0
            case_var=2;
            H_layers=sum(H_sun.*fLAIz,2).*LAIsun+sum(H_shade.*fLAIz,2).*LAIshade;             %H = [W/m2]
            H=sum(H_layers);
%             H_top=H_layers(size(H_layers,1));
            [TAz] = SIMILARITY_ADJ(FORCING, PARAMS, CONSTANTS, H, -9999, TAz, case_var);    % Esther: 10/26/18
        end
%         kinds = find(Km < .01);
%         Km (kinds) = .01;
%         if cnt_LW > 0
%             sinds = find( abs(Sh-Sh_prev) > .005);
%             if sum (sinds >0)
%                 sinds;
%                  if abs(Sh_soil) > .0001
%                     Sh_soil = .0001;
%                 end
%             end
%             Sh(sinds) = (Sh(sinds)+Sh_prev(sinds))./2;
%         end
       
        [TAz] = ORDER_1_CLOSURE_ALL ( TAz, znc, dzc, Km, Sh, Sh_soil, hcan );
        
        
 %% Remove Outliers - Statistical Variance - Last calc_set Values %Mere
 
 if (removing_turb_outliers)
        
        % Inputs
        starter_threshold_T = 55;
        starter_threshold_C = 1200;
        starter_threshold_E = 4;
        num_sd = 4;                 % number of standard deviations from mean to remove outliers for energy balance issues
        calc_set = 20000;             % size of set used to calculate mean/variance for removal 
        display_turb_tt = 0;        % display turbulence replacements in command window; [1] yes [0] no
        
        %Create Timestep n
        n = tt;
        set_i = calc_set - 1;
        layers = size(TAz,1); 
        
        
        % Initialize variables
        if n==1
            tinds = find(TAz > starter_threshold_T); %Mere
            TAz(tinds)=TAz_h(tinds);
            cinds = find(CAz > starter_threshold_C); %Mere
            CAz(cinds)=CAz_h(cinds);
            einds = find(EAz > starter_threshold_E); %Mere
            EAz(einds)=EAz_h(einds);
            outlier = outlier + length(tinds) + length(cinds) + length(einds);
            if (outlier > 0); disp([num2str(outlier),' Outliers @ Timestep:', num2str(tt)]); end
            
            sigma2_TAz = 0;
            sigma2_CAz = 0;
            sigma2_EAz = 0;
            xbar_TAz = TAz;
            xbar_CAz = CAz;
            xbar_EAz = EAz;
        elseif n==2
            for i=1:layers
                if abs(TAz(i)) > starter_threshold_T
%                     if (display_turb_tt); disp(['TAz replaced: ',num2str(TAz(i)),' @ Timestep:', num2str(tt)]); end
%                     TAz(i) = TAz_h(i,n-1); %Replace with previous timestep
                    TAz = TAz_h(:,n-1);
                    i=layers;
                    outlier = outlier + 1; %Mere Count the number of diverging timesteps
                end
            end
            for i=1:layers
                if CAz(i) > starter_threshold_C
                    if (display_turb_tt); disp(['CAz replaced: ',num2str(CAz(i)),' @ Timestep:', num2str(tt)]); end
                    CAz = CAz_h(:,n-1); %Replace with previous timestep
                    i=layers;
                    outlier = outlier + 1; %Mere Count the number of diverging timesteps
                end
            end
            for i=1:layers
                if EAz(i) > starter_threshold_E
                    if (display_turb_tt); disp(['EAz replaced: ',num2str(EAz(i)),' @ Timestep:', num2str(tt)]); end
                    EAz = EAz_h(:,n-1); %Replace with previous timestep
                    i=layers;
                    outlier = outlier + 1; %Mere Count the number of diverging timesteps
                end
            end
            
            sigma2_TAz = 1/(n-1).*(TAz-xbar_TAz).^2; %variance
            sigma2_CAz = 1/(n-1).*(CAz-xbar_CAz).^2;
            sigma2_EAz = 1/(n-1).*(EAz-xbar_EAz).^2;
            xbar_TAz = 1/(n).*(TAz + TAz_h(:,n-1)); %mean
            xbar_CAz = 1/(n).*(CAz + CAz_h(:,n-1));
            xbar_EAz = 1/(n).*(EAz + EAz_h(:,n-1));            
        elseif n>2 && n<=calc_set
            %check for outliers before variance calculated
            for i=1:layers
                if abs(TAz(i)) > starter_threshold_T
%                     if (display_turb_tt); disp(['TAz replaced: ',num2str(TAz(i)),' @ Timestep:', num2str(tt)]); end
%                     TAz(i) = TAz_h(i,n-1); %Replace with previous timestep
                    TAz = TAz_h(:,n-1);
                    i=layers;
                    outlier = outlier + 1; %Mere Count the number of diverging timesteps
                end
            end
            for i=1:layers
                if CAz(i) > starter_threshold_C
                    if (display_turb_tt); disp(['CAz replaced: ',num2str(CAz(i)),' @ Timestep:', num2str(tt)]); end
                    CAz = CAz_h(:,n-1); %Replace with previous timestep
                    i=layers;
                    outlier = outlier + 1; %Mere Count the number of diverging timesteps
                end
            end
            for i=1:layers
                if EAz(i) > starter_threshold_E
                    if (display_turb_tt); disp(['EAz replaced: ',num2str(EAz(i)),' @ Timestep:', num2str(tt)]); end
                    EAz = EAz_h(:,n-1); %Replace with previous timestep
                    i=layers;
                    outlier = outlier + 1; %Mere Count the number of diverging timesteps
                end
            end
            
            sigma2_TAz = (n-2)./(n-1).*sigma2_TAz + (1./(n-1).*(TAz - ((n-1)/n).*(xbar_TAz + 1/(n-1).*TAz)).^2); %variance function of previous
            sigma2_CAz = (n-2)./(n-1).*sigma2_CAz + (1./(n-1).*(CAz - ((n-1)/n).*(xbar_CAz + 1/(n-1).*CAz)).^2);
            sigma2_EAz = (n-2)./(n-1).*sigma2_EAz + (1./(n-1).*(EAz - ((n-1)/n).*(xbar_EAz + 1/(n-1).*EAz)).^2);
            xbar_TAz = (n-1)/n.*(xbar_TAz + 1/(n-1).*TAz); %mean as a function of previous mean
            xbar_CAz = (n-1)/n.*(xbar_CAz + 1/(n-1).*CAz);
            xbar_EAz = (n-1)/n.*(xbar_EAz + 1/(n-1).*EAz);    
        end
        
        % Remove Outliers
        if n>calc_set
            
            sigma_TAz = sqrt(sigma2_TAz); %standard deviation
            sigma_CAz = sqrt(sigma2_CAz);
            sigma_EAz = sqrt(sigma2_EAz);

            for i=1:layers
                % Outlier if > or < 99.7 percentile
               if TAz(i) > xbar_TAz + num_sd*sigma_TAz(i)
                  if (display_turb_tt); disp(['TAz replaced: ',num2str(TAz(i)),' @ Timestep:', num2str(tt)]); end
                  TAz(i) = TAz_h(i,n-1); % replace with previous timestep values
                  outlier = outlier + 1; %Mere Count the number of outliers removed
               end
               if TAz(i) < xbar_TAz - num_sd*sigma_TAz(i)
                  if (display_turb_tt); disp(['TAz replaced: ',num2str(TAz(i)),' @ Timestep:', num2str(tt)]); end
                  TAz(i) = TAz_h(i,n-1);
                  outlier = outlier + 1; 
               end
               
               if CAz(i) > xbar_CAz + num_sd*sigma_CAz(i) 
                  if (display_turb_tt); disp(['CAz replaced: ',num2str(CAz(i)),' @ Timestep:', num2str(tt)]); end
                  CAz(i) = CAz_h(i,n-1);
                  outlier = outlier + 1; 
               end
               if CAz(i) < xbar_CAz - num_sd*sigma_CAz(i) 
                  if (display_turb_tt); disp(['CAz replaced: ',num2str(CAz(i)),' @ Timestep:', num2str(tt)]); end
                  CAz(i) = CAz_h(i,n-1);
                  outlier = outlier + 1; 
               end
               
               if EAz(i) > xbar_EAz + num_sd*sigma_EAz(i) 
                  if (display_turb_tt); disp(['EAz replaced: ',num2str(EAz(i)),' @ Timestep:', num2str(tt)]); end
                  EAz(i) = EAz_h(i,n-1); 
                  outlier = outlier + 1; 
               end
               if EAz(i) < xbar_EAz - num_sd*sigma_EAz(i) 
                  if (display_turb_tt); disp(['EAz replaced: ',num2str(EAz(i)),' @ Timestep:', num2str(tt)]); end
                  EAz(i) = EAz_h(i,n-1); 
                  outlier = outlier + 1; 
               end
            end
            
            %Calculate variance and mean for next timestep
            for i=1:layers
                sumt=0;
                sumc=0;
                sume=0;
                dift=0;
                difc=0;
                dife=0;
                for h=n-set_i:n-1
                    sumt=TAz_h(i,h)+sumt;
                    sumc=CAz_h(i,h)+sumc;
                    sume=EAz_h(i,h)+sume;
                end
                xbar_TAz(i) = (sumt+TAz(i))/calc_set; %mean
                xbar_CAz(i) = (sumc+CAz(i))/calc_set;
                xbar_EAz(i) = (sume+EAz(i))/calc_set;
                for h=n-set_i:n-1
                    dift=dift+(TAz_h(i,h)-xbar_TAz(i))^2;
                    difc=difc+(CAz_h(i,h)-xbar_CAz(i))^2;
                    dife=dife+(EAz_h(i,h)-xbar_EAz(i))^2;
                end
                dift=dift+(TAz(i)-xbar_TAz(i))^2;
                difc=difc+(CAz(i)-xbar_CAz(i))^2;
                dife=dife+(EAz(i)-xbar_EAz(i))^2;
                sigma2_TAz(i) = dift/set_i; %variance
                sigma2_CAz(i) = difc/set_i;
                sigma2_EAz(i) = dife/set_i;
            end
        end
        
        % ASSIGN
        VARIABLES.CANOPY.xbar_TAz = xbar_TAz;   %Mere Need these to carry over for mean calcs
        VARIABLES.CANOPY.xbar_CAz = xbar_CAz;   
        VARIABLES.CANOPY.xbar_EAz = xbar_EAz;   
        VARIABLES.CANOPY.sigma2_TAz = sigma2_TAz;   %Mere Need these to carry over for variance calcs
        VARIABLES.CANOPY.sigma2_CAz = sigma2_CAz;   
        VARIABLES.CANOPY.sigma2_EAz = sigma2_EAz;
        VARIABLES.CANOPY.TAz_history(:,tt) = TAz; %Mere Need these to carry over to replace if needed
        VARIABLES.CANOPY.CAz_history(:,tt) = CAz; 
        VARIABLES.CANOPY.EAz_history(:,tt) = EAz; 

 end
        
        