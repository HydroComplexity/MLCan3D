% Kunxuan wang
function [rpp,rpp_wgt,krad,kax,dwat,smp,kboundary,klayer,...
    qlayer,layeruptake,layeruptake_all,mberrormm, type, hor_drainage,hor_drainage_lay,flux_Ss]=...
    ROOTSOIL_EXPLICIT(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, CONSTANTS, nspecies)

%=========================================================================
% This code solves the explicit function for water flow in the plant root 
% system and soil moisture where root pressure and conductivity are not 
% iterated with soil moisture model. Root and soil model based on 
% ROOTSOIL, not ROOT_HR and ROOT_NOHR
%
%------------------------- Input Variables -------------------------------
%       SWITCHES        % SWITCHES structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%       VARIABLES       % VARIABLES structure
%       CONSTANTS       % CONSTANTS structure
%       nspecies        % Number of species that are considered in the simulations
%------------------------- Output Variables ------------------------------
%       rpp             % [mm] root water potential 
%       rpp_wgt         % [mm] Wighted mean rpp over root uptake profile 
%       krad            % [/s] radial conductivity of the root 
%       kax             % [mm/s] axial conductivity of the root
%       dwat            % [%] Change in soil moisture
%       smp             % [mm] soil matric potential
%       kboundary       % [mm/s] Soil Hydraulic Conductivity at Interface of layers
%       klayer          % [mm/s] Soil Hydraulic Conductivity at Node Layers
%       qlayer          % [mm/s] Fluxes of water between layers in the soil 
%       layeruptake     % [mm/s] Fluxes of plant water uptake at each layer. Shows the total from all the species        
%       layeruptake_all % [mm/s] Fluxes of plant water uptake at each layer. Shows the fluxes for each sps separately        
%       mberrormm       % [mm/dtime] Mass balance error.
%       type            % 1.OK, 2. Supersaturation First Layer, 3. Super Saturation in between layers, thus hor_drainage >0 
%       hor_drainage    % [mm/s] Horizontal Drainage. Only happens when superaturation occurs.
%       hor_drainage_lay% [mm/s] Horizontal drainage if super saturation for all layers
% 
%========================================================================

    % De-reference Structure Values
    dtime=CONSTANTS.dtime;
    zmm=VERTSTRUC.znsmm;
    dzmm=VERTSTRUC.dzsmm;
    znsmm=VERTSTRUC.zhsmm;    
    zimm=VERTSTRUC.zhsmm;
    nl_soil = PARAMS.Soil.nl_soil;     

    n=VERTSTRUC.VanGen_n;
    alpha=VERTSTRUC.VanGen_alpha;
    theta_dry=VERTSTRUC.theta_dry;
    vanGen=SWITCHES.vanGen;
    thetar=VERTSTRUC.VanGen_thetar;
            
    for ii=1:1:nspecies
        if VARIABLES.CANOPY.TR_can_all(ii) > 0
            etr(ii) =  VARIABLES.CANOPY.TR_can_all(ii);
        else
            etr(ii) =  0;            
        end        
        rootfr(:,ii) = VERTSTRUC.rootfr(:,ii); 
        roottr(:,ii) = VERTSTRUC.roottr(:,ii);
        nl_root(ii) = VERTSTRUC.nl_root(ii);
        rpp(:,ii) = VARIABLES.ROOT.rpp(:,ii);
        K_rad(ii) = PARAMS.Soil.K_rad(ii);
        K_axs(ii) = PARAMS.Soil.K_axs(ii);                
    end
    smp = VARIABLES.SOIL.smp;    
    wliq = VARIABLES.SOIL.volliq;
    effporsl = VERTSTRUC.eff_poros;    
    hr= SWITCHES.HR_on;    
    rhc = SWITCHES.rhc;
    rtcond = SWITCHES.rtcond;    
    phi0 = VERTSTRUC.psi0;
    bsw=VERTSTRUC.bsw;
    hksati = VERTSTRUC.HKsat;
    pthr = VARIABLES.SOIL.qinfl;
    PLC = VARIABLES.ROOT.PLC;
    kpar_ax = PARAMS.Soil.kpar_ax;
    plants = SWITCHES.plants;     

    %========================================================================
    % Set type to [1, 0  0 ]. Type is a vector saying the type of numerical
    % colution for the root-soil solution
    type = [1 0 0];

    kpar_ax_spc = kpar_ax(1);
    for ii=1:1:nspecies
    [krad(:,ii),kax(:,ii)] = conductivities(nl_soil,zmm,dzmm, znsmm,rootfr(:,ii),...
                                roottr(:,ii),rpp(:,ii),smp,wliq,effporsl,...
                                   K_rad(ii),K_axs(ii),1,rhc, etr(ii), PLC(:,ii),kpar_ax_spc,rtcond); 

    end

    if plants == 1
        for ii=1:1:nspecies
            kpar_ax_spc = kpar_ax(ii);

            [rpp(:,ii)] = rootmodel(nl_soil,nl_root(ii),zmm,etr(ii),smp,krad(:,ii),kax(:,ii));


            [krad(:,ii),kax(:,ii),etr(ii)] = conductivities(nl_soil,zmm,dzmm,znsmm,rootfr(:,ii),...
                                roottr(:,ii),rpp(:,ii),smp,wliq,effporsl,...
                                   K_rad(ii),K_axs(ii),hr(ii),rhc,etr(ii), PLC(:,ii),kpar_ax_spc,rtcond); 

        end
    else 

       krad = zeros(nl_soil,nspecies);

    end

    [dwat,smp,kboundary,klayer,ft,fb,qlayer,layeruptake,layeruptake_all,mberrormm,type,hor_drainage,hor_drainage_lay,flux_Ss] ...
        = soilmodel(nl_soil,dtime,effporsl,phi0,bsw,hksati,zmm',dzmm',zimm',wliq,plants,rpp,krad,...
        pthr,nspecies, type, n, alpha, thetar,vanGen);

    if plants == 1        
        for ii=1:1:nspecies     
           [krad(:,ii),kax(:,ii),etr(ii)] = conductivities(nl_soil,zmm,dzmm,znsmm, rootfr(:,ii),...
                       roottr(:,ii),rpp(:,ii),smp,wliq,effporsl,...
                       K_rad(ii),K_axs(ii),hr(ii),rhc,etr(ii), PLC(:,ii),kpar_ax_spc,rtcond); 
        end                                         

    else
        % Dongkook: Start
        %krad = zeros(nl_soil,1);
        krad = zeros(nl_soil,nspecies);
        % Dongkook: End
    end

    for ii=1:1:nspecies
        rpp_wgt(ii) = sum(rpp(:,ii).*roottr(:,ii))/sum(roottr(:,ii));% 
    end
end