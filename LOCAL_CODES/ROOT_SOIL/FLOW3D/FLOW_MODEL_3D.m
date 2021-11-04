%% 3D soil moisture using mixed richards eq with AIADI
function [VARIABLES, dwat,smp,kboundary,klayer,qlayer,layeruptake,layeruptake_all, ...
    mberrormm, type, hor_drainage,hor_drainage_lay,flux_Ss] = ...
    FLOW_MODEL_3D(VERTSTRUC, PARAMS, HORSTRUC, CONSTANTS, VARIABLES, krad, rpp, type, nspecies)
%=========================================================================
% This function solves the 3D soil moisture model and 2D overland flow model
% using implicit schemes for a given time step. This function is a wrapper 
% functions that calls the 3D soil moisture module and the 2D overland flow 
% module, facilitate connection between the two processes, and aggregate
% results to 1D to fit the existing MLCan framework.
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       VERTSTRUC       % VERTSTRUC structure
%       HORSTRUC        % HORSTRUC structure
%       PARAMS          % PARAMS structure
%       VARIABLES       % VARIABLES structure
%       CONSTANTS       % CONSTANTS structure
%       FORCING         % CONSTANTS structure
%       SWITCHES        % SWITCHES structure
%       krad            % [1/s] Root Radial conductivities
%       rpp             % [mm] Root water potential 
%       type            % 1.OK, 2. Supersaturation First Layer, 3. Super Saturation in between layers, thus hor_drainage >0 
%------------------------- Output Variables -------------------------------
%       VARIABLES       % VARIABLES structure
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
%       flux_Ss         % [mm/s] Fluxes of elastic storage for mass balance
%=========================================================================
    %% De-reference Structure Values

    % model grid variables
    dzs = VERTSTRUC.dzs;    % vertical vector of soil layer thickness [m]
    Dzs_mat = PARAMS.SOIL3D.Dzs_mat;    % matrix version of dzs
    Topo_mat = PARAMS.SOIL3D.Topo_mat;
    nl_soil = PARAMS.Soil.nl_soil;  % number of soil layers
    BCond = PARAMS.SOIL3D.BCond;
    % dem = HORSTRUC.dem;
    topo = HORSTRUC.topo_norm;
    nx = HORSTRUC.nx;
    ny = HORSTRUC.ny;
    dx = HORSTRUC.dx;  % [m] horizontal x direction resolution 
    dy = HORSTRUC.dy;  % [m] horizontal y direction resolution 
    topo_sinx = HORSTRUC.topo_sinx; % sin and co of ropography slopes
    topo_cosx = HORSTRUC.topo_cosx;
    topo_siny = HORSTRUC.topo_siny;
    topo_cosy = HORSTRUC.topo_cosy;

    x_ind = 1:nx; y_ind = 1:ny; z_ind = 1:nl_soil; 

    % time
    dtime = CONSTANTS.dtime;  
    dt = dtime/3600; % [hr] timestep converted from seconds
    timestep = VARIABLES.timestep;  % [-] timestep index or count
    t_scale = PARAMS.SOIL3D.t_scale; % sub-timescale scaling factor
    dtsub = dt/t_scale;     % sub-timescale times step

    % tolerances and numerical parameters
    stop_tol = PARAMS.SOIL3D.stop_tol;
    max_iter = PARAMS.SOIL3D.max_iter; 
    iter_param = PARAMS.SOIL3D.iter_param;   % iteration parameter for ADADI
    New_weight = PARAMS.SOIL3D.New_weight;
    %Ksat_weight = 0.5;
    save_interval = PARAMS.SOIL3D.save_interval;
    output_folder = PARAMS.SOIL3D.output_folder;

    % Soil Parameters
    alpha   = PARAMS.SOIL3D.alpha;       % parameter related to the inverse of the air entry suction [1/cm]
    n       = PARAMS.SOIL3D.n;           % Pore-size distributions [-]
    m       = PARAMS.SOIL3D.m;           % [-]
    theta_R = PARAMS.SOIL3D.theta_R;     % Residual water content [-]
    theta_S = PARAMS.SOIL3D.theta_S;     % Saturated water content [-] 
    Ksat    = PARAMS.SOIL3D.Ksat;        % Saturated hydraulic conductivity [m/hr]
    k_ratio  = PARAMS.SOIL3D.k_ratio;      % Scalar scale factor for horizontal anisotrophic hydraulic conductivity []
    Ss      = PARAMS.SOIL3D.Ss;          % Storage 
    % Psimin = 0.0005;
    % air_dry = -15.0;
    
    % Soil model
    K_n1m = VARIABLES.SOIL3D.K_n1m;
    Theta_n1m = VARIABLES.SOIL3D.Theta_n1m;
    H_n1m = VARIABLES.SOIL3D.H_n1m;

    % Overland 
    s_min = PARAMS.OVERLAND2D.s_min;
    h_min = PARAMS.OVERLAND2D.h_min;
    mann = PARAMS.OVERLAND2D.mann;
    pond_depth = VARIABLES.OVERLAND2D.pond_depth; % need to maintain separate from qinfl so that there is overland flow

    %% Get forcings for this timestep
    % get root uptake flux
    Krad_mat = ScatterCol2Mat(zeros(nx,ny,nl_soil), sum(krad, 2)*3600, x_ind, y_ind, z_ind); % from [1/s] to [1/hr] sum on all species
    Rpp_mat = ScatterCol2Mat(zeros(nx,ny,nl_soil), sum(rpp, 2)/1000, x_ind, y_ind, z_ind);  % from [mm] to [m] sum on all species 
    Krad_mat(Theta_n1m <= theta_R) = 0;
    Qroot_mat = Krad_mat.*(H_n1m-Rpp_mat);  % [m/hr] root uptake flux 3D matrix based on root Krad, positive leaving soil
    
    % get surface infiltration flux
    qinfl = VARIABLES.SOIL.qinfl;  % [mm/s] Net Infiltration into the soil, includs Esoil
    qinfl_sub = qinfl * 3600 / 1000 * ones(nx,ny); % from [mm/s] to [m/hr] Infiltration into soil top layer 2D. qin does not /t_scale because it is a flux

    %% Subsurface and Overland model 
    % qinfl_sub_sum = sum(qinfl_sub, 'all') 
    [K_out, Theta_out, H_out, Iter_out, q_pond_sat, q_soil, Mbe] = ...
    SOILMOISTURE_3D(Topo_mat, topo_sinx, topo_cosx, topo_siny, topo_cosy, ...
        Qroot_mat, pond_depth, qinfl_sub, H_n1m, Theta_n1m,...
        alpha, theta_S, theta_R, n, m, Ksat, k_ratio, Ss, ...
        nx, ny, nl_soil, dx, dy, dzs, Dzs_mat, dtsub, timestep, t_scale, iter_param, stop_tol, max_iter, BCond);

    % set type var for supersaturated conditions
    if(sum(pond_depth, 'all') > 1e-4)
        type(2) = 1;
    end
    if(sum(q_pond_sat, 'all') > 1e-4)
        type(3) = 1;
    end
  
    pond_depth = pond_depth + q_pond_sat*dt; % [m] ponded water due to supersaturation
    q_pond = pond_depth/dt;
    
    pond_depth_out = ...
    OVERLAND_2D(topo, pond_depth, qinfl_sub, q_soil, mann, ...
            s_min, h_min, nx, ny, dx, dy, dtsub, timestep, [0,0,0,0]);

    % Update#
    VARIABLES.SOIL3D.K_n1m = K_out;
    VARIABLES.SOIL3D.Theta_n1m = Theta_out;
    VARIABLES.SOIL3D.H_n1m = H_out;
    VARIABLES.OVERLAND2D.pond_depth = pond_depth_out;
    VARIABLES.OVERLAND2D.q_pond = q_pond; % flux into overland flow

    %% End 3D part and aggregate for 1D outputs
    % aggregate to 1D
    dzsmm = VERTSTRUC.dzsmm;          % [mm] Vector with Layer thicknesses
    znsmm = VERTSTRUC.znsmm;           % [mm] Vector with depth of the different nodes in the layers  
    smp = GatherMat2Col(H_out - Topo_mat, BCond)*1000; % [mm] Soil matrix solution  
    smp_init = GatherMat2Col(H_n1m - Topo_mat, BCond)*1000; % [mm] Soil matrix solution  
    klayer = GatherMat2Col(K_out, BCond)*1000/3600;               % [mm/s] Hydraulic conductivity at middle of each soil layers
    kboundary = comave(klayer,znsmm,dzsmm,VERTSTRUC.zhsmm,nl_soil);          % [mm/s] Hydraulic conductivity at interface between soil layers, length n-1
    fb = klayer(nl_soil);           % [mm/s] Flux boundary condition at the botoom.
    theta_init = GatherMat2Col(Theta_n1m, BCond);
    theta_final = GatherMat2Col(Theta_out, BCond);
    theta_S_1D = VERTSTRUC.eff_poros;
    CC_Ss = diag(Ss./1000./theta_S_1D.*theta_final./dtime); % [-/s] diagonal matrix for elastic storage term     

    % type            % 1.OK, 2. Supersaturation First Layer, 3. Super Saturation in between layers, thus hor_drainage >0 

    [qlayer,layeruptake,layeruptake_all,fluxt,fluxb,dwat,flux_Ss,mberror, mberrormm] = ...
    compflux(smp, rpp, kboundary, krad, qinfl, fb, dzsmm, znsmm, nl_soil, type, theta_init, theta_final, CC_Ss, smp_init, dtime, nspecies);

    hor_drainage = GatherMat2Col(q_pond, BCond);   
    hor_drainage_lay = zeros(nl_soil,1);

    %% output files
    if mod(timestep, save_interval) == 0
        % save subsurface
        filename = sprintf('subsurface3D_%d.mat',timestep);
        save(fullfile(pwd,output_folder,filename), 'H_out', 'Theta_out', 'K_out', 'Mbe');

        % save overland
        filename = sprintf('overland2D_%d.mat',timestep);
        save(fullfile(pwd,output_folder,filename), 'q_pond_sat', 'pond_depth_out');
    end

end
% 
% %% util funtions
% % 
% function mat = ScatterLayer2Mat(mat, layer, z_ind, increment)
%     for k=z_ind
%         mat(:,:,k) = layer+(k-1)*increment;
%     end
% end
% 
% function vec = MatrixVerticalStdv(mat)
%     vec = squeeze(std(mat,0,[1,2]));
% end
% 
% function topo = GetTopo(dem)
%     topo = dem - min(dem,[],'all');
% end
% 
% function PlotSurface(surface, nx, ny)
%     [X,Y] = meshgrid(1:nx,1:ny);
%     surf(X,Y,surface,'FaceAlpha',0.5)
% end
