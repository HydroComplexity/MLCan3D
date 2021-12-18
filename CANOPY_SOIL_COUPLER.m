nl_soil=PARAMS.nl_soil;

load './Temps/temp_variable.mat'...
    'dat_root1' 'dat_root2' 'dat_root3' 'dat_root4'
zns = dat_root1(:,1);

% Soil layer thicknesses
dzs(1)  = 0.5*(zns(1)+zns(2));
dzs(nl_soil)= zns(nl_soil)-zns(nl_soil-1);
for j = 2:nl_soil-1
    dzs(j)= 0.5*(zns(j+1)-zns(j-1));
end
dzs=dzs';

% Soil layer interface depths from the surface [m]
zhs(nl_soil) = zns(nl_soil) + 0.5*dzs(nl_soil);
for j = 1:nl_soil-1
    zhs(j)= 0.5*(zns(j)+zns(j+1));
end
zhs=zhs';

znsmm = zns(:)*1000;      % [mm]
dzsmm = dzs(:)*1000;      % [mm]
zhsmm = zhs(:)*1000;      % [mm]

if (kspecies == 5 && Sim_species_con == 1 )
    rootfr = dat_root1(:,2);
    roottr = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root = sum(count_nl_root1 ~= 0);
elseif (kspecies == 5 && Sim_species_con == 2 )
    rootfr = dat_root2(:,2);
    roottr = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root = sum(count_nl_root2 ~= 0);
elseif (kspecies == 5 && Sim_species_con == 3 )
    rootfr = dat_root3(:,2);
    roottr = dat_root3(:,2);
    count_nl_root3=dat_root3(:,2);
    count_nl_root3(1)=1;
    nl_root = sum(count_nl_root3 ~= 0);
elseif (kspecies == 5 && Sim_species_con == 4 )
    rootfr = dat_root4(:,2);
    roottr = dat_root4(:,2);
    count_nl_root4=dat_root4(:,2);
    count_nl_root4(1)=1;
    nl_root = sum(count_nl_root4 ~= 0);
elseif (kspecies == 1)
    rootfr = dat_root1(:,2);
    roottr = rootfr;
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root = sum(count_nl_root1 ~= 0);
elseif (kspecies == 2)
    rootfr(:,1) = dat_root1(:,2);
    roottr(:,1) = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root(1) = sum(count_nl_root1 ~= 0);
    rootfr(:,2) = dat_root2(:,2);
    roottr(:,2) = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root(2) = sum(count_nl_root2 ~= 0);
elseif (kspecies == 3)
    rootfr(:,1) = dat_root1(:,2);
    roottr(:,1) = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root(1) = sum(count_nl_root1 ~= 0);
    rootfr(:,2) = dat_root2(:,2);
    roottr(:,2) = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root(2) = sum(count_nl_root2 ~= 0);
    rootfr(:,3) = dat_root3(:,2);
    roottr(:,3) = dat_root3(:,2);
    count_nl_root3=dat_root3(:,2);
    count_nl_root3(1)=1;
    nl_root(3) = sum(count_nl_root3 ~= 0);
elseif (kspecies == 4)
    rootfr(:,1) = dat_root1(:,2);
    roottr(:,1) = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root(1) = sum(count_nl_root1 ~= 0);
    rootfr(:,2) = dat_root2(:,2);
    roottr(:,2) = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root(2) = sum(count_nl_root2 ~= 0);
    rootfr(:,3) = dat_root3(:,2);
    roottr(:,3) = dat_root3(:,2);
    count_nl_root3=dat_root3(:,2);
    count_nl_root3(1)=1;
    nl_root(3) = sum(count_nl_root3 ~= 0);
    rootfr(:,4) = dat_root4(:,2);
    roottr(:,4) = dat_root4(:,2);
    count_nl_root4=dat_root4(:,2);
    count_nl_root4(1)=1;
    nl_root(4) = sum(count_nl_root4 ~= 0);
end

% ASSIGN
VERTSTRUC.zns = zns;
VERTSTRUC.dzs = dzs;
VERTSTRUC.zhs = zhs;
VERTSTRUC.znsmm = znsmm;
VERTSTRUC.dzsmm = dzsmm;
VERTSTRUC.zhsmm = zhsmm;
VERTSTRUC.rootfr = rootfr;
VERTSTRUC.roottr = roottr;
VERTSTRUC.nl_root = nl_root;
VERTSTRUC.nl_soil = nl_soil;
PARAMS.Soil.nl_soil = nl_soil;

[VERTSTRUC] = SOIL_PROPERTIES(PARAMS, VERTSTRUC);
theta_dry = VERTSTRUC.theta_dry;
porsl = VERTSTRUC.porsl;
TK_dry = VERTSTRUC.TK_dry;
TK_sol = VERTSTRUC.TK_sol;
HC_sol = VERTSTRUC.HC_sol;

% ALLOCATE STORAGE FOR MODELLED VARIABLES
ALLOCATE_STORAGE;

% INITIALIZE CANOPY STATES
VARIABLES.CANOPY.gsv_sun = 0.01*ones(nl_can,nspecies);
VARIABLES.CANOPY.gsv_shade = 0.01*ones(nl_can,nspecies);
VARIABLES.CANOPY.TR = zeros(length(znc),nspecies);
VARIABLES.CANOPY.Sh2o_prof = zeros(length(znc),1);
VARIABLES.CANOPY.Tl_prev_dt = Ta_in(1) * ones(nl_can,1);

% INITIALIZE SOIL STATES
snowIC = 0;         %Mere Assume zero snow initial conditions [0]; Manually input snow Initial conditions [1]
if (~snowIC)
    % Initialize snow moisture variables
    VARIABLES.SOIL.voltotsn = 1;
    VARIABLES.SOIL.voltotli = 1;
    VARIABLES.SOIL.volliqli = volliqliinit;
    VARIABLES.SOIL.voliceli = 0;
    VARIABLES.SOIL.volliqsn = 0;
    VARIABLES.SOIL.volicesn = 0;
    VARIABLES.SOIL.zliqsl = (VARIABLES.SOIL.dzlit_m*1000)*volliqliinit;
    VARIABLES.SOIL.zicesl = 0;
    VARIABLES.SOIL.zliqsl_prev = VARIABLES.SOIL.zliqsl;
    VARIABLES.SOIL.zicesl_prev = VARIABLES.SOIL.zicesl;
    VARIABLES.SOIL.zsn = VARIABLES.SOIL.zliqsl;
    VARIABLES.SOIL.wliqsl = (VARIABLES.SOIL.zliqsl/1000)*PARAMS.Soil.rho_liq;
    VARIABLES.SOIL.wicesl = 0;
    VARIABLES.SOIL.wsn = VARIABLES.SOIL.wliqsl;
    VARIABLES.SOIL.rhosn = 1000;
else
    % Document where this data is from: WCr 2003 spinup 2 multispecies output
    % Initialize snow moisture variables manually
    VARIABLES.SOIL.voltotsn = [0.697256924620478];
    VARIABLES.SOIL.voltotli = [1.15953221134674];
    VARIABLES.SOIL.volliqli = [2.32452945780892e-17];
    VARIABLES.SOIL.voliceli = [1.15953221134674];
    VARIABLES.SOIL.volliqsn = [1.39780011722062e-17];
    VARIABLES.SOIL.volicesn = [0.697256924620478];
    VARIABLES.SOIL.zliqsl = [6.97358837342677e-16];
    VARIABLES.SOIL.zicesl = [34.7859663404021];
    VARIABLES.SOIL.zliqsl_prev = VARIABLES.SOIL.zliqsl;
    VARIABLES.SOIL.zicesl_prev = VARIABLES.SOIL.zicesl;
    VARIABLES.SOIL.zsn = [49.8897395093442];
    VARIABLES.SOIL.wliqsl = [6.97358837342677e-16];
    VARIABLES.SOIL.wicesl = [31.8987311341487];
    VARIABLES.SOIL.wsn = [31.8987311341487];
    VARIABLES.SOIL.rhosn = [639.384599876978];

%     % Document where this data is from: WCr 2004 multispecies output
%     % Initialize snow moisture variables manually
%     VARIABLES.SOIL.voltotsn = [0.728301208399523];
%     VARIABLES.SOIL.voltotli = [2.23389090121064];
%     VARIABLES.SOIL.volliqli = 0;
%     VARIABLES.SOIL.voliceli = [2.23389090121064];
%     VARIABLES.SOIL.volliqsn = 0;
%     VARIABLES.SOIL.volicesn = [0.728301208399523];
%     VARIABLES.SOIL.zliqsl = 0;
%     VARIABLES.SOIL.zicesl = [67.0167270363192];
%     VARIABLES.SOIL.zliqsl_prev = VARIABLES.SOIL.zliqsl;
%     VARIABLES.SOIL.zicesl_prev = VARIABLES.SOIL.zicesl;
%     VARIABLES.SOIL.zsn = [92.0178715391557];
%     VARIABLES.SOIL.wliqsl = 0;
%     VARIABLES.SOIL.wicesl = [61.4543386923047];
%     VARIABLES.SOIL.wsn = [61.4543386923047];
%     VARIABLES.SOIL.rhosn = [667.852208102363];
end

% Fixing the constant soil layer problem
VARIABLES.SOIL.volice = zeros(nl_soil,1);
VARIABLES.SOIL.snow_tcount = 0;
VARIABLES.SOIL.volliqli = volliqliinit ;   % Initial value of litter soil moisture


VARIABLES.CANOPY.Sh2o_can_prev = 0;

% initialize ice content & soil temperature
Ts = Tsinit;
VARIABLES.SOIL.Tli = Ta_in(1);
VARIABLES.SOIL.Tsl = Tslint;
VARIABLES.SOIL.Tlprev = Ta_in(1);
VARIABLES.SOIL.TKsoil=VERTSTRUC.TK_sol;

%% INITIALIZE 3D FLOW MODEL STATES
if SWITCHES.soil3D
    % TOPO SLOPE CALCULATIONS
    % dem, nx, ny etc loaded in LOAD_SITE_INFO.m
    topo = double(dem - mean(dem(:)));    % Get normalized topography from DEM layer

    % calculate sin and cos of topography slopes
    % get height difference between cells
    topo_diff_x = zeros(nx+1,ny);
    topo_diff_y = zeros(nx,ny+1);
    topo_diff_x(2:nx,:) = topo(2:nx,:) - topo(1:nx-1,:);
    topo_diff_y(:,2:ny) = topo(:, 2:ny) - topo(:, 1:ny-1);
    % get hypotenues 
    hypo_x = sqrt(topo_diff_x.^2 + dx^2);
    hypo_y = sqrt(topo_diff_y.^2 + dy^2);
    % get slope
    topo_sinx = topo_diff_x./hypo_x;     % opposite over hypo
    topo_cosx = dx./hypo_x;  % adjacent over hypo
    topo_siny = topo_diff_y./hypo_y;     
    topo_cosy = dy./hypo_y;  
    
    HORSTRUC.topo_norm = topo;
    HORSTRUC.topo_sinx = topo_sinx;
    HORSTRUC.topo_cosx = topo_cosx;
    HORSTRUC.topo_siny = topo_siny;
    HORSTRUC.topo_cosy = topo_cosy;
    
    % 3D FLOW PARAMETERS
    load './Temps/temp_variable.mat' 'para_soil3D';
    if any(strcmp(para_soil3D(:,2), '2D/3D matrix'))
        load './Temps/temp_variable.mat' 'working_params3D'
        load(working_params3D);
    end
    
    row_num = 1; % alpha: parameter related to the inverse of the air entry suction [1/cm]
    var_choice = para_soil3D{row_num,2};
    if strcmp(var_choice,'Texture based')
        % use texture based values
        % set 3D values using 1D
        alpha = ScatterCol2Mat(zeros(nx,ny,nl_soil), VERTSTRUC.VanGen_alpha * 10, 1:nx, 1:ny, 1:nl_soil);        % [1/mm] to [1/cm]
    elseif strcmp(var_choice,'Single value')
        % use single value
        % take default from para_soil3D
        if isempty(para_soil3D{row_num,3})
            alpha = para_soil3D{row_num,4};
        else
            alpha = para_soil3D{row_num,3};
        end
        VERTSTRUC.VanGen_alpha(:) = alpha/10; % [1/cm] to [1/mm]
    else
        % else use user input 3D values, already loaded
        VERTSTRUC.VanGen_alpha(:) = MatrixVerticalAvg(alpha)/10; % [1/cm] to [1/mm]
    end
    PARAMS.SOIL3D.alpha = alpha;
    
    row_num = 2; % n: Pore-size distributions [-] 
    var_choice = para_soil3D{row_num,2};
    if strcmp(var_choice,'Texture based')
        n = ScatterCol2Mat(zeros(nx,ny,nl_soil), VERTSTRUC.VanGen_n, 1:nx, 1:ny, 1:nl_soil);
    elseif strcmp(var_choice,'Single value')
        if isempty(para_soil3D{row_num,3})
            n = para_soil3D{row_num,4};
        else
            n = para_soil3D{row_num,3};
        end
        VERTSTRUC.VanGen_n(:) = n;
    else
        % else use user input 3D values, already loaded
        VERTSTRUC.VanGen_n = MatrixVerticalAvg(n);        
    end
    PARAMS.SOIL3D.n = n;
    m       = (n-1)./n;              % [-]
    PARAMS.SOIL3D.m = m;
    
    row_num = 3; % theta_S = 0.4: Saturated water content [-] 
    var_choice = para_soil3D{row_num,2};
    if strcmp(var_choice,'Texture based')
        theta_S = ScatterCol2Mat(zeros(nx,ny,nl_soil), VERTSTRUC.porsl, 1:nx, 1:ny, 1:nl_soil);
    elseif strcmp(var_choice,'Single value')
        if isempty(para_soil3D{row_num,3})
            theta_S = para_soil3D{row_num,4};
        else
            theta_S = para_soil3D{row_num,3};
        end
        VERTSTRUC.eff_poros(:) = theta_S;
    else
        % else use user input 3D values, already loaded
        VERTSTRUC.eff_poros = MatrixVerticalAvg(theta_S); 
    end
    VERTSTRUC.porsl = VERTSTRUC.eff_poros;
    PARAMS.SOIL3D.theta_S = theta_S;
    
    row_num = 4; % theta_R = 0.05:  Residual water content [-]
    var_choice = para_soil3D{row_num,2};
    if strcmp(var_choice,'Texture based')
        theta_R = ScatterCol2Mat(zeros(nx,ny,nl_soil), VERTSTRUC.VanGen_thetar, 1:nx, 1:ny, 1:nl_soil);
    elseif strcmp(var_choice,'Single value')
        if isempty(para_soil3D{row_num,3})
            theta_R = para_soil3D{row_num,4};
        else
            theta_R = para_soil3D{row_num,3};
        end
        VERTSTRUC.theta_dry(:) = theta_R;
        VERTSTRUC.VanGen_thetar(:) = theta_R;
    else
        % else use user input 3D values, already loaded
        VERTSTRUC.theta_dry = MatrixVerticalAvg(theta_R); 
        VERTSTRUC.VanGen_thetar(:) = MatrixVerticalAvg(theta_R);
    end
    PARAMS.SOIL3D.theta_R = theta_R;
    
    row_num = 5; % Ksat = 0.03: Saturated hydraulic conductivity [m/hr] %
    var_choice = para_soil3D{row_num,2};
    if strcmp(var_choice,'Texture based')
        Ksat = ScatterCol2Mat(zeros(nx,ny,nl_soil), VERTSTRUC.HKsat*3600/1000 , 1:nx, 1:ny, 1:nl_soil); % [mm/s] to [m/hr]
    elseif strcmp(var_choice,'Single value')
        if isempty(para_soil3D{row_num,3})
            Ksat = para_soil3D{row_num,4};
        else
            Ksat = para_soil3D{row_num,3};
        end
        VERTSTRUC.HKsat(:) = Ksat/3600*1000;  % [m/hr] to [mm/s]
    else
        % else use user input 3D values, already loaded
        VERTSTRUC.HKsat = MatrixVerticalAvg(Ksat)/3600*1000; % [m/hr] to [mm/s]
    end
    PARAMS.SOIL3D.Ksat = Ksat;

    % 2D Overland flow
    row_num = 6; % mann = 0.025*((12*2.54/100)^(1/3))/3600;   
    % units hr/m^(1/3) Need to do unit conversion from look up tables (sec/ft^(1/3))!!!
    var_choice = para_soil3D{row_num,2};
    if strcmp(var_choice,'Texture based')
        % Soil texture based mannings constant
        if SWITCHES.litter
            mann = 0.025*((12*2.54/100)^(1/3))/3600; % for earth and veg resistance
        else 
            mann = 0.018*((12*2.54/100)^(1/3))/3600; % for bare earth
        end
        mann = mann * ones(nx,ny);
    elseif strcmp(var_choice,'Single value')
        if isempty(para_soil3D{row_num,3})
            mann = para_soil3D{row_num,4};
        else
            mann = para_soil3D{row_num,3};
        end
        mann = mann * ones(nx,ny);
    end
    PARAMS.OVERLAND2D.mann =  mann; 
    
    % INITIAL 3D CONDITIONS AND STATES
    Topo_mat = zeros(nx, ny, nl_soil);    % [m] topography based gravimetric potential matrix
    for k=1:nl_soil
        Topo_mat(:,:,k) = topo-zns(k);
    end
    PARAMS.SOIL3D.Topo_mat = Topo_mat;
    PARAMS.SOIL3D.Dzs_mat = ScatterCol2Mat(zeros(nx,ny,nl_soil), dzs, 1:nx, 1:ny, 1:nl_soil); % dz in 3D

    % 3D INITIAL CONDITIONS
    load './Temps/temp_variable.mat' 'Sub3DFlow_3d_pond_depth_init' 'Sub3DFlow_3d_Theta_init';
    if Sub3DFlow_3d_pond_depth_init
        load './Temps/temp_variable.mat' 'pond_depth_init';
    else
        pond_depth_init = zeros(nx,ny);
    end 
    VARIABLES.OVERLAND2D.pond_depth = double(pond_depth_init);
    VARIABLES.OVERLAND2D.q_pond = zeros(nx,ny); % flux into overland flow

    if Sub3DFlow_3d_Theta_init
        load './Temps/temp_variable.mat' 'Theta_init';
        % check if greater than 3D thetaR or porosl !!! TODO
        if any(Theta_init > theta_S)
            msgbox('Initial soil moisture is higher than saturated soil moisture!', 'MLCan Error','error');
            error('Initial soil moisture is higher than saturated soil moisture!')
        end
        if any(Theta_init <= theta_R)
            msgbox('Initial soil moisture is lower than residual soil moisture!','MLCan Error','error');
            error('Initial soil moisture is lower than residual soil moisture!')
        end
    else
        Theta_init = ScatterCol2Mat(zeros(nx,ny,nl_soil), volliqinit, 1:nx, 1:ny, 1:PARAMS.nl_soil); % don't need to error check , checked later
    end 

    Smp = vanGenuchtenInverse(Theta_init, alpha, theta_S, theta_R, n, m);  %matric potential !!!!! need to set vg params 
    [C_n1m,K_n1m,Theta_n1m] = vanGenuchten(Smp, alpha, theta_S, theta_R, n, m, Ksat);
    H_n1m = Smp+Topo_mat;
    
    VARIABLES.SOIL3D.K_n1m = K_n1m;
    VARIABLES.SOIL3D.Theta_n1m = Theta_n1m;
    VARIABLES.SOIL3D.H_n1m = H_n1m;

    if abs(sum(Theta_init-Theta_n1m, 'all')) > 1e-10  % test if matrices are approximately equal
        msgbox({'VanGenuchten forward and reverse functions do not results in same soil moisture', 'Solution: Check vanGenuchtenInverse() and vanGenuchten() '},'error');
        error('VanGenuchten forward and reverse functions wrong')
    end
    
    % 3D VanGen parameters
    load './Temps/temp_variable.mat' 'para_soil3D';
end
    
% initialize soil moisture
theta_dry = VERTSTRUC.theta_dry;
porsl = VERTSTRUC.porsl;
VARIABLES.SOIL.volliq = volliqinit;
VARIABLES.SOIL.smp = VERTSTRUC.psi0 .* (VARIABLES.SOIL.volliq ./ VERTSTRUC.porsl).^(-VERTSTRUC.bsw);

% Message box for soil moisture
if sum(volliqinit>porsl) > 0 %Mere; changed >= to >
    msgbox({'Initial soil moisture is higher than saturated soil moisture!', 'Solution: Modify initilal soil moisture or Increase % of sand.'},'error');
    error('Initial soil moisture is higher than saturated soil moisture!')
end
if sum(volliqinit<=theta_dry) > 0
    msgbox({'Initial soil moisture is lower than residual soil moisture!', 'Solution: Modify initilal soil moisture or Increase % of sand.'},'error');
    error('Initial soil moisture is lower than residual soil moisture!')
end

%% 
% INITIALIZE ROOT POTENTIAL
for ii=1:nspecies
    VARIABLES.ROOT.rpp_wgt(:,ii) =  VARIABLES.SOIL.smp(1);
    VARIABLES.ROOT.rpp(:,ii)= VARIABLES.SOIL.smp;
end

% PEDOTRANSFER FUNCTIONS
if SWITCHES.Pedofunctions % 0 by default
    [VERTSTRUC] = PEDOSOIL_PROPERTIES(PARAMS, VERTSTRUC, VARIABLES);
    porsl = VERTSTRUC.porsl;                                         % POROSITY
    psi0 = VERTSTRUC.psi0;                                           % MINIMUM SOIL SUCTION = SOIL POTENTIAL AT SATURATION [mm]
    bsw = VERTSTRUC.bsw;                                             % B PARAMETER BROKS AND COREY SHAPE PARAMETER
    Ksat = VERTSTRUC.HKsat;                                          % HYDRAULIC CONDUCTIVITY AT SATURATION [mm / s]
    eff_poros = VERTSTRUC.eff_poros;
end


% RUN CANOPY-ROOT-SOIL MODEL
% LOOP OVER EACH YEAR TO RE-INITIALIZE CANOPY/SOIL STATES FOR EACH YEAR

% TimeBar 1/3
tff=yendinds(end);
t00=ybeginds(1);
hh = timebar('Progress','MLCan Simulation');

tic
for yy = 1:length(Each_year)
    % compute the range of time steps in current year
    yy;
    ybind = ybeginds(yy);
    yeind = yendinds(yy);
    
    % LOOP OVER EACH TIME PERIOD IN YEAR yy
    %for tt = ybind:yeind
    for tt = ybind:1:yeind
        
        tt;
        % TimeBar 2/3
        timebar(hh,(tt-t00)/(tff-t00))

        timestep = tt-ybind + 1;
        VARIABLES.timestep = timestep;
        [VERTSTRUC VARIABLES rootfr] = ROOT_RESPONSE_DRY(VARIABLES,...
            SWITCHES, VERTSTRUC, CONSTANTS, PARAMS, doy, smp_store);
        
        % FORCING CONDITIONS
        FORCING.doy = doy(tt);
        FORCING.Rg = Rg_in(tt);
        FORCING.Pa = Pa_in;
        if PARAMS.LWcom == 1
            FORCING.LWdn = LWdn_in(tt);
        end
        FORCING.zen = ZEN_in(tt);
        FORCING.U = U_in(tt);
        FORCING.ppt = PPT_in(tt);    % [mm]
        FORCING.Ta = Ta_in(tt);
        FORCING.ea = ea_in(tt);
        FORCING.Ca = CO2base;
        FORCING.ELEV=ELEV;
        FORCING.Ustar = ustar_in(tt); %Mere
        
        if (~SWITCHES.soilheat_on)
            VARIABLES.SOIL.Ts = (Ta_in(tt)-5)*ones(size(zns));
        else
            VARIABLES.SOIL.Ts = Ts;
        end
        VARIABLES.SOIL.Ts = Ts;
        VARIABLES.SOIL.Tsurf=Ts(1);
        
        % CANOPY STRUCTURE
        if SWITCHES.plants
            for kk=1:1:nspecies
                if SWITCHES.LT == 1
                    LAILT = LAI_in(timestep,kk);
                    VERTSTRUC.LAIzall(:,kk) = LAILT*LADnorm_all(:,kk);
                else
                    VERTSTRUC.LAIzall(:,kk) = LAI_in(tt,kk)*LADnorm_all(:,kk);
                end
            end
        else
            VERTSTRUC.LAIzall(:,1:nspecies) = zeros(nl_can,nspecies);
        end
        LADnorm = sum(VERTSTRUC.LAIzall,2)/sum(LAI_in(tt,:));
        VERTSTRUC.LAIz = sum(VERTSTRUC.LAIzall,2);
        VERTSTRUC.LADz = VERTSTRUC.LAIz ./ dzc; % Total LAD distribution
        fLAIz =VERTSTRUC.LAIzall./(repmat(sum(VERTSTRUC.LAIzall,2),1,nspecies));
        fLAIz(isnan(fLAIz)) = 0; % Set to zero whenever there is not LAI at a given layer
        VERTSTRUC.fLAIz = fLAIz; % Fraction of LAI in each species at each relative height level

        % create vinds
        % 1. For the total canopy
        LADmax = (max(VERTSTRUC.LAIzall,[],2)); % Maximum LAD
        nvinds = find(LADmax<=0);
        vinds = find(LADmax>0);
        VERTSTRUC.vinds = vinds;
        VERTSTRUC.nvinds = nvinds;
        
        % 2. For All the species
        for kk=1:nspecies
            nvinds_all{kk} = find(VERTSTRUC.LAIzall(:,kk) <= 0);
            vinds_all{kk} = find(VERTSTRUC.LAIzall(:,kk) > 0);
        end
        VERTSTRUC.nvinds_all = nvinds_all;
        VERTSTRUC.vinds_all = vinds_all;
        
        % INITIALIZE CANOPY ENVIRONMENT
        VARIABLES.CANOPY.TAz = Ta_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.CAz = CO2base * ones(nl_can,1);
        VARIABLES.CANOPY.EAz = ea_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.PAz = Pa_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.Uz = U_in(tt) * ones(nl_can,1);
        
        VARIABLES.CANOPY.TR_sun = zeros(nl_can,nspecies);
        VARIABLES.CANOPY.TR_shade = zeros(nl_can,nspecies);
        
        % INITIALIZE CANOPY STATES
        VARIABLES.CANOPY.Tl_can_sun = VARIABLES.CANOPY.TAz;
        VARIABLES.CANOPY.Tl_can_shade = VARIABLES.CANOPY.TAz;
        VARIABLES.CANOPY.Tl_sun = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
        VARIABLES.CANOPY.Tl_shade = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
        VARIABLES.CANOPY.Ci_sun = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        VARIABLES.CANOPY.Ci_shade = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        
       if tt == 1 %Mere added loop so there is a variable to pass in CANOPY_MODEL
           VARIABLES_last_tt = VARIABLES;
           countingmaxedout=0; 
           diverging=0;
           largeremainder=0;
           outlier=0;
       end
       repeat_noturb = 0;
        
        % CANOPY MODEL SOLUTION      
        [An_can, Ph_can, LE_can, H_can, dHcan, Rnrad_can, TR_can, ...
            Fc_soil, LE_soil, H_soil, Rnrad_soil, G, Tsurf, remainsoil,remaincan,remaineco,...
            Rnrad_sun, Rnrad_shade, Rnrad_eco, ...
            An_sun, An_shade, LE_sun, LE_shade, H_sun, H_shade, TR_sun, TR_shade, ...
            Tl_sun, Tl_shade, psil_sun, psil_shade, gsv_sun, gsv_shade, fsvg_sun, fsvm_sun,...
            fsvg_shade, fsvm_shade, Ci_sun, Ci_shade, CAz, TAz, EAz, Uz, gbv_sun, gbh_sun, gbv_shade, gbh_shade, ...
            LAI_sun, LAI_shade, fsun, fshade, ...
            Ph_limit_sun, Jc_C3_sun, Jj_C3_sun, Js_C3_sun, Jc_C4_sun, Jj_C4_sun, Js_C4_sun, ...
            Ph_limit_shade, Jc_C3_shade, Jj_C3_shade, Js_C3_shade, Jc_C4_shade, Jj_C4_shade, Js_C4_shade, ...
            PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, SWout, ...
            LWabs_can, LWabs_soil, LWemit_soil, LWemit_can, LWemit_sun, LWemit_shade, LWout, LWoutM, RH_soil, fdiff, ...
            Sh2o_prof, Sh2o_can, ppt_ground, Ch2o_prof, Ch2o_can, Evap_prof, Evap_can, ...
            dryfrac, wetfrac, Vz, VARIABLES, FORCING,...
            SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out, SWsoildir_in, SWsoildir_out,...
            SWsoildif_in, SWsoildif_out, LWabs_canM, LWabs_soilM, LSshaCON, LSsunCON, countingmaxedout, diverging, largeremainder, outlier, repeat_noturb] = ...
            CANOPY_MODEL(SWITCHES, VERTSTRUC, FORCING, PARAMS, VARIABLES, CONSTANTS, tt, VARIABLES_last_tt, countingmaxedout, diverging, largeremainder, outlier, repeat_noturb); %Mere added tt
        
        if (repeat_noturb)
            % INITIALIZE CANOPY ENVIRONMENT
            VARIABLES.CANOPY.TAz = Ta_in(tt) * ones(nl_can,1);
            VARIABLES.CANOPY.CAz = CO2base * ones(nl_can,1);
            VARIABLES.CANOPY.EAz = ea_in(tt) * ones(nl_can,1);
            VARIABLES.CANOPY.PAz = Pa_in(tt) * ones(nl_can,1);
            VARIABLES.CANOPY.Uz = U_in(tt) * ones(nl_can,1);

            VARIABLES.CANOPY.TR_sun = zeros(nl_can,nspecies);
            VARIABLES.CANOPY.TR_shade = zeros(nl_can,nspecies);

            % INITIALIZE CANOPY STATES
            VARIABLES.CANOPY.Tl_can_sun = VARIABLES.CANOPY.TAz;
            VARIABLES.CANOPY.Tl_can_shade = VARIABLES.CANOPY.TAz;
            VARIABLES.CANOPY.Tl_sun = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
            VARIABLES.CANOPY.Tl_shade = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
            VARIABLES.CANOPY.Ci_sun = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
            VARIABLES.CANOPY.Ci_shade = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        
            [An_can, Ph_can, LE_can, H_can, dHcan, Rnrad_can, TR_can, ...
            Fc_soil, LE_soil, H_soil, Rnrad_soil, G, Tsurf, remainsoil,remaincan,remaineco,...
            Rnrad_sun, Rnrad_shade, Rnrad_eco, ...
            An_sun, An_shade, LE_sun, LE_shade, H_sun, H_shade, TR_sun, TR_shade, ...
            Tl_sun, Tl_shade, psil_sun, psil_shade, gsv_sun, gsv_shade, fsvg_sun, fsvm_sun,...
            fsvg_shade, fsvm_shade, Ci_sun, Ci_shade, CAz, TAz, EAz, Uz, gbv_sun, gbh_sun, gbv_shade, gbh_shade, ...
            LAI_sun, LAI_shade, fsun, fshade, ...
            Ph_limit_sun, Jc_C3_sun, Jj_C3_sun, Js_C3_sun, Jc_C4_sun, Jj_C4_sun, Js_C4_sun, ...
            Ph_limit_shade, Jc_C3_shade, Jj_C3_shade, Js_C3_shade, Jc_C4_shade, Jj_C4_shade, Js_C4_shade, ...
            PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, SWout, ...
            LWabs_can, LWabs_soil, LWemit_soil, LWemit_can, LWemit_sun, LWemit_shade, LWout, LWoutM, RH_soil, fdiff, ...
            Sh2o_prof, Sh2o_can, ppt_ground, Ch2o_prof, Ch2o_can, Evap_prof, Evap_can, ...
            dryfrac, wetfrac, Vz, VARIABLES, FORCING,...
            SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out, SWsoildir_in, SWsoildir_out,...
            SWsoildif_in, SWsoildif_out, LWabs_canM, LWabs_soilM, LSshaCON, LSsunCON, countingmaxedout, diverging, largeremainder, outlier, repeat_noturb] = ...
            CANOPY_MODEL(SWITCHES, VERTSTRUC, FORCING, PARAMS, VARIABLES, CONSTANTS, tt, VARIABLES_last_tt, countingmaxedout, diverging, largeremainder, outlier, repeat_noturb); %Mere added tt
        end
        
        % SOLUTION OF SNOW-LITTER PACK DYNAMICS
        if SWITCHES.litter
            [VARIABLES] = FLUXES_WATER_SOIL_LITTER (PARAMS, VARIABLES, CONSTANTS, FORCING, SWITCHES);
        else
            [VARIABLES] = FLUXES_WATER_SOIL (PARAMS, VARIABLES, CONSTANTS, FORCING, SWITCHES);
        end
        
        % assign
        qinfl = VARIABLES.SOIL.qinfl;
        qinflL = VARIABLES.SOIL.qinfl;
        net_qinflL = VARIABLES.SOIL.net_qinflL;
        drainlitter = VARIABLES.SOIL.drainlitter;
        volliqli = VARIABLES.SOIL.volliqli;
        
        % Implicit Solution
        if (SWITCHES.soil3D)
            [VARIABLES, rpp,rpp_wgt,krad,kax,dwat,smp,bk,hk,...
                qlayer,layeruptake,layeruptake_all,mberrormm, type, hor_drainage,hor_drainage_lay,flux_Ss] = ...
                ROOTSOIL_3D(SWITCHES, VERTSTRUC, HORSTRUC, PARAMS, VARIABLES, CONSTANTS, nspecies);
%             disp('In ROOTSOIL_3D');
            VARIABLES.ROOT.rpp = rpp;
            VARIABLES.ROOT.rpp_wgt = rpp_wgt;
            VARIABLES.ROOT.krad = krad;
            VARIABLES.SOIL.type = type;
            VARIABLES.SOIL.flux_Ss = flux_Ss;
        elseif (SWITCHES.ns)
            if tt == 321
                stop = 1;
            end
            [rpp,rpp_wgt,krad,kax,dwat,smp,bk,hk, ...
                qlayer,layeruptake,layeruptake_all,mberrormm,type, hor_drainage, hor_drainage_lay,flux_Ss]=ROOTSOIL(SWITCHES, VERTSTRUC,...
                PARAMS, VARIABLES, CONSTANTS, nspecies);
            VARIABLES.ROOT.rpp = rpp;
            VARIABLES.ROOT.rpp_wgt = rpp_wgt;
            VARIABLES.ROOT.krad = krad;
            VARIABLES.SOIL.type = type;
            VARIABLES.SOIL.flux_Ss = flux_Ss;
        else
            [rpp,rpp_wgt,krad,kax,dwat,smp,bk,hk,qlayer,layeruptake, ...
                layeruptake_all,mberrormm,type, hor_drainage, hor_drainage_lay,flux_Ss]=ROOTSOIL_EXPLICIT(SWITCHES, VERTSTRUC,...
                PARAMS, VARIABLES, CONSTANTS, nspecies);
            
            VARIABLES.ROOT.rpp = rpp;
            VARIABLES.ROOT.rpp_wgt = rpp_wgt;
            VARIABLES.ROOT.krad = krad;
            VARIABLES.SOIL.type = type;
            VARIABLES.SOIL.flux_Ss = flux_Ss;


%             if (SWITCHES.HR_on)
%                 [rpp, rpp_wgt, krad, kax] = ROOTS_HR(SWITCHES, VERTSTRUC, PARAMS, VARIABLES);
%             else
%                 [rpp, rpp_wgt, krad, kax] = ROOTS_NOHR(SWITCHES, VERTSTRUC, PARAMS, VARIABLES);
%             end
%             
%             VARIABLES.ROOT.rpp = rpp;
%             VARIABLES.ROOT.rpp_wgt = rpp_wgt;
%             VARIABLES.ROOT.krad = krad;
%             
%             % Soil Moisture Solution
% %             [dwat, smp, hk, smp_wgt, thsatfrac_wgt, qlayer] = ...
% %                 SOILMOISTURE(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, CONSTANTS);
% 
%             % De-reference Structure Values
%             dtime=CONSTANTS.dtime;
%             zmm=VERTSTRUC.znsmm;
%             dzmm=VERTSTRUC.dzsmm;   
%             zimm=VERTSTRUC.zhsmm;
%             nl_soil = PARAMS.Soil.nl_soil;     
% 
%             % Dongkook: Start
%             n=VERTSTRUC.VanGen_n;
%             alpha=VERTSTRUC.VanGen_alpha;
%             vanGen=SWITCHES.vanGen;
%             thetar=VERTSTRUC.VanGen_thetar;
%             wliq = VARIABLES.SOIL.volliq;
%             effporsl = VERTSTRUC.eff_poros; 
%             phi0 = VERTSTRUC.psi0;
%             bsw=VERTSTRUC.bsw;
%             hksati = VERTSTRUC.HKsat;
%             pthr = VARIABLES.SOIL.qinfl;
%             plants = SWITCHES.plants;   
%             
%             %[dwat,smp,bk,hk,ft,fb,flux_s,flux_sr,flux_sr_all,mberrormm,type, hor_drainage, hor_drainage_lay,flux_Ss]
%             %[dwat,smp,bk,hk,ft,fb,qlayer,layeruptake,layeruptake_all,mberrormm,type,hor_drainage,hor_drainage_lay,flux_Ss]
%             [dwat,smp,bk,hk,ft,fb,qlayer,layeruptake,layeruptake_all,mberrormm,type,hor_drainage,hor_drainage_lay,flux_Ss] ...
%             = soilmodel(nl_soil,dtime,effporsl,phi0,bsw,hksati,zmm',dzmm',zimm',wliq,plants,rpp,krad,...
%             pthr,nspecies, type, n, alpha, thetar,vanGen);
%             
%             VARIABLES.SOIL.type = type;  % TODO: update #########
%             VARIABLES.SOIL.flux_Ss = flux_Ss;
%             hor_drainage = nan;
%             layeruptake = (smp - rpp).*krad;
%             layeruptake_all = layeruptake;
%             hor_drainage = 0;
%             hor_drainage_lay = zeros(length(smp),1);
%             mberrormm = nan;
        end
        
        % Update Volumetric Liquid Content
        volliq = VARIABLES.SOIL.volliq;
        volliq = volliq(:) + dwat(:);
        volliq = max(theta_dry, volliq);
        volliq = min(VERTSTRUC.eff_poros, volliq);
        
        % ASSIGN
        VARIABLES.SOIL.dwat = dwat;
        VARIABLES.SOIL.volliq = volliq;
        VARIABLES.SOIL.smp = smp;
        VARIABLES.SOIL.qlayer = qlayer;
        VARIABLES.SOIL.hor_drainage = hor_drainage;
        VARIABLES.SOIL.hor_drainage_lay = hor_drainage_lay;
        VARIABLES.SOIL.layeruptake = layeruptake;
        VARIABLES.SOIL.layeruptake_all = layeruptake_all;
        
        % RECOMPUTE MASS BALANCE INCLUDING THE FLUX BACK FROM
        % INFILTRATION SOLUTION
        if SWITCHES.litter
            [VARIABLES] = FLUXES_WATER_SOIL_LITTER_BACK (VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, SWITCHES);
        else
            [VARIABLES] = FLUXES_WATER_SOIL_BACK (VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);
        end
        % assign
        qinfl = VARIABLES.SOIL.qinfl;
        qinflL = VARIABLES.SOIL.qinfl;
        net_qinflL = VARIABLES.SOIL.net_qinflL;
        drainlitter = VARIABLES.SOIL.drainlitter;
        volliqli = VARIABLES.SOIL.volliqli;
        
        if (SWITCHES.soilheat_on)
            Ginto = VARIABLES.SOIL.Gs;
            
            % Soil Temperature Solution
            volice = 0;
            [VARIABLES] = SOILHEAT (Ginto, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);
            Ts = VARIABLES.SOIL.Ts;
            cpv = VARIABLES.SOIL.cpv;
        end
        
        % ************************************************************************
        %     [PARAMS, VARIABLES] = Nitrogen_Plant(PARAMS, FORCING, VARIABLES, CONSTANTS,nspecies);
        % ************************************************************************
        if SWITCHES.soilCN_on
            % For nitrogen remobilization
            if tt == 1
                STORAGE.UP_nit_m2_store = zeros(nl_soil,nl_soil,1);
                STORAGE.UP_amm_m2_store = zeros(nl_soil,nl_soil,1);
            end

            [VARIABLES, SWITCHES, PARAMS, FORCING] = ...
                core_N(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE);
            
            CN_STORE_DATA ();
        end

        if (SWITCHES.entropy_on)            
            [SSresults] = ...
                COMPUENTROPY (SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
                SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
                SWout, fdiff,LWabs_canM, LWabs_soilM, LWemit_soil, LWemit_sun, LWemit_shade,...
                LWout, Tsurf, FORCING,SWITCHES, CONSTANTS, PARAMS,...
                VARIABLES, VERTSTRUC);
        end
        
        [VARIABLES] = MASS_BALANCE (VARIABLES, CONSTANTS, PARAMS, FORCING, SWITCHES, VERTSTRUC, tt);
        
        STORE_DATA;
        
        VARIABLES_last_tt = VARIABLES; %Mere added
        
    end
    
    disp(['Total outliers in turbulence calculations: ', num2str(outlier)]);%Mere
    ratio_outliers = outlier/(yeind*3)/nl_can;
    disp(['Percent of timesteps with outliers removed: ', num2str(ratio_outliers)]);%Mere
    disp(['Total timesteps diverging: ', num2str(diverging)]);%Mere
    disp(['Total timesteps reaching max iterations: ', num2str(countingmaxedout)]);%Mere
    disp(['Total timesteps with large remainder: ', num2str(largeremainder)]);%Mere
    ratioreplaced = (diverging + countingmaxedout + largeremainder)/yeind;
    disp(['Percent of timesteps replaced: ', num2str(ratioreplaced)]);%Mere
    
    toc
end
delete(gcp('nocreate'));    % close parallel pool if it is active
% Timebar 3/3
close(hh);
