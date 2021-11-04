function [U, Km, tau] = ORDER_1_CLOSURE_U(VERTSTRUC, PARAMS, Ucan)

%=========================================================================
% Augments current canopy water storage with intercepted precipitation
%
% Written By: Darren Drewry
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
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
  
    LAD = VERTSTRUC.LADz;                                                  %[m^2 leaf / m^2 ground] LAI at different layers 
    znc = VERTSTRUC.znc;                                                   %[m] Height of canopy levels 
    dzc = VERTSTRUC.dzc;                                                   % [m] Tickness of layers in canopy 
    
    hcan = PARAMS.CanStruc.hcan;                                           %[m] Height canopy  
    Cd = PARAMS.MicroEnv.Cd;                                               %[] Constant for cons. of momentum model   
    alph = PARAMS.MicroEnv.alph;                                           %[] mixing length parameter, Constant for cons. of momentum model

%*************************************************************************
 
% Correct for when wind speed is too low for convergence (%Mere)
% min_wind = 0.6;
% if Utop < min_wind
%     Utop = min_wind;
% end

% Mixing Length
l_mix = alph * hcan;

z = [0, znc];
Ulow = 0;

N=length(z);
U = linspace (Ulow, Ucan, N);   %Utop --> Ucan
LAD = [0; LAD(:)];
LAD = LAD';

%-------- Start iterative solution
err = 10^9;
epsilon = 10^-2;
while err > epsilon

    %-------- dU/dz
    y=[];
    dU=[];
    y(2:N) = diff(U) ./ dzc;
    y(1) = y(2);
    dU = y;

    %------Add model for diffusivity (Km)
    Km=( (l_mix).^2) .* abs(y);

    tau=-Km.*y;

    %------Set up coefficients for ODE
    a1=-Km;
    a2(2:N) = -diff(Km) ./ dzc;
    a2(1) = a2(2);
    a3 = 0.5 * Cd * LAD .* U;
    dx=dzc;

    %------ Set the elements of the Tri-diagonal Matrix
    upd=(a1 ./ (dx * dx) + a2 ./ (2 * dx));
    dia=(-a1 .* 2 / (dx * dx) + a3);
    lod=(a1 ./ (dx * dx)-a2./(2 * dx));
    co=ones(1,N)*0;
    aa=[];
    bb=[];
    cc=[];
    dd=[];
    co(1)=Ulow;
    co(N)=Ucan;     %Utop --> Ucan
    aa=lod;
    bb=dia;
    cc=upd;
    dd=co;
    aa(1)=0;
    aa(N)=0;
    cc(1)=0;
    cc(N)=0;
    bb(1)=1;
    bb(N)=1;

    %------Use the Thomas Algorithm to solve the tridiagonal matrix
    %Un = Thomas(aa, bb, cc, dd);
    Un = TRIDIAG(N,aa,bb,cc,dd);
    err = max( abs(Un - U) );

    %-----Use successive relaxations in iterations
    eps1 = 0.5;
    U=(eps1 * Un + (1 - eps1) * U);
    
end

% Remove added bottom node for surface no-slip
U = U(2:end);
U = U(:);
