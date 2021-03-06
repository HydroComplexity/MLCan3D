
function [rpp, rpp_weight, krad, kax] = ROOTS_NOHR( SWITCHES, VERTSTRUC, PARAMS, VARIABLES )

%=========================================================================
% This code solves the model for water flow in the plant root system. 
% The upper boundary condition is set to the transpiration rate while 
% the lower boundary is set to no flux.
%
%------------------------- Input Variables -------------------------------
%       dtime           % time step [s]
%       z               % layer depth [mm]
%       dz              % layer thickness [mm]
%       zi              % interface level below a "z" level [mm]
%       TR             % actual transpiration [mm/s]
%       rootfr          % root fraction of a layer [-]
%       rpp             % root pressure potential [mm]
%       smp             % soil matric potential [mm]
%       volliq         % soil water per unit volume [mm/mm]
%       eff_porosity    % effective porosity of soil
%       thetadry        % soil moisture content at dryness [-]
%       K_rad           % radial conductivity of the root system [s^-1]
%       K_axs           % axial specific conductivity of the root system [mm/s]
%       hr              % option for hydraulic redistribution [-]
%       rhc             % option for root hydraulic conductivity [-]
%
%------------------------- Output Variables ------------------------------
%       rpp             % root pressure potential [mm]
%       krad            % radial conductivity of the root [/s]
%       kax             % axial conductivity of the root [mm/s]
%                    
%-------------------------- local variables ------------------------------
%       amx             % "a" left off diagonal of tridiagonal matrix
%       bmx             % "b" diagonal column for tridiagonal matrix
%       cmx             % "c" right off diagonal tridiagonal matrix
%       dmx             % "d" forcing term of tridiagonal matrix
% 
%=========================================================================

    
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    rhc = SWITCHES.rhc;

    TR = VARIABLES.CANOPY.TR_can;
    
    smp = VARIABLES.SOIL.smp;
    volliq = VARIABLES.SOIL.volliq;
    
    z = VERTSTRUC.znsmm;
    dz = VERTSTRUC.dzsmm;
    rootfr = VERTSTRUC.rootfr;
    thetadry = VERTSTRUC.theta_dry;
    eff_porosity = VERTSTRUC.eff_poros;

    K_rad = PARAMS.Soil.K_rad;
    K_axs = PARAMS.Soil.K_axs;
%*************************************************************************
%*************************************************************************



N = length(z);  % # soil layers

% Compute the radial and axial conductivities of the roots
%-------------------------------------------------------------------------
% Root conductativities (both radial and axial) for each soil layer are
% obtained by weighting the conductivity of the root system by the root
% distribution within the layer. The effect of soil moisture on root
% conductivity is also taken into account.

for i = 1:N

    % for radial conductivity, root fraction is used as a weighting factor,
    % because the uptake from a layer is proportional to the total surface 
    % area of roots in that layer. Thus, 
    krad(i) = rootfr(i)*K_rad ...
              * volliq(i)/eff_porosity(i);

    % for axial conductivity, root density is used as a weighting factor,
    % because the flow is proportional to the total x-sectional area of  
    % the roots in that layer. Thus,                       
    kax(i) = (rootfr(i)/dz(i))*K_axs ...
             * volliq(i)/eff_porosity(i);                            

end
krad = krad(:);
kax = kax(:);

% For the case where the root hydraulic conductivity is allowed to increase
% with depth, a linear increasing effect is considered
if rhc == 1     % if conductivity is allowed to increases with depth
    krad = krad + (z/1000).*krad;   % linearly increasing 
    kax  = kax + (z/1000).*kax;
end


% Root Pressure Potential --> No HR
    if rhc == 1     % if increasing root hydraulic conductivity is allowed
        for i = 1:N
            rootr(i) = max(0, rootfr(i)*((volliq(i)-thetadry(i))...
                /(eff_porosity(i)-thetadry(i)))); % effect of soil moisture
            % include the effect of root conductivity. Here, a linearly
            % increasing function is considered
            rootr(i) = rootr(i) + (z(i)/1000)*rootr(i);  % effect of root conductivity
        end
    else
        for i = 1:N
            rootr(i) = max(0, rootfr(i)*((volliq(i)-thetadry(i))...
                /(eff_porosity(i)-thetadry(i))));  % effect of soil moisture
        end
    end
    if (sum(rootr) ~= 0), rootr = rootr/sum(rootr); end
    
    for i = 1:N
        et(i) = TR * rootr(i);              % layer contribution to transpiration
        rpp(i) = smp(i) - et(i)/krad(i);    % root pressure potential [mm]
    end
   

% Weighted mean smp over root uptake profile [mm]
rpp = rpp(:);
rinds = find(rootfr>0);
rpp_weight = sum(rpp(rinds).*rootfr(rinds)/sum(rootfr(rinds)));% * mmH2OtoMPA;


