function [K_out,Theta_out, H_n1m, iter_out, q_pond_sat, q_infil, Mbe] = SOILMOISTURE_3D(Topo_mat, ...
    topo_sinx, topo_cosx, topo_siny, topo_cosy, Qt_mat, pond_depth_out, q_ppt, H_n1m, Theta_n1m, ...
    alpha, theta_S, theta_R, n, m, Ksat, k_ratio, Ss, ...
    nx, ny, nz, dx, dy, dzs, Dzs_mat, dtsub, t, t_scale, iter_param, stop_tol, max_iter, BCond)
%=========================================================================
% This function solves the 3D subsurface flow model implicitly using AIADI
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%     Topo_mat        % [m] 3D matrix of normalized topography
%     topo_sinx, topo_cosx, topo_siny, topo_cosy   % [] slopes of topography
%     Qt_mat          % [m/hr] root uptake flux 3D matrix
%     pond_depth_out  % [m] ponded water
%     q_ppt           % [m/hr] Infiltration into soil top layer
%     H_n1m           % [m] soil total head 3D matrix
%     Theta_n1m       % [m] soil water content 3D matrix
%     alpha, theta_S, theta_R, n, m, Ksat  % vanGenuchten parameters
%     k_ratio         % [] scaling parameter for anistropic Ksat
%     Ss              % [1/m] Specific storage parameter 
%     iter_param      % [] interation parameter for AIADI numerical scheme
%     BCond           % [] boundary condition
%------------------------- Output Variables -------------------------------
%     K_out           % [m/hr] hydraulic conductivity
%     Theta_out       % [-] soil water content 3D matrix
%     H_n1m           % [m] soil total head 3D matrix
%     iter_out        % [] number of iterations to convergence
%     q_pond_sat      % [m/hr] saturation excess as flux
%     q_infil         % [m/hr] infiltration flux
%     Mbe             % [m^3] mass balance error
%=========================================================================

    iter_out = 0; 
    q_pond_sat = zeros(nx,ny);
    Mbe = zeros(nx,ny,nz);
  
    % limit infiltration
    storage_avail_flux = sum((theta_S - Theta_n1m(:,:,:)).*Dzs_mat, 3) / dtsub;
    top_head = q_ppt*dtsub + pond_depth_out+Topo_mat(:,:,1)+dzs(1);
    
    for tsub=1:t_scale
        err = 1;
        niter_sub = 0;
        H_n = H_n1m;
        Theta_n = Theta_n1m;
        H_n1m1 = zeros(nx,ny,nz);
        H_n1m2 = zeros(nx,ny,nz);
        H_n1m3 = zeros(nx,ny,nz);
        
        while (err > stop_tol && niter_sub<max_iter)
            % calc C & K
            [C_n1m,K_n1m,Theta_n1m] = vanGenuchten(H_n1m-Topo_mat, alpha, theta_S, theta_R, n, m, Ksat);
            
            % get special K for top of domain
            Ktop = (K_n1m(:,:,1)+Ksat(:,:,1))/2;

            % calc actual infiltration
            q_capa = Ktop.*(top_head - H_n1m(:,:,1))./dzs(1);   % flux by head bc
            q_avail = (top_head-Topo_mat(:,:,1)-dzs(1))/dtsub;   % flux due to available water 
            top_flux = min(q_capa, q_avail);
            top_flux = min(top_flux, storage_avail_flux);
            head_boundary = q_capa == top_flux; %ones(size(q_capa)); %

            % get mid-point K for x,y,z-direction
            Kx_mid = (K_n1m(1:nx-1,:,:)+K_n1m(2:nx,:,:))./2;    
            Ky_mid = (K_n1m(:,1:ny-1,:)+K_n1m(:,2:ny,:))./2;    
            Kz_mid = (K_n1m(:,:,1:nz-1)+K_n1m(:,:,2:nz))./2;   
            
            
            % get full K for each of the 6 directions
            Kn = cat(1, K_n1m(1,:,:), Kx_mid);
            Ks = k_ratio * cat(1, Kx_mid, K_n1m(nx,:,:));
            Kw = k_ratio * cat(2, K_n1m(:,1,:), Ky_mid);
            Ke = k_ratio * cat(2, Ky_mid, K_n1m(:,ny,:));
            Ku = k_ratio * cat(3, Ktop, Kz_mid);
            Kd = cat(3, Kz_mid, K_n1m(:,:,nz));

            K_bar = Kn+Ks+Kw+Ke+Ku+Kd;  % (An 2011)
            Im = iter_param^niter_sub;  % Iteration parameter, (An 2011; Weeks 2004)

            % get coefficients for iterations
            % 1st pass
            C_ImK = C_n1m/dtsub + Im*K_bar;
            ThSS = Theta_n1m.*Ss./theta_S./dtsub;
            T_T = (Theta_n1m - Theta_n)/dtsub;
            Kndx = Kn/(dx^2);
            Ksdx = Ks/(dx^2);
            Kedy = Ke/(dy^2);
            Kwdy = Kw/(dy^2);
            Kudz = Ku./(Dzs_mat.^2);
            Kddz = Kd./(Dzs_mat.^2);

            % start solving the equation
            % 1 pass iterate over horizontal columns (x direction) because
            % matlab is column major
            parfor j=1:ny
                for k=1:nz

                    % get needed time diff h vectors x2
                    h_m_ijk = H_n1m(:,j,k);
                    h_n_ijk = H_n(:,j,k);

                    % get space diff h vectors x4 w/ boundary conditions
                    if j==ny
                        h_m_ijp1k = H_n1m(:,j,k);   % side boundaries, no flow
                    else
                        h_m_ijp1k = H_n1m(:,j+1,k);
                    end

                    if j==1
                        h_m_ijm1k = H_n1m(:,j,k);
                    else
                        h_m_ijm1k = H_n1m(:,j-1,k);
                    end

                    if k==nz
                        h_m_ijkp1 = H_n1m(:,j,k);  % bottom boundary, no flow
                        if BCond(6)
                            h_m_ijkp1 = H_n1m(:,j,k)-dzs(nz);  % bottom boundary, free flow
                        end
                    else
                        h_m_ijkp1 = H_n1m(:,j,k+1);
                    end

                    if k==1
                          h_m_ijkm1 = H_n1m(:,j,k) + top_flux(:,j)*dzs(1)./Ku(:,j,k);
                          h_m_ijkm1(head_boundary(:,j)) = top_head(head_boundary(:,j));
                    else
                        h_m_ijkm1 = H_n1m(:,j,k-1);
                    end

                    % get coeffs
                    % lhs
                    la = C_ImK(:,j,k) + ThSS(:,j,k) + Kndx(:,j,k).*topo_cosx(1:nx,j) + Ksdx(:,j,k).*topo_cosx(2:nx+1,j);
                    lb = Ksdx(:,j,k).*topo_cosx(2:nx+1,j);
                    lc = Kndx(:,j,k).*topo_cosx(1:nx,j);

                    %rhs
                    ra = C_ImK(:,j,k) - (Kedy(:,j,k).*topo_cosy(:,j+1) + Kwdy(:,j,k).*topo_cosy(:,j)) - (Kudz(:,j,k) + Kddz(:,j,k));
                    rb = ThSS(:,j,k);
                    rc = T_T(:,j,k);
                    rd = Kedy(:,j,k).*topo_cosy(:,j+1);
                    re = Kwdy(:,j,k).*topo_cosy(:,j);
                    rf = Kddz(:,j,k);
                    rg = Kudz(:,j,k);
                    rsin = (Ks(:,j,k).*topo_sinx(2:nx+1,j) - Kn(:,j,k).*topo_sinx(1:nx,j))/dx...
                        + (Ke(:,j,k).*topo_siny(:,j+1) - Kw(:,j,k).*topo_siny(:,j))/dy;

                    % set up linear eq and solve
                    % set up A matrix
                    % join diagonoals into matrix
                    la(1) = la(1) - lc(1);
                    la(nx) = la(nx) - lb(nx);
                    b_sp = [0; lb(1:end-1)];    % spdiags will cut off first value of b
                    c_sp = [lc(2:end); 0];   % spdiags will cut off last value of c
                    A = spdiags([-c_sp, la, -b_sp],[-1,0,1],nx,nx);

                    % rhs
                    rhs = ra.*h_m_ijk + rb.*h_n_ijk - rc + ...
                            rd.*h_m_ijp1k + re.*h_m_ijm1k + ...
                            rf.*h_m_ijkp1 + rg.*h_m_ijkm1 + ...
                            rsin - ...
                            Qt_mat(:,j,k);

                    % solve equation
                    H_n1m1(:,j,k) = A\rhs;
                end  
            end

            % 2nd pass
            parfor i=1:nx
                for k=1:nz
                    % get needed h vectors
                    h_m1 = H_n1m1(i,:,k)';
                    h_m = H_n1m(i,:,k)';

                    % get coeffs
                    rb = Kedy(i,:,k).*topo_cosy(i,2:ny+1);
                    rc = Kwdy(i,:,k).*topo_cosy(i,1:ny);
                    ra = rb+rc;
                    rd = C_ImK(i,:,k) + ThSS(i,:,k);

                    la = rd+ra;
                    lb = rb;
                    lc = rc;
                    
                    % set up linear eq and solve
                    % lhs - set up A matrix
                    % join diagonoals into matrix
                    la(1) = la(1) - lc(1);
                    la(ny) = la(ny) - lb(ny);
                    b_sp = [0; lb(1:end-1)'];    % spdiags will cut off first value of b
                    c_sp = [lc(2:end)'; 0];   % spdiags will cut off last value of c

                    A = spdiags([-c_sp, la', -b_sp],[-1,0,1],ny,ny);

                    % rhs - set up B matrix
                    ra(1) = ra(1) - rc(1); % no flow boundary condition
                    ra(ny) = ra(ny) - rb(ny);
                    B = spdiags([-c_sp, ra', -b_sp],[-1,0,1],ny,ny);

                    % get rhs
                    rhs_y = rd'.*h_m1 + B*h_m;

                    % solve equation
                    H_n1m2(i,:,k) = A\rhs_y;
                end
            end

            % 3rd pass
            parfor j=1:ny
                for i=1:nx
                    % get needed h vectors
                    h_m2 = squeeze(H_n1m2(i,j,:));
                    h_m = squeeze(H_n1m(i,j,:));

                    % get coeffs
                    rb = squeeze(Kddz(i,j,:));
                    rc = squeeze(Kudz(i,j,:));
                    ra = rb+rc;
                    rd = squeeze(C_ImK(i,j,:) + ThSS(i,j,:));

                    la = rd+ra;
                    lb = rb;
                    lc = rc;

                    % set up linear eq and solve
                    % lhs - set up A matrix
                    % join diagonoals into matrix
                    la(1) = la(1)-lc(1);
                    la(nz) = la(nz) - lb(nz);
                    b_sp = [0; lb(1:end-1)];    % spdiags will cut off first value of b
                    c_sp = [lc(2:end); 0];   % spdiags will cut off last value of c

                    A = spdiags([-c_sp, la, -b_sp],[-1,0,1],nz,nz);
                    rhs_bc = zeros(size(h_m));  
                    % consideration of bounary conditions are 0 because they cancel out on lhs and rhs
                    

                    % rhs - set up B matrix
                    ra(1) = ra(1) - rc(1);
                    ra(nz) = ra(nz) - rb(nz);
                    B = spdiags([-c_sp, ra, -b_sp],[-1,0,1],nz,nz);

                    % get rhs
                    rhs = rd.*h_m2 + B*h_m + rhs_bc;

                    % solve equation
                    H_n1m3(i,j,:) = A\rhs;
                end
            end

            % check error
            delta = squeeze(mean((H_n1m3-H_n1m),[1,2]));
            err = max(abs(delta(1:nz)),[],'all');

            % update H and index m
            if(err > 2)
                New_weight = 0.1;
            else
                New_weight = 1;
            end
            H_n1m = (1-New_weight) * H_n1m + New_weight * H_n1m3;
            niter_sub = niter_sub + 1;

            if niter_sub==max_iter
                fprintf('Reached Max Iter in Richards Solution at tstep=%d\n', t)
            end
        end
        
        [C_out,K_out,Theta_out] = vanGenuchten(H_n1m-Topo_mat, alpha, theta_S, theta_R, n, m, Ksat);
        
        iter_out = iter_out + niter_sub;
        
        %%  check mass balance ##############
        qdflux = zeros(nx,ny,nz);   
        quflux = zeros(nx,ny,nz);
        qnflux = zeros(nx,ny,nz);
        qsflux = zeros(nx,ny,nz);
        qwflux = zeros(nx,ny,nz);
        qeflux = zeros(nx,ny,nz);
        
        % make all q fluxes into system
        qdflux(:,:,1:nz-1) = ((H_n1m(:,:,2:nz) - H_n1m(:,:,1:nz-1))./Dzs_mat(:,:,1:nz-1)).* Kd(:,:,1:nz-1);   % entering system

        if BCond(6)
            qdflux(:,:,nz) = -1 * Kd(:,:,nz);   % set only if free flow boundary condition; else no flow.
        end
        
        quflux(:,:,2:nz) = ((H_n1m(:,:,1:nz-1) - H_n1m(:,:,2:nz))./Dzs_mat(:,:,2:nz)).* Ku(:,:,2:nz);
        q_infil = top_flux; % ((H_n1m(:,:,1) - H_n1m(:,:,2))/dz).* Ktop;
        q_head_bc = ((top_head - H_n1m(:,:,1))/dzs(1)).* Ktop;
        q_infil(head_boundary) = q_head_bc(head_boundary);
        quflux(:,:,1) = q_infil; %((H_n1m(:,:,1) - H_n1m(:,:,2))/dz).* Ktop;
        
        
        
        qnflux(2:nx,:,:) = ((H_n1m(1:nx-1,:,:) - H_n1m(2:nx,:,:))/dx).* Kn(2:nx,:,:);   % 0 at boundaries
        qsflux(1:nx-1,:,:) = ((H_n1m(2:nx,:,:) - H_n1m(1:nx-1,:,:))/dx).* Ks(2:nx,:,:); 
        qwflux(:,2:ny,:) = ((H_n1m(:,1:ny-1,:) - H_n1m(:,2:ny,:))/dy).* Kw(:,2:ny,:);
        qeflux(:,1:ny-1,:) = ((H_n1m(:,2:ny,:) - H_n1m(:,1:ny-1,:))/dy).* Ke(:,1:ny-1,:);
        
        ssflux = (H_n1m - H_n).*(Ss./dtsub.*Theta_n1m./theta_S).*Dzs_mat;  % leaving system
        qtflux = Qt_mat;    % leaving system

        Mass_in = ((qdflux + quflux - qtflux - ssflux)*dx*dy + (qnflux + qsflux)*dy.*Dzs_mat + (qwflux + qeflux)*dx.*Dzs_mat) * dtsub; %[m^3]
        
        % add supersaturated amount to ponding input flux
        InfilAvail = (theta_S - Theta_n) * dx * dy .* Dzs_mat;  % [m^3]
        SuperSatAmount = Mass_in - InfilAvail; 
        SuperSatAmount(SuperSatAmount < 0) = 0;
        q_pond_sat = q_pond_sat + sum(SuperSatAmount, 3)/dx/dy/(dtsub*t_scale); %[m/hr]
        
        % correct mass balance
        dtheta = (Theta_out - Theta_n) * dx * dy .* Dzs_mat;  %additional mass in [m^3]
        Mbe_sub = dtheta - Mass_in;
        Mbe = Mbe + Mbe_sub;
        
    end
    iter_out = iter_out / t_scale;
    
end