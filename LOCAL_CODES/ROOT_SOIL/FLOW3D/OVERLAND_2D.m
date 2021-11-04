function water_depth_new = OVERLAND_2D(topo, water_depth, ppt_in, infil_out, ...
    mann, S_min, h_min, nx, ny, dx, dy, dt, t, bc)
%=========================================================================
% This function solves the 2D overland flow model using implicit schemes
% for a given time step. 
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       topo            % [m] normalized topography
%       water_depth     % [m] depth of ponded water
%       ppt_in          % [m/hr] Input into soil layers from canopy
%       infil_out       % [m/hr] Drainage into subsurface
%       mann            % [hr/m^(1/3)] Mannings coefficient
%       S_min           % [-] minimum slope for water movement
%       h_min           % [m] miminum head difference for water movement
%       t               % [] current timestep
%       bc              % [] boundary condition
%------------------------- Output Variables -------------------------------
%       water_depth_new  % [m] new depth of ponded water
%=========================================================================

    % get water elevation (topo+depth)
    H = topo + water_depth;
    H_new = 0*H;
    D0 = 0*H;
    
    % get D overland
    % add boundaries to H
    H_bc = zeros(nx+2, ny+2);
    H_bc(2:nx+1, 2:ny+1) = H;
    
    if bc(1) == 1
        H_bc(1, 2:ny+1) = topo(1,:); 
    else 
        H_bc(1, 2:ny+1) = H(1,:);
    end
    
    if bc(2) == 1
        H_bc(nx+2, 2:ny+1) = topo(nx, :);
    else
        H_bc(nx+2, 2:ny+1) = H(nx, :);
    end
    
    H_bc(:,1) = H_bc(:,2);
    if bc(3) == 1
        H_bc(2:nx+1,1) = topo(:,1);
    end
    
    H_bc(:,ny+2) = H_bc(:,ny+1);
    if bc(4) == 1
        H_bc(2:nx+1,ny+2) = topo(:,ny);
    end
    
    
    [h_im1, h_ip1, h_jm1, h_jp1, h_im1jm1, h_im1jp1, h_ip1jm1, h_ip1jp1] = getKernel_wBC(H_bc, nx, ny);
   
    % get Sn
    Sn_n = sqrt(((h_im1jp1+h_jp1 - h_im1jm1-h_jm1)./(4*dx)).^2 + ((H - h_im1)./dy).^2);
    Sn_s = sqrt(((h_ip1jp1+h_jp1 - h_ip1jm1-h_jm1)./(4*dx)).^2 + ((h_ip1 - H)./dy).^2);
    Sn_w = sqrt(((H - h_jm1)./dx).^2 + ((h_ip1jm1+h_ip1 - h_im1jm1-h_im1)./(4*dy)).^2);
    Sn_e = sqrt(((h_jp1 - H)./dx).^2 + ((h_ip1jp1+h_ip1 - h_im1jp1-h_im1)./(4*dy)).^2);

    % get D
    [mann_im1, mann_ip1, mann_jm1, mann_jp1] = getSpatialKernel_2D(mann);
    [mann_n, mann_s, mann_w, mann_e] = getNSWE(mann, mann_im1, mann_ip1, mann_jm1, mann_jp1);

    water_depth_bc = zeros(nx+2, ny+2);
    water_depth_bc(2:nx+1, 2:ny+1) = water_depth;
    
    if bc(1) == 1
        water_depth_bc(1, 2:ny+1) = 0;
    else
        water_depth_bc(1, 2:ny+1) = water_depth(1,:);
    end
    water_depth_bc(nx+2, 2:ny+1) = water_depth(nx, :);
    water_depth_bc(:,1) = water_depth_bc(:,2);
    water_depth_bc(:,ny+2) = water_depth_bc(:,ny+1);
    
    [wd_im1, wd_ip1, wd_jm1, wd_jp1, ph1, ph2, ph3, ph4] = getKernel_wBC(water_depth_bc, nx, ny);  % only need first 4, ph = place holders
    [wd_n, wd_s, wd_w, wd_e] = getNSWE(water_depth, wd_im1, wd_ip1, wd_jm1, wd_jp1);

    D_n = getDiffusivity(h_im1, H, wd_im1, wd_n, water_depth, mann_n, Sn_n, S_min, h_min, D0);
    D_s = getDiffusivity(h_ip1, H, wd_ip1, wd_s, water_depth, mann_s, Sn_s, S_min, h_min, D0);
    D_w = getDiffusivity(h_jm1, H, wd_jm1, wd_w, water_depth, mann_w, Sn_w, S_min, h_min, D0);
    D_e = getDiffusivity(h_jp1, H, wd_jp1, wd_e, water_depth, mann_e, Sn_e, S_min, h_min, D0);
    
    % set up linear system
    % set matrix diagonals
    bx = D_s/(dx^2);
    cx = D_n/(dx^2);
    by = D_e/(dy^2);
    cy = D_w/(dy^2);
    a = 1/dt + bx + cx + by + cy;
    rhs_bc = 0*topo;    % init for boundary conditions for rhs

    % set boundaries
    % north, i = 1
    if bc(1) == 1
        rhs_bc(1,:) = rhs_bc(1,:) + cx(1,:) .* H_bc(1,2:ny+1);
        cx(1,:) = 0;
    else
        a(1,:) = a(1,:) - cx(1,:);
        cx(1,:) = 0; % set to 0 to make sure it is not included later
    end

    % south, i = nx
    if bc(2) == 1
        rhs_bc(nx,:) = rhs_bc(nx,:) + bx(nx,:) .* H_bc(nx+2, 2:ny+1);
        bx(nx,:) = 0;
    else
        a(nx,:) = a(nx,:) - bx(nx,:);
        bx(nx,:) = 0;
    end

    % west, j = 1
    if bc(3) == 1
        rhs_bc(:,1) = rhs_bc(:,1) + cy(:,1) .* H_bc(2:nx+1, 1);
        cy(:,1) = 0;
    else
        a(:,1) = a(:,1) - cy(:,1);
        cy(:,1) = 0;
    end
    
    % east, j = ny
    if bc(4) == 1
        rhs_bc(:,ny) = rhs_bc(:,ny) + by(:,ny) .* H_bc(2:nx+1, ny+2);
        by(:,ny) = 0;
    else
        a(:,ny) = a(:,ny) - by(:,ny);
        by(:,ny) = 0;
    end
    
    % form matrix
    d = [-nx,-1,0,1,nx];
    nxny = nx*ny;
    % nxny+1: make m<n to reverse truncation method of off diagonals
    A = spdiags([-cy(:), -cx(:), a(:), -bx(:), -by(:)],d,nxny,nxny+1); 
    A = A(1:nxny, 1:nxny);  % truncate A to remove extra column due to nxny+1
    
    % rhs with boundary flux 
    rhs = H(:)./dt + ppt_in(:) - infil_out(:) + rhs_bc(:);

    % Solve
    H_new(:) = A\rhs;

    % get water depth
    water_depth_new = H_new-topo;
    
%     % correct for negative water depth
%     total_neg = -1*sum(water_depth_new(water_depth_new < 0), 'all');
%     total_pos = sum(water_depth_new(water_depth_new > 0), 'all');
%     if (total_neg > 0 && total_neg <= total_pos)
%         pos_scale = (total_pos - total_neg)/total_pos;
%         water_depth_new(water_depth_new < 0) = 0;
%         water_depth_new = water_depth_new * pos_scale;
%         fprintf('Correct for negative water depth at tstep=%d\n', t)
%     elseif (total_neg > 0)
%         fprintf('Mass Imbalance in Overland Solution at tstep=%d\n', t)
%     end
    
%     q_out = D_n(1,:) .* water_depth_new(1,:) ./ dy; % [m^2/hr] match output from mass balance
end


function [im1, ip1, jm1, jp1, im1jm1, im1jp1, ip1jm1, ip1jp1] = getKernel_wBC(H_bc, nx, ny)
    im1 = H_bc(1:nx, 2:ny+1);
    ip1 = H_bc(3:nx+2, 2:ny+1);
    jm1 = H_bc(2:nx+1, 1:ny);
    jp1 = H_bc(2:nx+1, 3:ny+2);
    
    im1jm1 = H_bc(1:nx, 1:ny);
    im1jp1 = H_bc(1:nx, 3:ny+2);
    ip1jm1 = H_bc(3:nx+2, 1:ny);
    ip1jp1 = H_bc(3:nx+2, 3:ny+2);
end


function [im1,ip1,jm1,jp1] = getSpatialKernel_2D(data)
    [nx, ny] = size(data);
    % direct neighbors
    im1 = cat(1, data(1,:), data(1:nx-1,:));
    ip1 = cat(1, data(2:nx,:), data(nx,:));
    jm1 = cat(2, data(:,1), data(:,1:ny-1));
    jp1 = cat(2, data(:,2:ny), data(:,ny));
end

function [n,s,w,e] = getNSWE(data, im1,ip1,jm1,jp1)
    n = (data + im1)/2;
    s = (data + ip1)/2;
    w = (data + jm1)/2;
    e = (data + jp1)/2;
end

function D = getDiffusivity(h_next, h_cen, wd_next, wd_side, wd_cen, mann, Sn, S_min, h_min, D0)
    D = D0;
    % conditions for having D
    % 1. (h_next > h_cen and wd_next > h_min) or (h_next <= h_cen and wd_cen > h_min)
    % 2. wd_side > h_min
    % 3. abs(Sn) > S_min
    
    have_val = ( ((h_next > h_cen & wd_next > h_min) | (h_next <= h_cen & wd_cen > h_min) )...
        & wd_side > h_min & abs(Sn) > S_min);
    D(have_val) = (wd_side(have_val).^(5/3))./ (mann(have_val) .* sqrt(Sn(have_val)));
end