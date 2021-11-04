function var_1d = GatherMat2Col(var_3d, BCond)
%=========================================================================
% This function averages each layer of a 3D matrix that is within the given
% boundary in the n, s, e, w directions
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       var_3d        % [] 3D matrix
%       BCond         % [] Boundary condition [top, n, s, e, w, bottom]
%------------------------- Output Variables ------------------------------
%       var_1d        % [] 1D vertical vector
%-------------------------------------------------------------------------  

    % take out edge cells based on BCond
    n_index = BCond(2)+1;
    s_index = size(var_3d, 1) - BCond(3);
    e_index = size(var_3d, 2) - BCond(4);
    w_index = BCond(5)+1;
    
    % clip 3d input and average
    var_3d_center = var_3d(n_index:s_index, w_index:e_index,:);
    var_1d = squeeze(mean(var_3d_center,[1,2]));
end