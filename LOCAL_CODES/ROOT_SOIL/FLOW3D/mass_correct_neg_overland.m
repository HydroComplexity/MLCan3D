function water_depth = mass_correct_neg_overland(water_depth, err_msg)
%=========================================================================
% This function corrects for numerically induced negative overland ponded
% water values by filling negative with positive values. Maintains mass
% balance. 
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       water_depth   % [m] ponded depth of water in overland flow
%       err_msg       % [] String containing an error message
%------------------------- Output Variables ------------------------------
%       water_depth   % [m] corrected ponded depth of water in overland flow
%-------------------------------------------------------------------------  
    tol = 1e-5; % tolerace for 0
    
    % correct for negative water depth
    total_neg = -1*sum(water_depth(water_depth < 0), 'all');
    total_pos = sum(water_depth(water_depth > 0), 'all');
    if (total_neg > 0 && total_neg <= total_pos + tol)
        pos_scale = (total_pos - total_neg)/total_pos;
        water_depth(water_depth < 0) = 0;
        water_depth = water_depth * pos_scale;
    elseif (total_neg > tol)  % have some negatice values but not less than positive
        %fprintf('Mass Imbalance in Overland Solution at tstep=%d\n', t)
        disp(err_msg);
    end
end