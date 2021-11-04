% distribute vertical column over spatial area
function mat = ScatterCol2Mat(mat, vec, x_ind, y_ind, z_ind)
%=========================================================================
% This function takes a vertical vector and set every column of a 3D matrix
% to the vector value
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       mat             % [] empty (or non-empty) matrix to be filled
%       vec             % [s] given vector values to fill matrix 
%       x_ind           % [] x-index values of matrix to fill
%       y_ind           % [] y-index values of matrix to fill
%       z_ind           % [] z-index values of matrix to fill
%------------------------- Output Variables ------------------------------
%       mat             % [] filled matrix
%-------------------------------------------------------------------------           
    for k=z_ind
        mat(x_ind, y_ind, k) = vec(k);
    end
end