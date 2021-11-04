function vec = MatrixVerticalAvg(mat)
%=========================================================================
% This function average each horizontal layer of a 3D matrix into a 1D
% vertical column vecter
%
% Written by Kunxuan Wang, UIUC, 2020
%
%------------------------- Input Variables -------------------------------
%       mat           % [] 3D matrix
%------------------------- Output Variables ------------------------------
%       vec           % [] vertical vector
%-------------------------------------------------------------------------  

    vec = squeeze(mean(mat,[1,2]));
end