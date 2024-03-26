function [dbc] = create_dbc_strct(ndof, dbc_idx, dbc_val)
%CREATE_DBC Create structure with Dirichlet boundary condition information.
%
% Input arguments
% ---------------
%   NDOF, DBC_IDX, DBC_VAL : See notation.m
%
% Output arguments
% ----------------
%   DBC : See notation.m

free_idx = setdiff((1:ndof)', dbc_idx);
dbc = struct('dbc_idx', dbc_idx, 'dbc_val', dbc_val, 'free_idx', free_idx);

end