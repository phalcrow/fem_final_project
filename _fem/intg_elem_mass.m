function [Me] = intg_elem_mass(transf_data, elem, elem_data)
%INTG_ELEM_MASS Integrate (single) element mass matrix
%
% Input arguments
% ---------------
%   TRANSF_DATA, ELEM, ELEM_DATA : See notation.m
%     TRANSF_DATA/ELEM_DATA in this function  is a *single* entry of the
%     TRANSF_DATA/ELEM_DATA structure array, i.e., TRANSF_DATA(e)/ELEM_DATA(e)
%     where e is the element number.
%
% Output arguments
% ----------------
%   ME : Array (NDOF_PER_ELEM, NDOF_PER_ELEM) : Element mass matrix

% Extract information from input
sz = size(elem_data.Tv_phys);
Tvar = reshape(elem_data.Tv_phys(:, :, 1, :), [sz(1), sz(2), sz(4)]);
[nvar_per_elem, nvar, nq] = size(Tvar);

% Default arguments
if nargin < 5, coeff = ones(nq, 1); end

% Compute modified weights, including volume scale and ceofficients
wq = elem.qrule.wq;
detG = transf_data.detG;
w_sqrt = sqrt(wq.*detG.*coeff);

% Modify bases to include weights
Tvar = bsxfun(@times, Tvar, reshape(w_sqrt, [1, 1, nq]));
Tvar = reshape(Tvar, [nvar_per_elem, nvar*nq]);

% Integrate mass matrix with single matmat
Me = Tvar*Tvar';

end