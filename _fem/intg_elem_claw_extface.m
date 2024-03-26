function [Re] = intg_elem_claw_extface(transf_data, elem, elem_data)
%INTG_ELEM_CLAW_EXTFACE Integrate element Galerkin form (boundary term) to
%form the boundary contribution to the element residual and Jacobian.
%
% Input arguments
% ---------------
%   TRANSF_DATA, ELEM, ELEM_DATA : See notation.m
%
% Output arguments
% ----------------
%   RE : Array (NDOF_PER_ELEM,) : Element residual (boundary contribution)

% Extract information from input : sizes
[nvar_per_elem, ~, ~, nqf, nf] = size(elem.Tvf_ref);
e2bnd = transf_data.e2bnd;

% Extract information from input : quadrature
wqf = elem.qrule.wqf;

% Extract information from input : isoparametric
sigf = transf_data.sigf;

% Preallocate element residual and Jacobian
Re = zeros(nvar_per_elem, 1);

% Integrate boundary term
% Code me!

end