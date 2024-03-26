function [elem] = create_elem_strct(eqn, qrule, lfcnsp)
%CREATE_ELEM_STRCT Create ELEM structure for system of PDEs where all
%components approximated with same basis (not mixed element).
%
% Input arguments
% ---------------
%   EQN, QRULE, LFCNSP : See notation.m
%
% Output arguments
% ----------------
%   ELEM : See notation.m

% Evaluate eqn-specific information (solution basis, face-to-element dof)
[Tv_ref, Tvf_ref] = create_elem_basis(eqn.nvar, lfcnsp.Qv, lfcnsp.Qvf);

% Create element structure
% [ndof_per_elem, nvar, ndimP1, ~] = size(Tv_ref); ndim = ndimP1-1;
elem = struct('eqn', eqn, 'qrule', qrule, ...
              'eval_elem_basis', @(z) eval_elem_basis(z, eqn.nvar, lfcnsp), ...
              'Tv_ref', Tv_ref, 'Tvf_ref', Tvf_ref);

end

function [Tv_ref] = eval_elem_basis(z, nvar, lfcnsp)
%EVAL_ELEM_BASIS Evaluate the solution basis given an arbitrary
%reference point (may not be quadrature node) for standard FE (not mixed).
%
% Input arguments
% ---------------
%   Z : Array (NDIM,) : Point at which to evaluate solution basis
%
%   NVAR, LFCNSP : See notation.m
%
% Output arguments
% ----------------
%   TV : Element basis and derivative evaluated at Z (see notation.m)

% Extract information from input
[nv, ndimP1, ~] = size(lfcnsp.Qv);
nf = size(lfcnsp.f2v, 2);
nz = size(z, 2);

% Evaluate geometry basis
Qv = lfcnsp.eval_basis_vol(z);
dum = zeros(nv, ndimP1, nz, nf);

% Evaluate solution basis
[Tv_ref, ~] = create_elem_basis(nvar, Qv, dum);

end