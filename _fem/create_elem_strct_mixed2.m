function [elem] = create_elem_strct_mixed2(eqn, qrule, nvar1, lfcnsp1, ...
                                                       nvar2, lfcnsp2)
%CREATE_ELEM_STRCT_MIXED2 Create ELEM structure for mixed element (two
%different local function spaces).
%
% Input arguments
% ---------------
%   EQN, QRULE : See notation.m
%
%   NVAR1, LFCNSP1 : NVAR, LFCNSP (see notation.m) for first local function
%     space used in mixed element
%
%   NVAR2, LFCNSP2 : NVAR, LFCNSP (see notation.m) for second local function
%     space used in mixed element
%
% Output arguments
% ----------------
%   ELEM : See notation.m

% Evaluate eqn-specific information (solution basis, face-to-element dof)
[Tv_ref, Tvf_ref] = create_elem_basis_mixed2(nvar1, lfcnsp1.Qv, lfcnsp1.Qvf, ...
                                             nvar2, lfcnsp2.Qv, lfcnsp2.Qvf);
                                                 
% Create element structure
elem = struct('eqn', eqn, 'qrule', qrule, ...
              'eval_elem_basis', @(z) eval_elem_basis(z, nvar1, lfcnsp1, nvar2, lfcnsp2), ...
              'Tv_ref', Tv_ref, 'Tvf_ref', Tvf_ref);
end

function [TV_ref] = eval_elem_basis(z, nvar1, lfcnsp1, nvar2, lfcnsp2)
%EVAL_ELEM_BASIS Evaluate the solution basis given an arbitrary
%reference point (may not be quadrature node) for a mixed (2) element.
%
% Input arguments
% ---------------
%   Z : Array (NDIM,) : Point at which to evaluate solution basis
%
%   NVAR1, LFCNSP1, NVAR2, LFCNSP2 : See main function
%
% Output arguments
% ----------------
%   TV : Element basis and derivative evaluated at Z (see notation.m)

% Extract information from input
[nv1, ndimP1, ~, nf] = size(lfcnsp1.Qvf);
nv2 = size(lfcnsp2.Qv, 1);
nz = size(z, 2);

% Evaluate geometry basis
Qv1 = lfcnsp1.eval_basis_vol(z);
Qv2 = lfcnsp2.eval_basis_vol(z);
dum1 = zeros(nv1, ndimP1, nz, nf);
dum2 = zeros(nv2, ndimP1, nz, nf);

% Evaluate solution basis
[TV_ref, ~] = create_elem_basis_mixed2(nvar1, Qv1, dum1, ...
                                       nvar2, Qv2, dum2);

end