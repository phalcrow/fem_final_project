function [Tv, Tvf] = create_elem_basis_mixed2(nvar1, Qv1, Qvf1, ...
                                              nvar2, Qv2, Qvf2)
%CREATE_ELEM_BASIS_MIXED2 Create element basis for mixed element (two
%different local function spaces) evaluated at point throughout volume (TV)
%and on each face (TVF) from basis of local function spaces evaluated at
%corresponding points (QV1, QVF1, QV2, QVF2). 
%
% Input arguments
% ---------------
%   NVAR1, QV1, QVF1 : NVAR, QV, QVF (see notation.m) for the first local
%     function space, e.g., velocity in Navier-Stokes.
%
%   NVAR2, QV2, QVF2 : NVAR, QV, QVF (see notation.m) for the second local
%     function space, e.g., pressure in Navier-Stokes.
%
% Output arguments
% ----------------
%   TV, TVF : See notation.m

% Extract information from input
ndim = size(Qv1, 2)-1;
[nv1, ~, nq] = size(Qv1);
[~, ~, nqf, nf] = size(Qvf1);
[nv2, ~] = size(Qv2);

% Element basis
% Code me!

end