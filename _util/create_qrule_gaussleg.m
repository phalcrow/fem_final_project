function [qrule] = create_qrule_gaussleg(etype, ndim, nquad_per_dim)
%CREATE_QRULE_GAUSSLEG Compute Gauss-Legendre quadrature rule.
%
%Input arguments
%---------------
%   ETYPE : See notation.m
%
%   NDIM : Number of dimensions of simplex
%
%   NQUAD_PER_DIM : Number of quadrature nodes per spatial dimension
%
%Output arguments
%----------------
%   QRULE : See notation.m

if strcmpi(etype, 'simp')
    qrule = create_qrule_gaussleg_simp(ndim, nquad_per_dim);
elseif strcmpi(etype, 'hcube')
    qrule = create_qrule_gaussleg_hcube(ndim, nquad_per_dim);
else
    error('Not supported');
end

end