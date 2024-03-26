function [qrule] = create_qrule_gaussleg_hcube(ndim, nquad_per_dim)
%CREATE_QRULE_GAUSSLEG_HCUBE Compute Gauss-Legendre quadrature rule for
%hypercube.
%
%Input arguments
%---------------
%   NDIM : Number of dimensions of simplex
%
%   NQUAD_PER_DIM : Number of quadrature nodes per spatial dimension
%
%Output arguments
%----------------
%   QRULE : See notation.m

% One-dimensional quadrature
[w0, z0] = create_quad_onedim_gaussleg(nquad_per_dim);

% Tensor product quadrature rule for volume
[wq, zq] = create_quad_hcube_from_onedim(ndim, w0, z0);

% Tensor product quadrature rule in (ndim-1)-dimensions for faces
[wqf, rq] = create_quad_hcube_from_onedim(ndim-1, w0, z0);

% Create quadrature rule structure
qrule = struct('wq', wq, 'zq', zq, 'wqf', wqf, 'rq', rq);

end