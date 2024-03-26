function [qrule] = create_qrule_gaussleg_simp(ndim, nquad_per_dim)
%CREATE_QRULE_GAUSSLEG_SIMP Compute Gauss-Legendre quadrature rule for
%simplex.
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

% Tensor product quadrature rule for hypercube
[w_hcube, z_hcube] = create_quad_hcube_from_onedim(ndim, w0, z0);
[wf_hcube, r_hcube] = create_quad_hcube_from_onedim(ndim-1, w0, z0);

% Map hypercube to simplex for volume
[wq, zq] = create_quad_simp_from_hcube(w_hcube, z_hcube);
[wqf, rq] = create_quad_simp_from_hcube(wf_hcube, r_hcube);

% Create quadrature rule structure
qrule = struct('wq', wq, 'zq', zq, 'wqf', wqf, 'rq', rq);

end