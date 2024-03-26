function [lfcnsp] = create_polysp_nodal_simp(ndim, porder, z, r)
%CREATE_POLYSP_NODAL_SIMP Create structure defining master simplex element.
%
%Input arguments
%---------------
%   NDIM, PORDER : See notation.m
%
%   Z : Array (NDIM, NPT) : Points in master domain (\Omega_\square)
%     at which to evaluate basis functions
%
%   R : Array (NDIM, NPT_FC) : Point in master boundary (\Gamma_\square)
%     at which to evaluate boundary parametrization
%
%Output arguments
%----------------
%   LFCNSP : See notation.m

% Create geometry
[zk, f2v, N] = create_nodes_bndy_refdom_simp(ndim, porder);
[rk, ~, ~] = create_nodes_bndy_refdom_simp(ndim-1, porder);
if nargin<3 || isempty(z), z = zk; end
if nargin<4 || isempty(r), r = rk; end
nv = size(zk, 2); nf = size(f2v, 2); nr = size(r, 2);

% Basis functions, volume and faces
Qv = eval_interp_simp_lagrange(zk, z);

% Mapping from face parametrized domain (ndim-1) to actual domain (ndim)
[zk_p1, f2v_p1, ~] = create_nodes_bndy_refdom_simp(ndim, 1);
[rk_p1, ~, ~] = create_nodes_bndy_refdom_simp(ndim-1, 1);
Qf = eval_interp_simp_lagrange(rk_p1, r);
Qf = reshape(Qf, [size(Qf, 1), ndim*nr]);
r2z = zeros(ndim, ndim, nr, nf);
for f = 1:nf
    r2z(:, :, :, f) = reshape(zk_p1(:, f2v_p1(:, f))*Qf, [ndim, ndim, nr]);
end

% Basis functions, faces
Qvf = zeros(nv, ndim+1, nr, nf);
for f = 1:nf
    Qvf(:, :, :, f) = eval_interp_simp_lagrange(zk, squeeze(r2z(:, 1, :, f)));
end

% Create element structures
lfcnsp = struct('desc', 'P', 'porder', porder, 'N', N, ...
                'zk', zk, 'rk', rk, 'f2v', f2v, ...
                'Qv', Qv, 'Qvf', Qvf, 'r2z', r2z, ...
                'eval_basis_vol', @(z_) eval_interp_simp_lagrange(zk, z_), ...
                'eval_basis_face', @(r_) eval_interp_simp_lagrange(rk, r_));

end