function [lfcnsp] = create_polysp_nodal_hcube(ndim, porder, z, r)
%CREATE_POLYSP_NODAL_HCUBE Create structure defining master
%hypercube element.
%
%Input arguments
%---------------
%   NDIM, PORDER : See notation.m
%
%   Z : Array (NDIM, NPT) : Points in master domain (\Omega_\square)
%     at which to evaluate basis functions
%
%   R : Array (NDIM, NPT_FC) : Point on master boundary (\Gamma_\square)
%     at which to evaluate boundary parametrization
%
%Output arguments
%----------------
%   LFCNSP : See notation.m

% Create geometry
[zk, f2v, N] = create_nodes_bndy_refdom_hcube(ndim, porder);
[rk, ~, ~] = create_nodes_bndy_refdom_hcube(ndim-1, porder);
if nargin<3 || isempty(z), z = zk; end
if nargin<4 || isempty(r), r = rk; end
zk0 = unique(zk(1, :));
nv = size(zk, 2); nf = size(f2v, 2); nr = size(r, 2);

% Basis functions, volume
Qv = eval_interp_hcube_lagrange(zk0, z);

% Mapping from face parametrized domain (ndim-1) to actual domain (ndim)
[zk_p1, f2v_p1, ~] = create_nodes_bndy_refdom_hcube(ndim, 1);
zk0_p1 = unique(zk_p1(1, :));
Qf = reshape(eval_interp_hcube_lagrange(zk0_p1, r), [2^(ndim-1), ndim*nr]);
r2z = zeros(ndim, ndim, nr, nf);
for f = 1:nf
    r2z(:, :, :, f) = reshape(zk_p1(:, f2v_p1(:, f))*Qf, [ndim, ndim, nr]);
end

% Basis functions, faces
Qvf = zeros(nv, ndim+1, nr, nf);
for f = 1:nf
    Qvf(:, :, :, f) = eval_interp_hcube_lagrange(zk0, squeeze(r2z(:, 1, :, f)));
end

% Create local function space structure
lfcnsp = struct('desc', 'Q', 'porder', porder, ...
                'zk', zk, 'rk', rk, 'f2v', f2v, 'N', N, ...
                'Qv', Qv, 'Qvf', Qvf, 'r2z', r2z, ...
                'eval_basis_vol', @(z_) eval_interp_hcube_lagrange(zk0, z_), ...
                'eval_basis_face', @(r_) eval_interp_hcube_lagrange(zk0, r_));

end