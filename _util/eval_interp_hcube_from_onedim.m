function [Q] = eval_interp_hcube_from_onedim(ndim, Q1d)
%EVAL_INTERP_HCUBE_FROM_ONEDIM Evaluate interpolation functions for
%hypercube from interpolation functions in one-dimension. The nodes of
%the interpolation functions (and points at which they are evaluated) are
%inherited from the tensor product structure.
%
%Input arguments
%---------------
%   NDIM : Number of spatial dimensions
%
%   Q1D : Array (NV, 2, NX) : One-dimensional interpolation functions (NV)
%     and their derivative evaluated at NX points, X (array of size (NX,)).
%
%Output arguments
%----------------
%   Q : Array (NV^NDIM, NDIM+1, NX^NDIM) : Interpolation functions and
%     their derivatives for the hypercube element evaluated at each point
%     inherited from the 1D element.

% Extract information from input
dQ1d = Q1d(:, 2, :);
Q1d = Q1d(:, 1, :);
[nv, nx] = size(Q1d);

% Handle case where ndim == 0 as special case, assuming nv = nx = 1
if ndim == 0
    Q = ones(1, 1, 1);
    return;
end

% Prepare permute operation: basis fcns first, eval points last
sz = zeros(1, 2*ndim);
perm = zeros(1, 2*ndim);
for i=1:ndim
    sz(2*i-1) = nv;
    sz(2*i) = nx;
    perm(i) = 2*i-1;
    perm(ndim+i) = 2*i;
end

% Tensor product of 1D basis functions and fix ordering
Q_ = tensprod_scalar_from_onedim_unif(Q1d(:), ndim);
Q_ = reshape(Q_, sz);
Q_ = permute(Q_, perm);

Q = zeros(nv^ndim, ndim+1, nx^ndim);
Q(:, 1, :) = reshape(Q_, [nv^ndim, nx^ndim]);

% Tensor product of 1d basis function derivatives and fix ordering
v = cell(1, ndim);
for k = 1:ndim
    [v{:}] = deal(Q1d(:));
    v{k} = dQ1d(:);
    dQ_ = tensprod_scalar_from_onedim(v);
    dQ_ = reshape(dQ_, sz);
    dQ_ = permute(dQ_, perm);
    Q(:, 1+k, :) = reshape(dQ_, [nv^ndim, nx^ndim]);
end

end