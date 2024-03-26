function [Q] = eval_interp_hcube_lagrange(xk0, x)
%EVAL_INTERP_HCUBE_LAGRANGE Evaluate interpolation functions from
%interpolation nodes in one-dimension (interpolation nodes of hypercube is
%the tensor product of these nodes).
%
%Input arguments
%---------------
%   XK0 : Array (NV0,) : Nodes at which Lagrange polynomials interpolate
%    values (nodes defining Lagrange polynomials) in each dimension.
%
%   X : Array (NDIM, NX) : Points at which Lagrange polynomials are
%     evaluated. 
%
%Output arguments
%----------------
%   Q : Array (NV, NDIM+1, NX) : Lagrange interpolation functions (NV) and
%     their derivative evaluated at each point in X. Q(i, 1, j) = value of
%     ith interpolation function evaluated at X(:, j), i.e., Q(i, 1, j) =
%     phi_i(X(j)). Q(i, 1+k, j) = kth partial derivative of ith
%     interpolation function evaluated at X(:, j), i.e., Q(i, 1+k, j) =
%     d(phi_i)/d(x_k))(X(:, j).

% Extract information from input
[ndim, nx] = size(x);
nv0 = numel(xk0);
nv = nv0^ndim;

if ndim == 0
    Q = ones(1, 1, nx);
    return;
end

% Construct 1d basis functions at each point (all dimensions)
[xr, ~, Ixr] = unique(x(:));
Q0_ = eval_interp_onedim_lagrange(xk0, xr);
Q0 = reshape(Q0_(:, 1, Ixr), [nv0, ndim, nx]);
dQ0 = reshape(Q0_(:, 2, Ixr), [nv0, ndim, nx]);

% Preallocate
Q = zeros(nv, ndim+1, nx);

% Tensor product of 1D basis functions
for k = 1:nx
    Q0_ = num2cell(Q0(:, :, k), 1);
    Q_ = tensprod_scalar_from_onedim(Q0_);
    Q(:, 1, k) = Q_(:);
    for j = 1:ndim
        dQ0_ = Q0_;
        dQ0_{j} = dQ0(:, j, k);
        dQ_ = tensprod_scalar_from_onedim(dQ0_);
        Q(:, 1+j, k) = dQ_(:);
    end
end

end