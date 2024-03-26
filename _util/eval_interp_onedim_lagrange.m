function [Q] = eval_interp_onedim_lagrange(xk, x)
%EVAL_INTERP_ONEDIM_LAGRANGE Evaluate Lagrange interpolation functions on
%interval (XK(1), XK(end)) associated with nodes XK (number of nodes = NV)
%at points X (number of points = NX).
%
%Input arguments
%---------------
%   XK : Array (NV,) : Nodes at which Lagrange polynomials interpolate
%    values (nodes defining Lagrange polynomials).
%
%   X : Array (NX,) : Points at which Lagrange polynomials are evaluated.
%
%Output arguments
%----------------
%   Q : Array (NV, 2, NX) : Lagrange interpolation functions (NV) and their
%     derivative evaluated at each point in X. Q(i, 1, j) = value of ith
%     interpolation function evaluated at X(j), i.e., Q(i, 1, j) =
%     phi_i(X(j)).  Q(i, 2, j) = derivative of ith interpolation function
%     evaluated at X(j), i.e., Q(i, 2, j) = (d(phi_i)/dx)(X(j)).

% Extract information from input, ensure xk and x appropriately shaped
nv = numel(xk);
nx = numel(x);
xk = xk(:); x = x(:);

% Evaluate Lagrange polynomials
tmp = ones(1, nx);
Q = ones(nv, nx);
dQ = zeros(nv, nx);
for j=1:nv
    for i=1:nv
        if i==j, continue; end
        Q(j, :) = Q(j, :) .* (x(:)'-xk(i))/(xk(j)-xk(i));
        
        tmp = 1 + 0*tmp;
        for m=1:nv
            if (m==j||m==i), continue; end
            tmp = tmp .* (x(:)'-xk(m))/(xk(j)-xk(m));
        end
        dQ(j, :) = dQ(j, :) + tmp/(xk(j)-xk(i));
    end
end
Q = reshape([Q; dQ], [nv, 2, nx]);
end