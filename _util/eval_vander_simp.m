function [V, dV] = eval_vander_simp(porder, x)
%VANDER_SIMP Compute NDIM-dimensional Vandermonde matrix of order PORDER
%for a  simplex and its derivative (NDIM determined from shape of X). The
%Vandermonde matrix, V, is the NX x M matrix  of multinomial terms and its
%derivative, dV, is the NX x M x NDIM matrix of the partial derivatives of
%multinomial terms, where M is the number of multinomial terms required for
%polynomial completeness (all combinations such that the sum of the
%exponents is <= porder), and X are the evaluation points.
%
%Input arguments
%---------------
%   PORDER : Polynomial degree of completeness
%
%   X : Array (NDIM, NX) : Points at which to evaluate multinomials in
%     definition of Vandermonde matrix 
%
%Output arguments
%----------------
%   V : Array (NX, M) : Vandermonde matrix
%
%   dV : Array (NX, M, NDIM) : Derivative of Vandermonde matrix

% Extract information from input
[ndim, nx] = size(x);

% Code me!

% Permissible exponents for basis functions
nv = nchoosek(ndim + porder, ndim);
p_ = (0:porder)';
a = cell(1, ndim); [a{:}] = ndgrid(1:porder + 1);
U0 = zeros(ndim, (porder + 1)^ndim);
for i = 1:ndim
    U0(i, :) = p_(a{i}, :)';
end
idx = find(sum(U0, 1) > porder);
Upsilon = U0(:, setdiff(1:size(U0, 2), idx));
if size(Upsilon, 2) ~= nv; error("Something is wrong."); end

% Vandermonde matrix (and its Jacobian)
V = ones(nx, nv);
dV = zeros(nx, nv, ndim);
for j = 1:nv
    for s = 1:ndim
        Usj = Upsilon(s, j);
        V(:, j) = V(:, j) .* x(s, :)'.^Usj;
        if Usj ~= 0
            dV(:, j, s) = Usj * x(s, :) .^ (Usj - 1);
            one2dim_but_s = setdiff(1:ndim, s);
            for m = one2dim_but_s
                dV(:, j, s) = dV(:, j, s) .* x(m, :).^Upsilon(m, j);
            end
        end
    end
end

end