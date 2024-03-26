function [zk, f2v, N] = create_nodes_bndy_refdom_simp(ndim, porder)
%CREATE_NODES_BNDY_REFDOM_SIMP Create nodal distribution and boundary of
%NDIM-dimensional simplex element of order PORDER.
%
%Input arguments
%---------------
%   NDIM, PORDER : See notation.m
%
%Output arguments
%----------------
%   ZK, F2V, N : See notation.m

% Treat 0-dimensional case (boundary of 1D element) as special case
if ndim == 0
    zk = zeros(0, 1); f2v = []; N = [];
    return;
end

% Extract information from input
nf = ndim+1;

% Code me!
nv = nchoosek(ndim + porder, ndim); % number of nodes per element
nvf = nchoosek(ndim + porder - 1, ndim - 1); % number of nodes per face

% Normal vectors
N = -eye(ndim);
N(:, nf) = ones(ndim, 1) / sqrt(ndim);

% Positions of the element nodes
if porder ~= 0
    v = linspace(0, 1, porder + 1); v = v(:); a = cell(1, ndim);
    [a{:}] = ndgrid(1:length(v)); zk = zeros(length(v)^ndim , ndim);
    for i = 1:ndim
        zk(:, i) = v(a{i}, :);
    end
    zk = zk';
    idx = sum(zk, 1) <= 1;
    zk = zk(:, idx); 
else
    zk = zeros(ndim, 1); f2v = [ones(1, ndim), NaN];
    return
end

% Face to vertices
f2v = zeros(nvf, nf);
for i = 1:ndim
    nds = find(zk(i, :) == 0);
    f2v(:, i) = nds';
end
f2v(:, nf) = find(sum(zk, 1) == 1);

% Checks
if size(zk, 2) ~= nv; error("Number of points of zk is incorrect."); end
if size(zk, 1) ~= ndim; error("Dimension of zk is incorrect."); end

end






















