function [w_simp, z_simp] = create_quad_simp_from_hcube(w_hcube, z_hcube)
%CREATE_QUAD_SIMP_FROM_HCUBE Gaussian quadrature for a simplex from
%Gaussian quadrature for hypercube (collapse face).
%
%Input arguments
%---------------
%   W_HCUBE : Array (NQ,) : Quadrature weights of hypercube
%
%   Z_HCUBE : Array (NDIM, NQ) : Quadrature points of hypercube
%
%Output arguments
%----------------
%   W_SIMP : Array (NQ,) : Quadrature weights of simplex
%
%   Z_SIMP : Array (NDIM, NQ) : Quadrature points of simplex

% Transform quadrature over [-1, 1]^d hcube to [0, 1]^d hcube
ndim = size(z_hcube, 1);
w_hcube = w_hcube * 0.5^ndim;
z_hcube = 0.5*(z_hcube+1);

% Initialize simplex quadrature
w_simp = w_hcube;
z_simp = z_hcube;

% Map nodes to simplex and adjust weights
for k=1:size(z_hcube, 2)
    for j=1:ndim
        w_simp(k) = w_simp(k) * (1-z_hcube(j, k))^(ndim-j);
    end
    for j=ndim:-1:1
        for i=1:j-1
            z_simp(j, k) = z_simp(j, k)*(1-z_simp(i, k));
        end
    end
end