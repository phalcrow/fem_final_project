function [w, z] = create_quad_hcube_from_onedim(ndim, w0, z0)
%CREATE_QUAD_HCUBE_FROM_ONEDIM Compute NQ0^NDIM-point Gauss-Legendre
%quadrature rule in NDIM-dimensions, i.e., tensor product of
%one-dimensional Gaussian quadrature rules corresponding to Legendre
%orthogonal polynomials.
%
%Input arguments
%---------------
%   NDIM : Number of dimensions
%
%   W0 : Array (NQ0,) : 1D quadrature weights
%
%   Z0 : Array (NQ0,) : 1D quadrature nodes
%
%Output arguments
%----------------
%   W : Array (NQ,) : Quadrature weights
%
%   Z : Array (NDIM, NQ) : Quadrature points

% Treat 0-dimensional case (single point) as special case
if ndim == 0
    w = 1; z = zeros(0, 1);
    return;
end

% General case using tensor products
w = tensprod_scalar_from_onedim_unif(w0, ndim);
z = tensprod_vector_from_onedim_unif(z0, ndim);

end