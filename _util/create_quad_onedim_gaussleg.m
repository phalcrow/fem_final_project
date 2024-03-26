function [w, z] = create_quad_onedim_gaussleg(N)
%CREATE_QUAD_ONEDIM_GAUSSLEG Compute N-point Gauss-Legendre quadrature rule
%in one-dimension, i.e., Gaussian quadrature corresponding to Legendre
%orthogonal polynomials. Legrendre polynomials: (a) weighting function = 1,
%(b) interval = (-1, 1), (c) recursion: (2j-1)/j*x*p_{j-1}(x) -
%(j-1)/j*p_{j-2}(x).
%
%Input arguments
%---------------
%   N : Number of points in quadrature rule
%
%Output arguments
%----------------
%   W : Array (NQ,) : Quadrature weights
%
%   Z : Array (NQ,) : Quadrature points

mu0 = 2.0; % int_[-1, 1] 1 * dx = 2
j = (1:N)';
a = (2*j-1)./j; c = (j-1)./j;
b = zeros(N, 1);
[w, z] = quad_onedim_gauss_from_recursion(a, b, c, mu0);

end

function [w, z] = quad_onedim_gauss_from_recursion(a, b, c, mu0)
%QUAD_ONEDIM_GAUSS_FROM_RECURSION Compute N-point Gaussian quadrature rule
%from orthogonal polynomial recursion:
%    p_j(x) = (a_j*x+b_j)*p_{j-1}(x) - c_j*p_{j-2}(x)     (j = 1, ..., N)
%and integral of weighting function over domain: mu0 = int_a^b w(x) dx.
%
%Input arguments
%---------------
%   A, B, C, MU0 : Defined in recursion formula above.
%
%Output arguments
%----------------
%   W : Array (NQ,) : Quadrature weights
%   Z : Array (NQ,) : Quadrature points

d0 = -b./a; d1 = sqrt(c(2:end)./a(1:end-1)./a(2:end));
A = diag(d0, 0) + diag(d1, -1) + diag(d1, 1);
[v, lam] = eig(A);
z = diag(lam);
w = mu0*v(1, :).^2;

end