function [U, info] = solve_fem(femsp, U0, tol, maxit)
%SOLVE_FEM_STATIC Solve nonlinear finite element problem given complete
%description of mesh, PDE, and boundary conditions (finite element space).
%
% Input arguments
% ---------------
%   FEMSP: See notation.m
%
%   U0 : Array (NDOF,) : Initial guess for vector of primary
%     variables at free indices (where essential boundary condition not
%     given).
%
%   TOL, MAXIT : See SOLVE_NEWTRAPH
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Approximate solution of PDE (primary variables) for
%     all global degrees of freedom (assembled).
%
%   INFO : See SOLVE_NEWTRAPH

% Extract information from input
ndof = max(femsp.ldof2gdof(:));

% Default arguments for u0, tol, maxit
if nargin < 2 || isempty(U0), U0 = zeros(ndof, 1); end
if nargin < 3, tol = 1.0e-8; end
if nargin < 4, maxit = 10; end

% Create finite element residual/Jacobian function
fcn = @(u_) create_fem_resjac(u_, femsp);

% % Check Jacobian with finite differences
% Ut = rand(size(Uf0));
% [R, dR] = fcn(Ut);
% dR_fd = compute_jac_findiff1(fcn, Ut, 1e-6);

% Solve nonlinear system using Newton-Raphson
Uu0 = U0(femsp.dbc.free_idx);
[Uu, info] = solve_newtraph(fcn, Uu0, tol, maxit);

% Assemble complete solution vector
U = zeros(ndof, 1);
U(femsp.dbc.dbc_idx) = femsp.dbc.dbc_val;
U(femsp.dbc.free_idx) = Uu;

end