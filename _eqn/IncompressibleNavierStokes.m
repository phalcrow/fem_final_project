function [claw] = IncompressibleNavierStokes(ndim)
%IncompressibleNavierStokes A class containing relevant information for the
%incompressible Navier-Stokes equation.
%
%Input arguments
%---------------
%   NDIM : See notation.m
%
%Output arguments
%----------------
%   CLAW : See notation.m 

claw = struct('nvar', ndim+1, 'ndim', ndim, ...
              'npars', 2, 'srcflux', @eval_ins_srcflux);

end

function [S, dSdU, dSdQ, F, dFdU, dFdQ] = eval_ins_srcflux(U, Q, pars)
%EVAL_INS_SRCFLUX Evaluate incompressible Navier-Stokes source term, flux
%function and their partial derivatives at a point.
%
%Input arguments
%---------------
%   U : Array (NVAR, 1) : PDE state vector at a point
%
%   Q : Array (NVAR, NDIM) : Gradient of PDE state vector at a point
%
%   PARS : Array (NPARS, 1) : Vector of parameters at a point
%
%Output arguments
%----------------
%   S, dSdU, dSdQ, F, dFdU, dFdQ : See notation.m

% Extract information from input
u = U; q = Q;
ndim = size(u, 1)-1;

% Define information regarding size of the system
neqn = ndim+1; ncomp = ndim+1;

% Extract parameters, primary variables, and gradient
rho = pars(1);
nu = pars(2);
v = u(1:ndim); p = u(end);
dv = q(1:ndim, :);

% Preallocate
S = zeros(neqn, 1);
dSdU = zeros(neqn, nvar);
dSdQ = zeros(neqn, nvar, ndim);

F = zeros(neqn, ndim);
dFdU = zeros(neqn, ndim, nvar);
dFdQ = zeros(neqn, ndim, nvar, ndim);

% Code me!

end