function [claw] = LinearElasticity(ndim)
%LinearElasticity A class containing relevant information for the
%linear elasticity PDE.
%
%Input arguments
%---------------
%   NDIM : See notation.m
%
%Output arguments
%----------------
%   CLAW : See notation.m 

claw = struct('nvar', ndim, 'ndim', ndim, ...
              'npars', ndim+2, 'srcflux', @eval_linelast_srcflux);

end

function [S, dSdU, dSdQ, F, dFdU, dFdQ] = eval_linelast_srcflux(U, Q, pars)
%EVAL_LINELAST_SRCFLUX Evaluate linear elasticity source term, flux
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
ndim = size(Q, 1);

% Define information regarding size of the system
neqn = ndim; nvar = ndim;

% Extract parameters
lam = pars(1);
mu = pars(2);
f = pars(3:3+ndim-1);

% Preallocate
S = zeros(neqn, 1);
dSdU = zeros(neqn, nvar);
dSdQ = zeros(neqn, nvar, ndim);

F = zeros(neqn, ndim);
dFdU = zeros(neqn, ndim, nvar);
dFdQ = zeros(neqn, ndim, nvar, ndim);

% Code me!
end