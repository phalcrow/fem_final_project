function [U, e, info] = solve_pde0(nel, porder, pltit)
%SOLVE_PDE Solve PDE0.
%
% Input arguments
% ---------------
%   NEL : Array (2,) : Number of elements in each direction
%
%   PORDER : See notation.m
%
%   PLTIT : bool : Whether to plot solution
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Global (assembled) finite element solution
%
%   E : number : L2 error using mass matrix
%
%   INFO : See NEWTRAPH

nvar = 1;

% Create finite element mesh
msh = create_mesh_hcube('hcube', [0, 1], nel, porder);
xcg = msh.xcg; e2vcg = msh.e2vcg;

% Setup equation parameters and natural boundary conditions
prob.eqn = Pde0();
prob.vol_pars_fcn = @(x) x^2;
prob.bnd_pars_fcn = @(x, bnd) -1*(x>=1-1e-12);

% Extract indices and set values of dirichlet boundary conditions
dbc_idx = [1];
dbc_val = 0*dbc_idx;
dbc = create_dbc_strct(size(xcg, 2)*nvar, dbc_idx, dbc_val);

% Evaluate exact solution on mesh
Ue = (2*cos(1-xcg(1,:))-sin(xcg(1,:)))/(cos(1))+xcg(1,:).^2-2; Ue = Ue(:);

% Create finite element space
femsp = create_femsp_cg(prob, msh, porder, e2vcg, dbc);

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, [], tol, maxit);
                         
% Compute L2 error
M = eval_assembled_mass(femsp.transf_data, femsp.elem, femsp.elem_data, femsp.spmat);
E = U-Ue; e = sqrt(E'*M*E);

% Plot solution, if requested
if pltit
    ax = visualize_fem([], msh, U(e2vcg), struct('plot_elem', true));
    visualize_fem(ax, msh, Ue(e2vcg), struct('plot_elem', true));
end

end