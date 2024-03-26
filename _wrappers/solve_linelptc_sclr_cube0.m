function [U, info] = solve_linelptc_sclr_cube0(nel, porder, pltit)
%SOLVE_LINELPTC_SCLR_CUBE0 Solve Poisson equation on unit cube, specifics
%of problem described in project handout.
%
% Input arguments
% ---------------
%   NEL : Array (3,) : Number of elements in each direction
%
%   PORDER : See notation.m
%
%   PLTIT : bool : Whether to plot solution
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Global (assembled) finite element solution
%
%   INFO : See NEWTRAPH

nvar = 1; ndim = 3;

% Create finite element mesh
msh = create_mesh_hcube('hcube',[0, 1; 0, 1; 0, 1], nel, porder);
xcg = msh.xcg; e2vcg = msh.e2vcg; e2bnd = msh.e2bnd;

% Setup equation parameters and natural boundary conditions
prob.eqn = LinearEllipticScalar(ndim);
prob.vol_pars_fcn = % TODO
prob.bnd_pars_fcn = % TODO

% Extract indices and set values of dirichlet boundary conditions
[~, f2v, ~] = create_nodes_bndy_refdom_hcube(ndim, porder);

dbc_idx = % TODO
dbc_val = % TODO
dbc = create_dbc_strct(size(xcg, 2)*nvar, dbc_idx, dbc_val);

% Create finite element space
femsp = create_femsp_cg(prob, msh, porder, e2vcg, dbc);

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, [], tol, maxit);

if pltit
    % Evaluate FEM solution throughout domain
    [Yeval, Xeval] = meshgrid(linspace(0, 1, 100));
    xeval = [Xeval(:), Yeval(:), 0.5+0*Xeval(:)]';
    Ux = eval_fem_soln(U(femsp.ldof2gdof), xeval, msh, femsp.elem);

    visualize_fem3d([], msh, U, [], 1:6);
    colorbar;
    
    figure;
    surf(Xeval, Yeval, reshape(Ux(1, 1, :), 100, 100));
end

end