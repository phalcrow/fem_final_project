function [U, info] = solve_linelast_hcyl0(c, r1, r2, h, nel, porder, pltit)
%SOLVE_LINELAST_HCYL0 Solve linear elasticity equation for deformation of
%hollow cylinder using FEM with ETYPE elements of polynomial completeness
%PORDER.
%
% Input arguments
% ---------------
%   C : Array (2,) : Center
%
%   R1, R2 : number : Inner/outer radius
%
%   H : number : Height
%
%   NEL : Array (NDIM,) : Number of elements in mesh in each direction
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

c = c(:);
ndim = 3;

% Create finite element mesh
msh = create_mesh3d_hollow_cylinder(c, r1, r2, h, nel, porder);
xcg = msh.xcg; e2vcg = msh.e2vcg; e2bnd = msh.e2bnd;
nnode = size(xcg, 2); f2v = msh.lfcnsp.f2v;

% Setup equation parameters and natural boundary conditions
prob.eqn = LinearElasticity(ndim);
prob.vol_pars_fcn = % TODO
prob.bnd_pars_fcn = % TODO

% Create finite element space
femsp = create_femsp_cg(prob, msh, porder, e2vcg);
ldof2gdof = femsp.ldof2gdof;

% Extract indices and set values of dirichlet boundary conditions
dbc_idx = % TODO
dbc_val = % TODO
dbc = create_dbc_strct(max(ldof2gdof(:)), dbc_idx, dbc_val);
femsp.dbc = dbc;

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, [], tol, maxit);

if pltit
    xcg_def = xcg+1*reshape(U, [ndim, nnode]);
    msh_def = create_mesh_strct('hcube', xcg_def, e2vcg, e2bnd);
    uabs = sqrt(U(1:3:end).^2+U(2:3:end).^2+U(3:3:end).^2);
    visualize_fem3d([], msh_def, uabs, struct('plot_elem', true), 1:4); colorbar;
    
    % Evaluate FEM solution throughout domain
    xeval = % TODO
    Ux = eval_fem_soln(U(ldof2gdof), xeval, msh, femsp.elem);

    figure;
    plot(xeval(3, :), squeeze(Ux(1, 1, :)), 'b-', 'linewidth', 2); hold on;
    plot(xeval(3, :), squeeze(Ux(2, 1, :)), 'k-', 'linewidth', 2);
    plot(xeval(3, :), squeeze(Ux(3, 1, :)), 'r-', 'linewidth', 2);
end

end