function [U, info] = solve_linelast_beam0(xlims, ylims, etype, nel, porder, pltit)
%SOLVE_LINELAST_BEAM0 Solve linear elasticity equation for deformation of
%multimaterial beam using FEM with ETYPE elements of polynomial
%completeness PORDER.
%
% Input arguments
% ---------------
%   XLIMS, YLIMS : Array (2,) : Domain limits
%
%   NEL : Array (NDIM,) : Number of elements in mesh in each direction
%     (must be even in x-direction for multimaterial beam).
%
%   ETYPE, PORDER : See notation.m
%
%   PLTIT : bool : Whether to plot solution
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Global (assembled) finite element solution
%
%   INFO : See NEWTRAPH

ndim = 2;

% Check input
if rem(nel(1), 2)~=0, error('nelem(1) must be multiple of 2'); end

% Create finite element mesh
msh = create_mesh_hcube(etype, [xlims(:)'; ylims(:)'], nel, porder);
xcg = msh.xcg; e2vcg = msh.e2vcg; e2bnd = msh.e2bnd;
f2v = msh.lfcnsp.f2v;

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

% Visualize FEM solution over domain and along line
if pltit
    nnode = size(xcg, 2);
    xcg_def = xcg+1*reshape(U, [ndim, nnode]);
    msh_def = create_mesh_strct(etype, xcg_def, e2vcg, e2bnd);
    uabs = sqrt(U(1:2:end).^2+U(2:2:end).^2);
    visualize_fem([], msh_def, uabs(e2vcg), struct('nref', 2, 'plot_elem', true));
    colorbar;
    
    % Evaluate FEM solution throughout domain
    xeval = % TODO
    Ux = eval_fem_soln(U(ldof2gdof), xeval, msh, femsp.elem);

    figure;
    plot(xeval(1, :), squeeze(Ux(2, 1, :)), 'k-', 'linewidth', 2);
end

end