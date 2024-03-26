function [femsp] = create_femsp_cg(prob, msh, porder, e2vcg_porder, dbc)
%CREATE_FEMSP_CG Create finite element space (standard element; not mixed)
%
%Input arguments
%---------------
%   PROB, MSH, DBC : See notation.m
%
%   PORDER : Polynomial degree of solution (may be different than msh)
%
%   E2VCG_PORDER : E2VCG (see notation.m) for solution (may be different
%     than msh, but usually not, e.g., E2VCG_PORDER=MSH.E2VCG).
%
%  Note: If dbc requires information generated in this function, do not
%  pass it in (will be assigned [] value) and then manually set femsp.dbc
%  outside function.
%
%Output arguments
%----------------
%   FEMSP : See notation.m

% Default argument
if nargin < 4, e2vcg_porder = msh.e2vcg; end
if nargin < 5, dbc = []; end

% Extract information from input
nvar = prob.eqn.nvar;
etype = msh.etype;
ndim = size(msh.xcg, 1);

% Create quadrature rule
qrule = create_qrule_gaussleg(etype, ndim, 3*porder);

% Create local (scalar) function space for mesh and transformation data
lfcnsp_msh = create_polysp_nodal(etype, ndim, msh.porder, qrule.zq, qrule.rq);
transf_data = create_transf_data_ndim(lfcnsp_msh, msh.xcg, msh.e2vcg, msh.e2bnd);

% Master and physical elements for particular application
lfcnsp_soln = create_polysp_nodal(etype, ndim, porder, qrule.zq, qrule.rq);
elem = create_elem_strct(prob.eqn, qrule, lfcnsp_soln);
elem_data = create_elem_data(elem.Tv_ref, elem.Tvf_ref, ...
                             transf_data, prob.vol_pars_fcn, ...
                             prob.bnd_pars_fcn);

% Create  ldof2gdof and sparsity structure
ldof2gdof = create_ldof2gdof_cg(nvar, e2vcg_porder);
spmat = create_sparsity_strct(ldof2gdof);

% Store in class
femsp = struct('qrule', qrule, 'lfcnsp_msh', lfcnsp_msh, ...
               'lfcnsp_soln', lfcnsp_soln, 'transf_data', transf_data, ...
               'elem', elem, 'elem_data', elem_data, 'spmat', spmat, ...
               'ldof2gdof', ldof2gdof, 'dbc', dbc);

end