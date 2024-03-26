function [femsp] = create_femsp_cg_mixed2(prob, msh, ...
                                          nvar1, porder1, e2vcg1, ...
                                          nvar2, porder2, e2vcg2, dbc)
%CREATE_FEMSP_CG_MIXED2 Create finite element space (mixed (2) element)
%
%Input arguments
%---------------
%   PROB, MSH, DBC : See notation.m
%
%   NVAR1, PORDER1, E2VCG1 : Number of variables, polynomial degree, and
%     connectivity associated with first local function space
%
%   NVAR2, PORDER2, E2VCG2 : Number of variables, polynomial degree, and
%     connectivity associated with second local function space
%
%  Note: If dbc requires information generated in this function, do not
%  pass it in (will be assigned [] value) and then manually set femsp.dbc
%  outside function.
%
%Output arguments
%----------------
%   FEMSP : See notation.m

% Default argument
if nargin < 9, dbc = []; end

% Extract information from input
etype = msh.etype;
ndim = size(msh.xcg, 1);

% Create quadrature rule
qrule = create_qrule_gaussleg(etype, ndim, 3*porder1);

% Create local (scalar) function space for mesh and transformation data
lfcnsp_msh = create_polysp_nodal(etype, ndim, msh.porder, qrule.zq, qrule.rq);
transf_data = create_transf_data_ndim(lfcnsp_msh, msh.xcg, msh.e2vcg, msh.e2bnd);

% Master and physical elements for particular application
lfcnsp_soln1 = create_polysp_nodal(etype, ndim, porder1, qrule.zq, qrule.rq);
lfcnsp_soln2 = create_polysp_nodal(etype, ndim, porder2, qrule.zq, qrule.rq);
elem = create_elem_strct_mixed2(prob.eqn, qrule, nvar1, lfcnsp_soln1, ...
                                                 nvar2, lfcnsp_soln2);
elem_data = create_elem_data(elem.Tv_ref, elem.Tvf_ref, ...
                             transf_data, prob.vol_pars_fcn, ...
                             prob.bnd_pars_fcn);

% Create ldof2gdof and sparsity structure
ldof2gdof = create_ldof2gdof_cg_mixed2(nvar1, e2vcg1, nvar2, e2vcg2);
spmat = create_sparsity_strct(ldof2gdof);

% Store in class
femsp = struct('qrule', qrule, 'lfcnsp_msh', lfcnsp_msh, ...
               'lfcnsp_soln1', lfcnsp_soln1, 'lfcnsp_soln2', lfcnsp_soln2, ...
               'transf_data', transf_data, 'elem', elem, ...
               'elem_data', elem_data, 'spmat', spmat, ...
               'ldof2gdof', ldof2gdof, 'dbc', dbc);

end