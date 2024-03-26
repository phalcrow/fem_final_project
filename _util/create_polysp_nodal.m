function [lfcnsp] = create_polysp_nodal(etype, ndim, porder, z, r)
%CREATE_POLYSP_NODAL Create structure defining master element.
%
%Input arguments
%---------------
%   ETYPE, NDIM, PORDER : See notation.m
%
%   Z : Array (NDIM, NPT) : Points in master domain (\Omega_\square)
%     at which to evaluate basis functions (defaults to ZK).
%
%   R : Array (NDIM, NPT_FC) : Point on master boundary (\Gamma_\square)
%     at which to evaluate boundary parametrization (defaults to RK).
%
%Output arguments
%----------------
%   LFCNSP : See notation.m

if nargin < 4, z = []; end
if nargin < 5, r = []; end

if strcmpi(etype, 'simp')
    lfcnsp = create_polysp_nodal_simp(ndim, porder, z, r);
elseif strcmpi(etype, 'hcube')
    lfcnsp = create_polysp_nodal_hcube(ndim, porder, z, r);
else
    error('Not supported');
end

end