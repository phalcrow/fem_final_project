function [varargout] = load_mesh(mshnm, etype, nref, porder)
%LOAD_MESH Load mesh from database.
%
% Input arguments
% ---------------
%   MSHNM : string : Prefix of mesh file
%
%   NREF : number : Level of refinement (usually = 0)
%
%   ETYPE, PORDER : See notation.m
%
%   DIRE: string : Path to _meshes directory. Default is .\_mesh\_meshes\
%     (if PC) and ./_mesh/_meshes/ (otherwise).
%
% Output arguments
% ----------------
%   MSH : See notation.m

delim = '/';
if ispc, delim='\'; end

dire=sprintf('_mesh%s_meshes%s', delim, delim);
fname = sprintf('%s%s-%s-nref%dp%d.mat', dire, mshnm, etype, nref, porder);
dat = load(fname);

% Return either msh or [xcg, e2vcg, e2bnd]
if nargout == 1
    msh = create_mesh_strct(etype, dat.xcg, dat.e2vcg, dat.e2bnd);
    varargout = {msh};
elseif nargout == 3
    varargout = {dat.xcg, dat.e2vcg, dat.e2bnd};
else
    error('Number of output arguments must be 1 (msh) or 3 (xcg, e2vcg, e2bnd)');
end

end