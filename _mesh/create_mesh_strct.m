function [msh] = create_mesh_strct(etype, xcg, e2vcg, e2bnd)
%CREATE_MESH Create mesh structure from the element type, nodes,
%connectivity, and boundary information.
%
% Input arguments
% ---------------
%   ETYPE, XCG, E2VCG, E2BND : See notation.m
%
% Output arguments
% ----------------
%   MSH : See notation.m

% Extract information from input
ndim = size(xcg, 1);
nnode_per_elem = size(e2vcg, 1);

% Determine polynomial order from number of nodes per element and etype
% and create geometry
if strcmpi(etype, 'hcube')
    porder = round((nnode_per_elem)^(1/ndim)-1);
elseif strcmpi(etype, 'simp')
    porder = 0;
    while true
        nv0 = nchoosek(porder+ndim, ndim);
        if nv0 == nnode_per_elem, break; end
        porder = porder+1;
    end
else
    error('Element not supported');
end
lfcnsp = create_polysp_nodal(etype, ndim, porder);

% Create geometry transformation data, setup face information
% (face-to-node, element-to-face, face-to-element), and create
% mesh structure
msh = struct('ndim', ndim, 'porder', porder, 'etype', etype, ...
             'xcg', xcg, 'e2vcg', e2vcg, 'e2bnd', e2bnd, 'lfcnsp', lfcnsp);

end