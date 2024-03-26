function [zk, f2v, N] = create_nodes_bndy_refdom(etype, ndim, porder)
%CREATE_NODES_BNDY_REFDOM Create nodal distribution and boundary of
%NDIM-dimensional element of type ETYPE of order PORDER.
%
%Input arguments
%---------------
%   ETYPE, NDIM, PORDER : See notation.m
%
%Output arguments
%----------------
%   ZK, F2V, N : See notation.m

if strcmpi(etype, 'simp')
    [zk, f2v, N] = create_nodes_bndy_refdom_simp(ndim, porder);
elseif strcmpi(etype, 'hcube')
    [zk, f2v, N] = create_nodes_bndy_refdom_hcube(ndim, porder);
else
    error('Not supported');
end

end