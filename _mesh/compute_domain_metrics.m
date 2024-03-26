function [v, c, sa] = compute_domain_metrics(transf_data, qrule)
%COMPUTE_DOMAIN_METRICS Compute the volume, centroid, and surface
%area of a domain (approximated) by the mesh described by TRANSF_DATA.
%
%Input arguments
%---------------
%  TRANSF_DATA, QRULE : See notation.m
%
%Output arguments
%----------------
%  V : number : Volume of domain
%
%  C : Array (NDIM,) : Centroid of domain
%
%  SA : number : Surface area of domain

% Extract relevant quantities
ndim = size(transf_data(1).xe, 1);
nf = size(transf_data(1).sigf, 2);
nelem = numel(transf_data);

% Initialize volume, centroid, surface area
c = zeros(ndim, 1);
v = 0.0; sa = 0.0;

% Code me!

end