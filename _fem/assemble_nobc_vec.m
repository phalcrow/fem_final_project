function [F] = assemble_nobc_vec(Fe, ldof2gdof)
%ASSEMBLE_NOBC_VEC Assemble element vector into a global vector without
%applying Dirichlet/essential boundary conditions.
%
%Input arguments
%---------------
%  FE : Array (NDOF_PER_ELEM, NELEM): Element vector for all
%    elements in mesh.
%
%  LDOF2GDOF : See notation.m
%
%Output arguments
%----------------
%   F : Array (NDOF,) : Assembled vector PRIOR to static condensation

% Extract quantities
ndof = max(ldof2gdof(:));
nelem = size(Fe, 2);

% Preallocate K, F
F = zeros(ndof, 1);

% Assemble over each element
for e = 1:nelem    
    % Update entries in global force vector
    idx = ldof2gdof(:, e);
    F(idx) = F(idx) + Fe(:, e);
end

end