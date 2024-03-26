function [M] = assemble_nobc_mat(Me, cooidx, lmat2gmat)
%ASSEMBLE_NOBC_MAT Assemble element matrices into a global matrix without
%applying Dirichlet/essential boundary conditions; use sparse format for
%matrix. 
%
%Input arguments
%---------------
%  ME : Array (NDOF_PER_ELEM, NDOF_PER_ELEM, NELEM): Element matrix for
%    all elements in mesh
%
%  COOIDX, LMAT2GMAT : See notation.m
%
%Output arguments
%----------------
%   M : Array (NDOF, NDOF) : Assembled matrix PRIOR to static condensation

% Extract quantities
nnz = size(cooidx, 1);

% Preallocate M
Mval = zeros(nnz, 1);

% Assemble over each element
nelem = size(Me, 3);
for e = 1:nelem    
    % Update entries in global stiffness matrix
    idx = lmat2gmat(:, :, e);
    Mval(idx) = Mval(idx) + Me(:, :, e);
end
M = sparse(cooidx(:, 1), cooidx(:, 2), Mval);

end