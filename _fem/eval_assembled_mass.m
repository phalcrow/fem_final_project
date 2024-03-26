function [M] = eval_assembled_mass(transf_data, elem, elem_data, spmat)
%EVAL_ASSEMBLED_MASS Evaluate assembled mass matrix
%
%Input arguments
%---------------
%   TRANSF_DATA, ELEM, ELEM_DATA, SPMAT : See notation.m
%
%Output arguments
%----------------
%   M : Array (NDOF, NDOF) : Assembled mass matrix PRIOR to static condensation

Me = eval_unassembled_mass(transf_data, elem, elem_data);
M = assemble_nobc_mat(Me, spmat.cooidx, spmat.lmat2gmat);

end