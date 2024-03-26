function  [Me] = eval_unassembled_mass(transf_data, elem, elem_data)
%EVAL_UNASSEMBLED_MASS Evaluate/store element mass matrix for each element. 
%
% Input arguments
% ---------------
%   TRANSF_DATA, ELEM, ELEM_DATA : See notation.m 
%
% Output arguments
% ----------------
%   ME : Array (NDOF_PER_ELEM, NDOF_PER_ELEM, NELEM): Element mass matrix
%     for all elements in mesh

% Extract relevant variables from element
nelem = numel(elem_data);
nvar_per_elem = size(elem.Tv_ref, 1);
    
% Preallocate and assemble
Me = zeros(nvar_per_elem, nvar_per_elem, nelem);
for e = 1:nelem
    Me(:, :, e) = intg_elem_mass(transf_data(e), elem, elem_data(e));
end

end