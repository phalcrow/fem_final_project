function [v] = tensprod_vector_from_onedim_unif(v1, N)
%TENSPROD_VECTOR_FROM_ONEDIM_UNIF Compute a "vector" tensor product of a
%vector with itself N times, i.e., v(:, i, j, ...) = (v1(i), v1(j), ...).
%
%Input arguments
%---------------
%   V1 : Array (m,) : Vector to use in tensor product
%
%   N : int : Number of times to tensor product v1 with itself
%
%Output arguments
%----------------
%   V : Array (m, m*...*m) : Array resulting from tensor product

in = cell(1, N);
[in{:}] = deal(v1);
v = tensprod_vector_from_onedim(in);

end