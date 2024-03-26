function [v] = tensprod_vector_from_onedim(vk)
%TENSPROD_VECTOR_FROM_ONEDIM Compute a "vector" tensor product of N
%vectors, i.e., v(:, i, j, ...) = (v{1}(i), v{2}(j), ...).
%
%Input arguments
%---------------
%   VK : Cell array of 1D arrays : Vectors to use in tensor product
%
%Output arguments
%----------------
%   V : Array resulting from tensor product

N = numel(vk);
out = cell(1, N);
[out{:}] = ndgrid(vk{:});
out = cellfun(@(M) M(:), out, 'UniformOutput', false);

v = [out{:}]';

end