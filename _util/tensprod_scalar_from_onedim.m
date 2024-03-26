function [v] = tensprod_scalar_from_onedim(vk)
%TENSPROD_SCALAR_FROM_ONEDIM Compute a "scalar" tensor product from
%a sequence of vectors, i.e., v = vk{1} x ... x vk{end} or
%vk(i, j, ...) = vk{1}(i)*vk{2}(j)*...
%
%Input arguments
%---------------
%   VK : Cell array of 1D arrays : Vectors to use in tensor product
%
%Output arguments
%----------------
%   V : Array : Array resulting from tensor product

v = vk{1}; v = v(:);
for k = 2:numel(vk)
    vk_ = vk{k}; vk_ = vk_(:);
    v = v*vk_'; v = v(:);
end

end