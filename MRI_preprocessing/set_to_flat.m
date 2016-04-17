function [out] = set_to_flat(tc,thresh)

if ~exist('thresh','var')
    thresh = 0.1;
end
out = tc;
% Find flat tcs, set to zero
flattc = var(tc)<thresh;
numflat = find(flattc);
out(:,flattc) = zeros(size(tc,1),size(numflat,2));
