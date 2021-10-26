function [normA] = normalize_beta(A)  


normA = A - min(A(:));
normA = normA ./ max(normA(:));


end
