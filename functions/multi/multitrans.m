function A = multitrans(A)
% Generalisation of A.' if A is an ND array
%

A = permute(A, [ 2 1 3:ndims(A)]);