function b = multi_linsolve(A,b,opts)
% Generalisation of linsolve(A,b,opts) if A is a square matrix and b is an ND array
%

sA = size(A);
sb = size(b);
if ~ ( ismatrix(A) && sA(1)==sA(2) && sA(2) == sb(1))
    error('wrong size');
end

b = reshape(b, [sb(1) prod(sb(2:end))]);
if nargin == 2
    b = linsolve(A,b);
else
    b = linsolve(A,b,opts);
end
clear A opts;
b = reshape(b, sb);
