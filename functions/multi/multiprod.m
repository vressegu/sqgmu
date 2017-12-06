function A = multiprod(A,B)
% Generalisation of A * B if A and B are ND array
%

sA=size(A);
nA=length(sA);
sB=size(B);

if all( sA(2:nA) == sB([1 3:ndims(B)]) )
    A = permute(A, [2:nA 1]);
    B = permute(B, [1 3:(ndims(B)+1) 2]);
    A = bsxfun(@times,A,B);
    clear B;
    A=sum(A,1);
    A= permute(A, [nA nA+1 2:(nA-1) 1]);
    
else
    error('wrong size');
end