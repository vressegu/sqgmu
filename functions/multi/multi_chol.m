function L = multi_chol(A)
% Generalisation of chol(A) if A is an ND array
% Cholesky-Crout algorithm (column by column)
% L' * L = A
%

idx =[];
for k=3:ndims(A)
   idx = [idx ',:']; 
end

sA = size(A);
if sA(1) ~= sA(2)
    error('wrong size');
end
L = zeros(sA);

d=sA(1);
eval(['L(1,1' idx ') = real(sqrt( A(1,1' idx ')));']);
eval(['L(2:end,1' idx ') =  A(2:end,1' idx ');']);
eval(['den = 1./L(1,1' idx ');']);
eval(['L(2:end,1' idx ') =  bsxfun(@times, L(2:end,1' idx '), den);']);
for j=2:(d-1)
    eval(['L(j,j' idx ') = real(sqrt( A(j,j' idx ') - multiprod( ' ...
        'L(j,1:(j-1)' idx ') , multitrans ( L(j,1:(j-1)' idx ') ) ) ));']);
    eval(['L((j+1):end,j' idx ') = A((j+1):end,j' idx ') - multiprod( ' ...
        'L((j+1):end,1:(j-1)' idx ') , multitrans ( L(j,1:(j-1)' idx ') ) );']);
    eval(['den = 1./L(j,j' idx ');']);
    eval(['L((j+1):end,j' idx ') =  bsxfun(@times, L((j+1):end,j' idx '), den);']);
end
eval(['L(d,d' idx ') = real(sqrt( A(d,d' idx ') - multiprod( ' ...
    'L(d,1:(d-1)' idx ') , multitrans ( L(d,1:(d-1)' idx ') ) )) );']);
clear A;
L = multitrans(L);

