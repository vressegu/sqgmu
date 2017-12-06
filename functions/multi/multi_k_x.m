function A = multi_k_x(A)
% Generalisation of k x A (cross product)  if A is an ND array
%

idx =[];
for k=2:ndims(A)
   idx = [idx ',:']; 
end

eval(['A = A([2 1]' idx ');']);
eval(['A(1' idx ')  = - A(1' idx ');']);