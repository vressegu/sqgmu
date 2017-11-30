function nabla_f = gradient_mat_2(f,dX)
% Compute the gradient of a matrix field f
% We must have size(f) = [Mx My n1 n2]
% where Mx My are the space dimensions of R^2.
% The result will be of size [ Mx My d n2 n1]
%

f=permute(f,[3 4 1 2]); % n1 n2 Mx My 
nabla_f = gradient_mat(f,dX); % n1 n2 Mx My 2
nabla_f=permute(nabla_f,[3 4 5 2 1]); % Mx My 2 n2 n1

end