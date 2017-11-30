function vort = vort_mat(v,dX)
% Compute the 2D vorticity of a matrix field v
% We must have size(v) = [Mx My d n2]
% where Mx My are the space dimensions of R^2.
% The result will be of size [ Mx My 1 n2]
%

nabla_v = gradient_mat_2(v,dX); % [ Mx My d n2 d]
vort = nabla_v(:,:,1,:,2) - nabla_v(:,:,2,:,1);

end