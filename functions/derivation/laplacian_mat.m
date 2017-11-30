function d2f = laplacian_mat(f,dX)
% Compute the laplacian of a matrix field f
% We must have size(f) = [d1 d d_{2+1} ... d_{2+d}]
% where d_{2+1} ... d_{2+d} are the space dimensions of R^d.
% The result will be of the same size
%

% warning(['Limit conditions are not taken into account.' ...
%     'So, the laplacian linearly extrapolate the interior values to compute the boundary values']);

siz = size(f);
d=length(siz)-2;
d2f=zeros([siz d]);
if d == 2
    % Compute the gradient along x
    dx=dX(1);
    d2f(:,:,2:end-1,:,1) = 1/dx^2 * (f(:,:,3:end,:)+f(:,:,1:end-2,:)-2*f(:,:,2:end-1,:));
    d2f(:,:,1,:,1) = 2 * d2f(:,:,2,:,1) - d2f(:,:,3,:,1);
    d2f(:,:,end,:,1) = -d2f(:,:,end-2,:,1) + 2 * d2f(:,:,end-1,:,1);
    
    % Compute the gradient along y
    dy=dX(2);
    d2f(:,:,:,2:end-1,2) = 1/dy^2 * (f(:,:,:,3:end)+f(:,:,:,1:end-2)-2*f(:,:,:,2:end-1));
    d2f(:,:,:,1,2) = 2 * d2f(:,:,:,2,2) - d2f(:,:,:,3,2);
    d2f(:,:,:,end,2) = -d2f(:,:,:,end-2,2) + 2 * d2f(:,:,:,end-1,2);
    
else
    % Compute the gradient along x
    dx=dX(1);
    d2f(:,:,2:end-1,:,:,1) = 1/dx^2 * (f(:,:,3:end,:,:)+f(:,:,1:end-2,:,:)-2*f(:,:,2:end-1,:,:));
    d2f(:,:,1,:,:,1) = 2 * d2f(:,:,2,:,:,1) - d2f(:,:,3,:,:,1);
    d2f(:,:,end,:,:,1) = -d2f(:,:,end-2,:,:,1) + 2 * d2f(:,:,end-1,:,:,1);
    
    % Compute the gradient along y
    dy=dX(2);
    d2f(:,:,:,2:end-1,:,2) = 1/dy^2 * (f(:,:,:,3:end,:)+f(:,:,:,1:end-2,:)-2*f(:,:,:,2:end-1,:));
    d2f(:,:,:,1,:,2) = 2 * d2f(:,:,:,2,:,2) - d2f(:,:,:,3,:,2);
    d2f(:,:,:,end,:,2) = -d2f(:,:,:,end-2,:,2) + 2 * d2f(:,:,:,end-1,:,2);
    
    % Compute the gradient along z
    dz=dX(3);
    d2f(:,:,:,:,2:end-1,3) = 1/dz^2 * (f(:,:,:,:,3:end)+f(:,:,:,:,1:end-2)-2*f(:,:,:,:,2:end-1));
    d2f(:,:,:,:,1,3) = 2 * d2f(:,:,:,:,2,3) - d2f(:,:,:,:,3,3);
    d2f(:,:,:,:,end,3) = -d2f(:,:,:,:,end-2,3) + 2 * d2f(:,:,:,:,end-1,3);
end

d2f = sum(d2f,3+d);
end