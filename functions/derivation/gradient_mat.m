function nabla_f = gradient_mat(f,dX)
% Compute the gradient of a matrix field f
% We must have size(f) = [m d Mx My (Mz)]
% where Mx My (Mz) are the space dimensions of R^d.
% The result will be of size [ size(f) d]
%

%%

siz = size(f);
d=length(siz)-2;
nabla_f=zeros([siz d]);

if d == 2
    % Compute the gradient along x
    dx=dX(1);
    nabla_f(:,:,2:end-1,:,1) = 1/(2*dx) * (f(:,:,3:end,:)-f(:,:,1:end-2,:));
    nabla_f(:,:,1,:,1) = 1/(dx) * (f(:,:,2,:)-f(:,:,1,:));
    nabla_f(:,:,end,:,1) = 1/(dx) * (f(:,:,end,:)-f(:,:,end-1,:));
    
    % Compute the gradient along y
    dy=dX(2);
    nabla_f(:,:,:,2:end-1,2) = 1/(2*dy) * (f(:,:,:,3:end)-f(:,:,:,1:end-2));
    nabla_f(:,:,:,1,2) = 1/(dy) * (f(:,:,:,2)-f(:,:,:,1));
    nabla_f(:,:,:,end,2) = 1/(dy) * (f(:,:,:,end)-f(:,:,:,end-1));
    
else
    % Compute the gradient along x
    dx=dX(1);
    nabla_f(:,:,2:end-1,:,:,1) = 1/(2*dx) * (f(:,:,3:end,:,:)-f(:,:,1:end-2,:,:));
    nabla_f(:,:,1,:,:,1) = 1/(dx) * (f(:,:,2,:,:)-f(:,:,1,:,:));
    nabla_f(:,:,end,:,:,1) = 1/(dx) * (f(:,:,end,:,:)-f(:,:,end-1,:,:));
    
    % Compute the gradient along y
    dy=dX(2);
    nabla_f(:,:,:,2:end-1,:,2) = 1/(2*dy) * (f(:,:,:,3:end,:)-f(:,:,:,1:end-2,:));
    nabla_f(:,:,:,1,:,2) = 1/(dy) * (f(:,:,:,2,:)-f(:,:,:,1,:));
    nabla_f(:,:,:,end,:,2) = 1/(dy) * (f(:,:,:,end,:)-f(:,:,:,end-1,:));
    
    % Compute the gradient along z
    dz=dX(3);
    nabla_f(:,:,:,:,2:end-1,3) = 1/(2*dz) * (f(:,:,:,:,3:end)-f(:,:,:,:,1:end-2));
    nabla_f(:,:,:,:,1,3) = 1/(dz) * (f(:,:,:,:,2)-f(:,:,:,:,1));
    nabla_f(:,:,:,:,end,3) = 1/(dz) * (f(:,:,:,:,end)-f(:,:,:,:,end-1));
end


end