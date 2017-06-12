function [fft_u, fft_psi ] = SQG_large_UQ (model, fft_b)
% Compute the streamfunction and the velocity from the
% buoyancy according to the SQG_MU model
% or according to usual SQG if a_H = 0 i.e. k_c = inf
%

%% Fourier grid
if ~isfield(model.grid,'k')
    PX=model.grid.MX/2;
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1];
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    [kx,ky]=ndgrid(kx,ky);
    kx=2*pi/model.grid.dX(1)*kx;
    ky=2*pi/model.grid.dX(2)*ky;
    k=sqrt(kx.^2+ky.^2);
    k(PX(1)+1,:)=0;
    k(:,PX(2)+1)=0;
    
    %% Specific operators
    on_k = 1./k;
    on_k ( k==0 ) = 0;
else
    kx = model.grid.k.kx;
    ky = model.grid.k.ky;
    on_k = model.grid.k.on_k;
end

%% Streamfunction
fft_psi = bsxfun(@times, ...
                 on_k./model.physical_constant.buoyancy_freq_N, ...
                 fft_b) ;

%% Velocity
fft_u(:,:,1,:) = bsxfun( @times, (-1i).*ky , fft_psi ) ;
fft_u(:,:,2,:) = bsxfun( @times, 1i.*kx , fft_psi ) ;