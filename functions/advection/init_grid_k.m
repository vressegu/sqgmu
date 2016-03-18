function model = init_grid_k (model)
% Create a grid in the Fourier space
%

if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

%%  Damped Fourier grid
meth_anti_alias='deriv_LowPass';
PX=model.grid.MX/2;
if strcmp(meth_anti_alias,'deriv_LowPass')
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
        .* fct_unity_approx5(model.grid.MX(1));
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1] ...
        .* fct_unity_approx5(model.grid.MX(2));
else
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
end
[kx,ky]=ndgrid(kx,ky);
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
k2=kx.^2+ky.^2;
k2(PX(1)+1,:)=0;
k2(:,PX(2)+1)=0;
k=sqrt(k2);

% Specific operators
on_k = 1./k;
on_k ( k==0 ) = 0;

%% Save
model.grid.k.kx=kx;
model.grid.k.ky=ky;
model.grid.k.k2=k2;
model.grid.k.k=k;
model.grid.k.on_k=on_k;
clear k kx ky

%% Unstable Fourier grid for Hyper-viscosity
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
[kx,ky]=ndgrid(kx,ky);
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
k2=kx.^2+ky.^2;
k2(PX(1)+1,:)=0;
k2(:,PX(2)+1)=0;

%% Save
model.grid.k_HV.k2=k2;

end

function t = fct_unity_approx5(N_t)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

slop_size_ratio=6;

t=ones(1,N_t);
P_t=N_t/2;
sslop=ceil(N_t/slop_size_ratio);
t((P_t-sslop+1):P_t)= (-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
t((P_t+2):(P_t+1+sslop))= (tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;

t(P_t+1)=0;

end

