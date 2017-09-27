function model = init_grid_k (model)
%% function model = init_grid_k (model)
% Create a grid in the Fourier space
%
% From struct "model", expecting:
%   - grid.MX (grid size);
%   - grid.dX (spatial sampling step);
%   - grid.dealias_method (de-aliasing method for pseudo-spectral codes)
%
% Modified by P. DERIAN 2017-06-05:
%   - added various dealiasing methods (model.grid.dealias_method parameter).
% Modified by P. DERIAN 2017-06-13:
%   - added ZM = PX+1 as a member of model.grid.k.

% check grid size is even
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX = model.grid.MX/2;
ZM = PX + 1; %index of the single high-freq mode to be zero'ed out.

%% "Unstable" Fourier grid
% for homogeneous diffusion Laplacian(b), hyper-viscosity Laplacian^p(b),
% SQG relationship...
nx = [ 0:(PX(1)-1) 0 (1-PX(1)):-1]; %NB: the central single high-freq is zero'ed-out
ny = [ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx = (1./model.grid.MX(1)) .* nx; 
ky = (1./model.grid.MX(2)) .* ny;
% the 2D grid
[kx,ky] = ndgrid(kx,ky);
kx = (2.*pi/model.grid.dX(1)) .* kx; %as wavenumbers
ky = (2.*pi/model.grid.dX(2)) .* ky;
k2 = kx.^2+ky.^2;
k2(ZM(1),:) = 0.; %de-alias the single high freq
k2(:,ZM(2)) = 0.;
k = sqrt(k2); %the modulus
% Specific operators
over_k = 1./k;
over_k( k==0 ) = 0.;

%%  Anti-aliased grid
% for non-linear operations, e.g. buoyancy advection w'*grad(b)
% and non-homogeneous stochastic diffusion div(a*grad(b))
if strcmp(model.grid.dealias_method, '2/3')
    % usual 2/3 rule: zero-out the 1/3 highest freqs
    maskx = double( abs(nx) < (2./3.)*PX(1) );
    masky = double( abs(ny) < (2./3.)*PX(2) );
elseif strcmp(model.grid.dealias_method, 'exp')
    % see "New Numerical Results for the Surface Quasi-Geostrophic
    % Equation", Constantin et al., J. Sci. Comput. (2012).
    alpha = 36.;
    order = 19.;
    maskx = exp(-alpha*( (2./model.grid.MX(1)).*abs(nx) ).^order);
    masky = exp(-alpha*( (2./model.grid.MX(2)).*abs(ny) ).^order);
elseif strcmp(model.grid.dealias_method, 'lowpass')
    % from legacy SQGMU code by V. Resseguier.
    maskx = fct_unity_approx5(model.grid.MX(1));
    masky = fct_unity_approx5(model.grid.MX(1));
else
    error('SQGMU:init_grid_k:invalidParameter', ...
          'Unknown de-aliasing method "%s"', model.grid.dealias_method);
end
maskxy = maskx'*masky;
maskxy(ZM(1),:) = 0.; %de-alias the single high freq
maskxy(:,ZM(2)) = 0.;

%% Save
% the "unstable" grid
model.grid.k.ZM = ZM; %indices of the single mode to force to zero
model.grid.k.kx = kx;
model.grid.k.ky = ky;
model.grid.k.k2 = k2;
model.grid.k.k = k;
model.grid.k.on_k = over_k;
model.grid.k.kx_over_ksqr = kx.*(over_k.^2);
model.grid.k.ky_over_ksqr = ky.*(over_k.^2);
% the "anti-aliased" grid
model.grid.k_aa.ikx = 1i.*maskxy.*kx; %precomputations for de-aliased gradient
model.grid.k_aa.iky = 1i.*maskxy.*ky;
model.grid.k_aa.mask = maskxy;

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

