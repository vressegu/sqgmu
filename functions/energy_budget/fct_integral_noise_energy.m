function [kern,M_tilde] = fct_integral_noise_energy(model,k1x,k1y,kpx,kpy)
% Compute the integral involved in the computation of M
%

% If there are no kp, we compute the digonal element, i.e. kpx=kpy=0
if nargin < 4
    kpx = 0;
    kpy = 0;
end

%% Useful wave vectors

iiikpx = kpx>0 ;
% Positive part of kpx
kpxp = kpx .* iiikpx;
% Negative part of kpx
kpxm = - kpx .* (~iiikpx);

iiikpy = kpy>0 ;
% Positive part of kpy
kpyp = kpy .* iiikpy;
% Negative part of kpy
kpym = - kpy .* (~iiikpy);

% Largest wave number of the simulation
k_inf = pi/sqrt(prod(model.grid.dX));

% Smallest wave number of sigma dBt
kmin = model.sigma.kappamin_on_kappamax * k_inf;

%% Integration of the covariance of the sigma dBt in Fourrier space
% To simplify the form of the domain, we consider two 2-dimensional
% integrals I1 and I2

% Integral I1
KXP = + k_inf - k1x - kpxp;
KXM = - k_inf - k1x + kpxm;
KYP = + k_inf - k1y - kpyp;
KYM = - k_inf - k1y + kpym;

I1 = +  fct_integrand_noise_energy(model, KXP , KYP) ...
     -  fct_integrand_noise_energy(model, KXP , KYM) ...
     -  fct_integrand_noise_energy(model, KXM , KYP) ...
     +  fct_integrand_noise_energy(model, KXM , KYM);
   
% Integral I2
KXP = min( kmin, KXP);
KXM = max(-kmin, KXM);
KYP = min( kmin, KYP);
KYM = max(-kmin, KYM);

I2 = +  fct_integrand_noise_energy(model, KXP , KYP) ...
     -  fct_integrand_noise_energy(model, KXP , KYM) ...
     -  fct_integrand_noise_energy(model, KXM , KYP) ...
     +  fct_integrand_noise_energy(model, KXM , KYM);
   
clear KXP KXM KYP KYM 

% Difference of the two integral
M_tilde = I1 - I2; clear I1 I2;

%% The symmetry of the matrix M and the pi/2 rotation J gives M_tilde
% the two other coefficient 

% M_tilde
M(1,2,:,:) = - M(1,2,:,:); 
M(2,1,:,:) = M(1,2,:,:);
M(2,2,:,:) = M(1,1,:,:);
% M(2,1,:,:) = M(1,2,:,:);
% M(2,2,:,:) = M(1,1,:,:);

% Kernel
kern =   kx .* (kx + kpx) .* M(1,1,:,:) ...
     + kx .* (ky + kpy) .* M(1,2,:,:) ...
     + ky .* (kx + kpx) .* M(2,1,:,:) ...
     + ky .* (ky + kpy) .* M(2,2,:,:) ;

