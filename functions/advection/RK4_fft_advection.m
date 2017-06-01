function fft_b_adv = RK4_fft_advection(model, fft_b, w, varargin)
% Time integration of the advection equation of b 
% by 4th order Runge-Kutta method with the speed w
%
% Modified by P. DERIAN 2017-05-16:
%   - added varargin to pass on optional args to deriv_fft_advection()

dt = model.advection.dt_adv;
k1 = deriv_fft_advection(model, fft_b, w, varargin{:});
k2 = deriv_fft_advection(model, fft_b + k1*dt/2., w, varargin{:});
k3 = deriv_fft_advection(model, fft_b + k2*dt/2., w, varargin{:});
k4 = deriv_fft_advection(model, fft_b + k3*dt, w, varargin{:});

fft_b_adv = fft_b + (dt/3.)*(k1/2. + k2 + k3 + k4/2.);

end

