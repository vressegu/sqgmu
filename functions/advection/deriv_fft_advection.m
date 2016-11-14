function d_fft_b_adv = deriv_fft_advection(model, fft_b, w, varargin)
% Compute the Fourier transform of the partial derivative along time 
% of the tracer according to the (stochastic or deterministic) transport
% equation with a hyperviscosity term
%
% Modified by P. DERIAN 2016-08-22
%   * added varargin and parse_inputs() to handle tensor A;
%   * added velocity drift correction accrodingly;
%   * changed grad(b);
%   * 

%% Check inputs
% mostly this is intended to check A matches the noise model
A = parse_inputs(model, varargin);

%% Grid of wave vectors
PX = model.grid.MX/2;
kx = model.grid.k.kx;
ky = model.grid.k.ky;
k2 = model.grid.k.k2;

%% Gradient of b
% Gradient of b in Fourier space
adv1x = 1i*kx.*fft_b;
adv1y = 1i*ky.*fft_b;
% Gradient of b in physical space
gradb(:,:,1) = real(ifft2(adv1x));
gradb(:,:,2) = real(ifft2(adv1y));

%% Advection term
% Drift correction of the velocity: w* = w - 0.5*(div(a))'
% with symmetric tensor a = [ [a_xx, a_xy]; [axy, ayy] ]
if ~isscalar(A)
    % Here A(:,:,1), A(:,:,2) and A(:,:,3) are a_xx, a_yy, a_xy, resp.
    % div(a)' is [dx(a_xx) + dy(a_xy); dx(a_xy) + dy(a_yy)]
    a11 = A(:,:,1);
    a22 = A(:,:,2);
    a12 = A(:,:,3);
    d1_a11 = reshape(model.sigma.D1*a11(:), model.grid.MX);
    d2_a22 = reshape(model.sigma.D2*a22(:), model.grid.MX);
    d1_a12 = reshape(model.sigma.D1*a12(:), model.grid.MX);
    d2_a12 = reshape(model.sigma.D2*a12(:), model.grid.MX);
    w(:,:,1) = w(:,:,1) - 0.5*(d1_a11 + d2_a12);
    w(:,:,2) = w(:,:,2) - 0.5*(d1_a12 + d2_a22);
    % remove the irrotational (i.e. with divergence) component
    if model.sigma.divfree_projection
        wir = helmholtz_decomp(model, w, 'irr');
        w = w - wir;
    end;
end
% Advective term in physical space -(w* + sigma_dBt/dt).grad(b)
% NB : in the stochastic case, w included both continuous (w*) and
% non-continuous (sigma.dBt/dt) components of the velocity
wgradT = -sum(bsxfun(@times,w,gradb),3);
% Advective term in Fourier space
adv1 = fft2(wgradT);
% Remove aliasing
adv1(PX(1)+1,:) = 0;
adv1(:,PX(2)+1) = 0;

%% Diffusion term (from the stochastic material derivative)
% Stochastic model
if model.is_stochastic
    if isscalar(A)
        % Homogeneous and Isotropic noise: 
        %   0.5*a*Laplacian(b), with a constant.
        % Note: here, computed directly in the spectral space (hence the minus).
        adv2 = (-0.5*A)*(k2.*fft_b);
    else
        % Non-homogeneous, anisotropic noise
        %   0.5*div(a.grad(b)), with a tensor.
        % Here we compute a.grad(b) in the physical space
        % then move to the spectral space for the div.
        % Note: a11, a22 and a12 already available from the advection term
        % computation above.
        a_dot_gradb_1 = a11.*gradb(:,:,1) + a12.*gradb(:,:,2); % first component
        a_dot_gradb_2 = a12.*gradb(:,:,1) + a22.*gradb(:,:,2); % second component
        adv2 = (0.5*1i)*( kx.*fft2(a_dot_gradb_1) + ky.*fft2(a_dot_gradb_2) );
    end
else
    adv2 = 0.;
end

%% Hyperviscosity
% Norm of wave vector
k2 = model.grid.k_HV.k2;
% Hyperviscosity term
adv4 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b;

%% Summing terms
% advection, stochastic diffusion, hyperviscosity
d_fft_b_adv = adv1 + adv2 + adv4;
end

function A = parse_inputs(model, v)
% parse inputs for deriv_ff_advection().
%
% Written by P. DERIAN 2016-08-22.

% If the chosen model requires an A
if model.is_stochastic
    % check there is something
    if isempty(v)
        error('SQGMU:deriv_fft_advection:InvalidInputs', ...
              'Expecting variance tensor/scalar A with stochastic model');
    end
    A = v{1};
    % make sure it is the right type
    switch model.sigma.type_noise
        % Spectral noise: constant scalar A=a0
        case 'Spectrum'
            if ~isscalar(A)
                 error('SQGMU:deriv_fft_advection:InvalidInputs', ...
                       'Expecting scalar A with "Spectrum" noise model');
            end
        % SVD: A [Mx Mx 3] matrix
        case {'SVD', 'Learning_SVD'}
            % make sure its size is correct
            if any(size(A)~=[model.grid.MX(1) model.grid.MX(2) 3])
                error('SQGMU:deriv_fft_advection:InvalidInputs', ...
                  'Expecting A matrix of shape [%d, %d, 3] for "SVD" noise model', ...
                  model.grid.MX(1), model.grid.MX(2));
            end
    end
% Deterministic case    
else
    % there should not be an a
    if ~isempty(v)
        warning('SQGMU:deriv_fft_advection:InvalidInputs', ...
                'ignoring A provided with the deterministic model');
    end
    A = 0;
end
end
