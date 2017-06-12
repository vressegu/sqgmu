function d_fft_b_adv = deriv_fft_advection(model, fft_b, w, varargin)
% Compute the Fourier transform of the partial derivative along time 
% of the tracer according to the (stochastic or deterministic) transport
% equation with a hyperviscosity term
%
% Modified by P. DERIAN 2016-08-22
%   * added varargin and parse_inputs() to handle tensor A;
%   * added velocity drift correction accrodingly;
%   * changed grad(b);
% Modified by P. DERIAN 2017-06-06
%   * changed de-aliasing.

    %% Check inputs
    % mostly this is intended to check A matches the noise model
    A = parse_inputs(model, varargin);
    
    %% Grid of wave vectors
    PX = model.grid.MX/2;
    kx = model.grid.k.kx; %"normal" grid
    ky = model.grid.k.ky;
    k2 = model.grid.k.k2;
    ikx_aa = model.grid.k_aa.ikx; %anti-aliased grid for gradient
    iky_aa = model.grid.k_aa.iky;
    mask_aa = model.grid.k_aa.mask; %anti-aliasing mask

    %% Gradient of b, anti-aliased
    % in Fourier space, de-aliased, then in physical space.
    gradbx_aa = real(ifft2( ikx_aa.*fft_b ));
    gradby_aa = real(ifft2( iky_aa.*fft_b ));    
    
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
    % dealisasing of the velocity: FT, apply mask, iFT
    ft_w = fft2( w );
    w_aa = real(ifft2( bsxfun(@times, ft_w, mask_aa) ));
    % Advective term in physical space -(w* + sigma_dBt/dt).grad(b)
    % NB : in the stochastic case, w included both continuous (w*) and
    % non-continuous (sigma.dBt/dt) components of the velocity
    wgradT = -( w_aa(:,:,1).*gradbx_aa + w_aa(:,:,2).*gradby_aa );
    % Advective term in Fourier space
    adv1 = fft2( wgradT );

    %% Diffusion term (from the stochastic material derivative)
    adv2 = 0.; %default value
    % Stochastic model
    if model.is_stochastic
        if isscalar(A)
            if A>0.
                % Homogeneous and Isotropic noise: 
                %   0.5*a*Laplacian(b), with a constant.
                % Note: here, computed directly in the spectral space (hence the minus).
                adv2 = (-0.5.*A).*( k2.*fft_b );
            end;
        else
            % Non-homogeneous, anisotropic noise
            %   0.5*div(a.grad(b)), with a tensor.
            % Here we compute a.grad(b) in the physical space
            % then move to the spectral space for the div.
            % Note: we need to de-alias A (and grad(b), but it's done
            % already)
            ft_A = fft2( A );
            A_aa = real(ifft2( bsxfun(@times, ft_A, mask_aa) ));
            a_dot_gradb_1 = A_aa(:,:,1).*gradbx_aa + A_aa(:,:,3).*gradby_aa; %first component
            a_dot_gradb_2 = A_aa(:,:,3).*gradbx_aa + A_aa(:,:,2).*gradby_aa; %second component
            adv2 = 0.5*( ikx_aa.*fft2(a_dot_gradb_1) + iky_aa.*fft2(a_dot_gradb_2) );
        end
    % Or deterministic with (isotropic, homogeneous) diffusion
    elseif isscalar(A) && A>0.
        %   0.5*a*Laplacian(b), with a constant.
        % Note: here, computed directly in the spectral space (hence the minus).
        adv2 = (-0.5*A).*( k2.*fft_b );
    end

    %% Hyperviscosity
    % Hyperviscosity term
    adv4 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b;
    
    %% Summing terms
    % advection, stochastic diffusion, hyperviscosity
    d_fft_b_adv = adv1 + adv2 + adv4;
    % remove aliasing on single high-freq
    % to make sure b is real.
    d_fft_b_adv(PX(1)+1,:) = 0.;
    d_fft_b_adv(:,PX(2)+1) = 0.;
end

function A = parse_inputs(model, v)
% parse inputs for deriv_ff_advection().
%
% Written by P. DERIAN 2016-08-22.

    % defualt value
    A = 0.;
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
            % SVD: A [MX 3] matrix
            case {'SVD', 'SVDfull'}
                % make sure its size is correct
                if any(size(A)~=[model.grid.MX  3])
                    error('SQGMU:deriv_fft_advection:InvalidInputs', ...
                      'Expecting A matrix of shape [%d, %d, 3] for "SVD" noise model', ...
                      model.grid.MX(1), model.grid.MX(2));
                end
            otherwise
                warning('SQGMU:deriv_fft_advection:UnknownModel',...
                        'The model %s was not recognized - continuing hopefully with given diffusion tensor A.',...
                        model.sigma.type_noise);
        end
    % Deterministic case    
    else
        % [WIP] enabling diffusion for deterministic case
        if ~isempty(v)
    %         warning('SQGMU:deriv_fft_advection:InvalidInputs', ...
    %                 'ignoring diffusion tensor A provided with the deterministic model');
            A = v{1};
        end
    end
end
