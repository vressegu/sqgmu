function [fft_buoy_part, model] = fct_fft_advection_sto(model,  fft_buoy_part)
% Advection of buoyancy using SQG or SQG_MU model
%
% Modified by P. DERIAN 2016-08-18--?
%   - changed stochastic/deterministic tests;
%   - added support for several noise models (switches);
%   - improved verbose mode;
%   - improved output path generation for cross-platform compatibility;
%   - merged diffusion coefficient with noise models pre-computations;
%   - removed model.advection.coef_diff for homogeneization purpose;


%% Folder to save plots and files
if model.is_stochastic % Stochastic case
    model.output.folder_simu = fullfile('images', 'SQG_MU', model.type_data, model.sigma.type_noise);
else % Deterministic case  
    model.output.folder_simu = fullfile('images', 'usual_SQG', model.type_data);
end
% Create the folders
fct_create_folder_plots(model)
% Colormap
load('BuYlRd.mat');
model.output.colormap = BuYlRd; clear BuYlRd
% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.output.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Grid
% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% Grid in Fourier space
model = init_grid_k(model);

%% Initialisation of the spatial fields
% Initial large-scale velocity
fft_w = SQG_large_UQ(model, fft_buoy_part);
w = real(ifft2(fft_w));
% Ensemble size
N_ech = model.advection.N_ech;
% Create several identical realizations of the intial buoyancy
fft_buoy_part = repmat(fft_buoy_part(:,:,1),[1 1 1 N_ech]);

%% Diffusion coeff
% [TODO] move within spectral. bound for CFL "bound1"?
model.sigma.a0 = 2.*model.physical_constant.f0/model.sigma.k_c^2; 

%% Hyperviscosity
% Root mean square Lyapunov
[dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s = (dxUx+dyUy)/2;
d = (dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
clear s d dxUx dxUy dyUx dyUy lambda
% Order of the hyperviscosity
model.advection.HV.order = 8;
% Hyperviscosity coefficient
model.advection.HV.val = ...
    40 * model.advection.lambda_RMS * ...
    (mean(model.grid.dX)/pi)^model.advection.HV.order;

%% Choice of time step : CFL
% CFL of the diffusion (or CFL of the white noise advection)
dX2 = (model.grid.dX /pi).^2;
bound1 = 2./model.sigma.a0*prod(dX2)/sum(dX2);
% CFL of the (large-scale) advection
dX = permute(model.grid.dX,[1 3 2]);
bound2 = sum(bsxfun(@times,abs(w),pi./dX),3);
bound2 = max(bound2(:));
bound2 = 1./bound2/4.;
% CFL of the hyperviscosity
bound3 = 1./model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
         (model.advection.HV.order/2);
clear dX dX2
% Minimum of the CFL
dt = min([bound1 bound2 bound3]);
clear bound1 bound2 bound3
if ~ isinf(model.sigma.k_c)
    dt = dt/2.;
    % Further constraint on dt due to the use of a (simple) Euler scheme 
    % for the SPDE
end
model.advection.dt_adv = dt;

%% Pre-computations for noise models
if model.is_stochastic 
    switch model.sigma.type_noise
        case 'Spectrum'
            % Fourier transform of the kernel \tilde sigma up to a multiplicative
            % constant
            [sigma_on_sq_dt, ~, missed_var_small_scale_spectrum ] ...
                = fct_sigma_spectrum(model,fft_w);
            % sigma_on_sq_dt will be used to simulate sigma d B_t
            % missed_var_small_scale_spectrum will be used to set the mulitplicative
            % constant
            % Muliplicative constant of the kernel \tilde sigma
            model.sigma.a0_on_dt = model.sigma.a0/dt;
            sigma_on_sq_dt = sqrt(2*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
                * sigma_on_sq_dt;
            % the factor d=2 is for the dimension d of the space R^d
            % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
            clear missed_var_small_scale_spectrum
        case 'SVD'
            % Scaling factor for the basis
            % to compute a_l from large scale variance E[(u_L - E[u_L])^2]
            % Note: * sqrt(dt) is a normalization for the BM (variance of BM ~ dt);
            %       * sqrt(P^(-2/3)) is an amplitude scaling for the change of
            %           scale (observation=>noise);
            model.sigma.scaling = sqrt(model.advection.dt_adv*double(model.sigma.P)^(-2./3.));
            % Reshaping matrix
            [R, shape_field] = reshape_mat(model.grid.MX(1)/model.sigma.P, ...
                                           model.grid.MX(2)/model.sigma.P, ...
                                           model.sigma.P);
            model.sigma.R = R;
            model.sigma.shape_field = shape_field;
            clear R shape_field
            % Periodic gradient operators (in the physical space)
            % these are centered O(2) finite differences with periodic BC,
            % used to compute div(a).
            [D1, D2] = gradient_periodic(model.grid.MX(1), model.grid.MX(2), ...
                                         model.grid.dX(1), model.grid.dX(2));
            model.sigma.D1 = D1;
            model.sigma.D2 = D2;
            clear D1 D2
            % Pre-allocations
            A = zeros([model.grid.MX 3 N_ech]); % Diffusion tensor data
            sigma_dBt_dt = zeros([model.grid.MX 2 N_ech]); % noise velocity
        otherwise
            error('SQGMU:fft_advection_sto:ValueError', '"%s" is not a valid noise model', ...
                  model.sigma.type_noise);
    end
end

%% Loop over time
% Used for the first plot
tt_last = -inf;
% Number of time step
N_t = ceil(model.advection.advection_duration/dt);
% Print some information
if model.verbose
    if model.is_stochastic
        fprintf(1, '\nStochastic simulation: %s ("%s")\n', ...
                model.sigma.long_name, model.sigma.type_noise);
    else
        fprintf(1, 'Deterministic simulation\n');
    end
    fprintf(1, 'Initial condition: %s\n', model.type_data);
    fprintf(1, 'Output folders: %s\n', model.output.folder_simu);
    fprintf(1, '1/k_c is equal to %.2f m\n', 1./model.sigma.k_c);
    fprintf(1, 'Time step: %.2f seconds\n', dt);
    fprintf(1, 'Total time of advection: %f days\n', N_t*dt/3600./24.);
    fprintf(1, 'Ensemble size: %d realization(s)\n', N_ech);
    fprintf(1, 'Random generator seed: %d\n', model.seed);
    fprintf(1, '\nBeginning simulation...\n');
end
% Now iterate
for t=1:N_t   
    %% Time-correlated velocity
    fft_w = SQG_large_UQ(model, fft_buoy_part);
    w = real(ifft2(fft_w));
    clear fft_w
    
    %% Time-uncorrelated velocity 
    % [TODO] check sigma_dBt_dt means "sigma_dBt over dt"
    if model.is_stochastic % Stochastic case
        switch model.sigma.type_noise
            % Spectral noise (isotropic and homogeneous in space)
            case 'Spectrum'
                % Fourier transform of white noise
                dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
                % Multiplication by the Fourier transform of the kernel \tilde \sigma
                fft_sigma_dBt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
                clear dBt_C_on_sq_dt
                % Homogeneous velocity field
                sigma_dBt_dt = real(ifft2(fft_sigma_dBt));
                clear fft_sigma_dBt
                % Stochastic diffusion (scalar here)
                A = model.sigma.a0;
            % Pseudo-observation noise (anisotropic, inhomogeneous in space)
            case 'SVD'
                % for each realization
                parfor sampl=1:N_ech
                    % build basis for noise
                    [tmp_U, tmp_S, tmp_A, ~] = svd_noise_basis(w(:,:,1,sampl), ...
                                                               w(:,:,2,sampl), ...
                                                               model.sigma.P, ...
                                                               model.sigma.N_obs, ...
                                                               model.sigma.scaling);
                    % Scale tmp_S
                    % [TODO] check this with Valentin
                    % 1/dt scaling of tmp_S is due to implementation specificities, 
                    % see further down section "Transport of tracer".
                    % It should not be applied to tmp_A, only to tmp_S,
                    % therefore it is not included in svd_noise_basis()
                    % scaling factor.
                    tmp_S = tmp_S./model.advection.dt_adv;
                    % Draw random coefficients
                    % Note: N_obs-1 is S dimension
                    tmp_gamma = randn(model.sigma.N_obs-1, model.sigma.P^2);
                    % Multiply by the normalized eigenvalues
                    tmp_S = bsxfun(@times, tmp_gamma, tmp_S);
                    % Rebuild the noise sigma.dBt.(1/dt) 
                    tmp_sigma_dBt_dt = tmp_U*tmp_S;
                    % Reshape
                    sigma_dBt_dt(:,:,:,sampl) = reshape(model.sigma.R*tmp_sigma_dBt_dt(:), ...
                                                        model.sigma.shape_field);
                    % Interpolate A
                    % to match simulation resolution
                    A(:,:,:,sampl) = imresize(tmp_A, model.grid.MX, ...
                                              model.sigma.interp_method);
                end
                clear tmp_*
            % Error
            otherwise
                error('SQGMU:fft_advection_sto:ValueError', ...
                      '"%s" is not a valid noise model', ...
                      model.sigma.type_noise);
        end
    else % Deterministic case
        sigma_dBt_dt = 0.;
    end
    
    %% Adding time-correlated and time decorrelated velocity
    w = w + sigma_dBt_dt;
    
    %% Transport of tracer
    if model.is_stochastic
        parfor sampl=1:N_ech
            if isscalar(A)
                tmp_A = A;
            else
                tmp_A = A(:,:,:,sampl);
            end
            % Euler scheme
            % Note: dt in factor requires appropriate compensation on
            % sigma_dBt above.
            fft_buoy_part(:,:,:,sampl) = fft_buoy_part(:,:,:,sampl) ...
                + deriv_fft_advection( ...
                model, fft_buoy_part(:,:,:,sampl), w(:,:,:,sampl), tmp_A)*dt;
        end
    else
        % Runge-Kutta 4 scheme
        fft_buoy_part = RK4_fft_advection(model,fft_buoy_part, w);
    end
    
    %% Discard particles which have blown up
    iii = isnan(fft_buoy_part) | isinf(abs(fft_buoy_part));
    if any(iii(:))
        iii=any(any(any(iii,3),2),1);
        if all(iii(:))
            error('the simulation has blown up');
        end
        nb_dead_pcl = sum(iii);
        warning([ num2str(nb_dead_pcl) ' particle(s) on ' num2str(N_ech) ...
            ' have(s) blown up and' ...
            ' are(is) resampled uniformly on the set of the others particles']);
        N_ech_temp = N_ech - nb_dead_pcl;
        fft_buoy_part(:,:,:,iii)=[];
        iii_sample = randi(N_ech_temp,nb_dead_pcl,1);
        for k=1:nb_dead_pcl
            fft_buoy_part(:,:,:,end+1) = fft_buoy_part(:,:,:,iii_sample(k));
        end
    end
    clear iii
    
    %% Plots and save
    tt = floor(t*dt/(3600*24)); % Number of days
    if tt > tt_last
        tt_last = tt;
        fprintf(1, '%.0f day(s) of advection\n', t*dt/(24*3600));
        day_str = num2str(floor(t*dt/24/3600));
        % Plots
        fct_plot(model,fft_buoy_part,day_str)
        % Save file
        output_file = fullfile(model.output.folder_simu, 'files', day_str);
        % first generic variables
        save(output_file, 'model','t','fft_buoy_part','w','sigma_dBt_dt')
        % then model-specific variables
        if model.is_stochastic
            switch model.sigma.type_noise
                % Spectral noise
                case 'Spectrum'
                    save(output_file, 'sigma_on_sq_dt', '-append');
                case 'SVD'
                    save(output_file, 'A', '-append');
            end
        end
        % Print message
        if model.verbose
            fprintf(1, '\tsaved mat archive: %s\n', output_file);
        end
    end
end
end