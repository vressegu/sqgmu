function [fft_buoy_part, model] = tmp_fct_fft_advection_sto(model,  fft_buoy_part)
% Advection of buoyancy using SQG or SQG_MU model
%
% Modified by P. DERIAN 2016-08-18
%   - changed stochastic/deterministic tests;
%   - added support for several noise models (switches);
%   - improved verbose mode;
%   - improved output path generation for cross-platform compatibility;
%   - merged diffusion coefficient with noise models pre-computations;
%   - removed model.advection.coef_diff for homogeneization purpose;
%   - renamed "sigma_dBt_dt" as "sigma_dBt_over_dt" for clarity.
% Modified by P. DERIAN 2016-10-11
%   - made plotting optional, enabled custom output root directory.
% Modified by P. DERIAN 2016-10-[24--28]
%   - function flow changed to boost parallelization;
%   - moved sigma_on_sq_dt to model.sigma to comply with parfor
%   - adjusting dt so that a day (86400 s) is a multiple of dt; 
%   - created a (nested) output() function for graphics and data;
%   - 

%% Folder to save plots and files
% Generate the case name (initial condition + forcing if any)
config_name = model.type_data;
if ~strcmp(model.advection.forcing.type_forcing, 'None') % Append name of forcing
    config_name = [config_name, '_', model.advection.forcing.type_forcing];
end
% Then the path, depending whether sotchastic or not
if model.is_stochastic % Stochastic case
    model.output.folder_simu = fullfile(model.output.folder_root, 'SQG_MU', ...
                                        config_name, model.sigma.type_noise);
else % Deterministic case  
    model.output.folder_simu = fullfile(model.output.folder_root, 'usual_SQG', ...
                                        config_name);
end
% Create the folders
fct_create_folder_plots(model); %[TODO] change "plots" to "output"?
% Colormap
tmp = load('BuYlRd.mat');
model.output.colormap = tmp.BuYlRd; clear tmp;
%model.output.colormap = load('BuYlRd.mat');
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
% Forcing (stored in model)
model.advection.forcing.F = forcing(model); 
% Create several identical realizations of the intial buoyancy
N_ech = model.advection.N_ech; % ensemble size
fft_buoy_part = repmat(fft_buoy_part(:,:,1),[1 1 1 N_ech]);

%% Diffusion coeff
% [TODO] move within spectral. bound for CFL "bound1"?
model.sigma.a0 = 2.*model.physical_constant.f0/model.sigma.k_c^2; 

%% Hyperviscosity
% Root mean square Lyapunov
tmp_w = w + model.advection.forcing.F;
[dxUx,dyUx] = gradient(tmp_w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(tmp_w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s = (dxUx+dyUy)/2;
d = (dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
clear s d tmp_w dxUx dxUy dyUx dyUy lambda
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
bound2 = sum(bsxfun(@times,abs(w)+abs(model.advection.forcing.amplitude),pi./dX),3);
bound2 = max(bound2(:));
bound2 = 1./bound2/4.;
% CFL of the hyperviscosity
bound3 = 1./model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
         (model.advection.HV.order/2);
clear dX dX2
% Minimum of the CFL
dt = min([bound1 bound2 bound3]);
clear bound1 bound2 bound3
if model.is_stochastic
    dt = dt/2.;
    % Further constraint on dt due to the use of a (simple) Euler scheme 
    % for the SPDE
end
% [WIP] adjust dt so that a day (86400 s) is a multiple of dt
dt = floor(dt);
while mod(86400, dt)>0 && dt>2 % brute force...
    dt = dt - 1;
end
model.advection.dt_adv = double(dt);
clear dt;
% Number of time step
N_t = ceil(model.advection.advection_duration/model.advection.dt_adv);

%% Pre-computations for noise models
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
        model.sigma.a0_on_dt = model.sigma.a0/model.advection.dt_adv;
        model.sigma.sigma_on_sq_dt = sqrt(2*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
            * sigma_on_sq_dt;
        % the factor d=2 is for the dimension d of the space R^d
        % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
        clear missed_var_small_scale_spectrum
    case 'SVD'
        % Scaling factor for the basis
        % to compute a_l from large scale variance E[(u_L - E[u_L])^2]
        % Note: * sqrt(P^(-2/3)) is an amplitude scaling for the change of
        %           scale (observation=>noise);
        %       * the time scaling is computed on the fly later.
        model.sigma.scaling = sqrt(double(model.sigma.P)^(-2./3.));
        % Reshaping matrix
        [R, shape_field] = reshape_mat(model.grid.MX(1)/model.sigma.P, ...
                                       model.grid.MX(2)/model.sigma.P, ...
                                       model.sigma.P);
        model.sigma.R = R;
        model.sigma.shape_field = shape_field;
        % Periodic gradient operators (in the physical space)
        % these are centered O(2) finite differences with periodic BC,
        % used to compute div(a).
        [D1, D2] = gradient_periodic(model.grid.MX(1), model.grid.MX(2), ...
                                     model.grid.dX(1), model.grid.dX(2));
        model.sigma.D1 = D1;
        model.sigma.D2 = D2;
        % Clean
        clear D1 D2 R shape_field
    case 'SVDfull'
        % Patch indices computed once and for all, for the correct
        % boundary condition
        model.sigma.idxpatch = patch_indices(model.grid.MX, model.sigma.P, ...
                                             model.sigma.boundary_condition); 
        % Scaling factor for the basis
        % to compute a_l from large scale variance E[(u_L - E[u_L])^2]
        % Note: * sqrt(P^(-2/3)) is an amplitude scaling for the change of
        %           scale (observation=>noise);
        %       * the time scaling is computed on the fly later.
        model.sigma.scaling = sqrt(double(model.sigma.P)^(-2./3.));
        % Periodic gradient operators (in the physical space)
        % these are centered O(2) finite differences with periodic BC,
        % used to compute div(a).
        [D1, D2] = gradient_periodic(model.grid.MX(1), model.grid.MX(2), ...
                                     model.grid.dX(1), model.grid.dX(2));
        model.sigma.D1 = D1;
        model.sigma.D2 = D2;
        % Clean
        clear D1 D2
    case 'Learning_SVD'
        if model.verbose
            fprintf('Indexing available observations...\n');
        end
        % index of observations
        obs_index_file = fullfile(model.sigma.observation_dir, ...
                                  'observations_index.mat');
        if exist(obs_index_file, 'file')
            load(obs_index_file);
        else     
            index = build_observation_index(...
                model.sigma.observation_dir, model.sigma.observation_pattern);
            save(obs_index_file, index);
        end
        model.sigma.observation_time = [index.time];
        model.sigma.observation_file = {index.filename};
        % Periodic gradient operators (in the physical space)
        % these are centered O(2) finite differences with periodic BC,
        % used to compute div(a).
        [D1, D2] = gradient_periodic(model.grid.MX(1), model.grid.MX(2), ...
                                     model.grid.dX(1), model.grid.dX(2));
        model.sigma.D1 = D1;
        model.sigma.D2 = D2;
        clear D1 D2 observation_index
    case 'None'
        % nothing to do here.
    otherwise
        error('SQGMU:fft_advection_sto:ValueError', '"%s" is not a valid noise model', ...
              model.sigma.type_noise);
end

%% Allocations for output variables
% Note: w, fft_buoy_part were created before, see "Initialisation of the spatial fields"
w = zeros([model.grid.MX 2 N_ech]); % time-correlated velocity
sigma_dBt_over_dt = zeros([model.grid.MX 2 N_ech]); % time-uncorrelated velocity
switch model.sigma.type_noise
    case {'SVD', 'SVDfull'}
        A = zeros([model.grid.MX 3 N_ech]); % Diffusion tensor data
end
% output frequency in time-step
% [TMP] half day for gisela 86400=>43200
output_frequency = ceil(86400/model.advection.dt_adv); % [WIP] the output frequency
model.output.frequency = output_frequency;

%% Display information
% Print some information
if model.verbose
    if model.is_stochastic
        fprintf(1, '\nStochastic simulation: %s ("%s")\n', ...
                model.sigma.long_name, model.sigma.type_noise);
    else
        fprintf(1, 'Deterministic simulation\n');
    end
    fprintf(1, 'Initial condition: %s\n', model.type_data);
    fprintf(1, 'Velocity forcing: %s\n', model.advection.forcing.type_forcing);
    fprintf(1, 'Output folders: %s\n', model.output.folder_simu);
    fprintf(1, 'Output frequency: %d steps (%.3f days)\n', model.output.frequency, ...
            model.output.frequency*model.advection.dt_adv/(3600*24.));
    if ~model.output.plot_results
        fprintf(1, '\t[!] Plotting disabled (-nojvm compatible mode)\n');
    end
    fprintf(1, '1/k_c is equal to %.2f m\n', 1./model.sigma.k_c);
    fprintf(1, 'Time step: %.2f seconds\n', model.advection.dt_adv);
    fprintf(1, 'Total time of advection: %f days\n', N_t*model.advection.dt_adv/3600./24.);
    fprintf(1, 'Ensemble size: %d realization(s)\n', N_ech);
    fprintf(1, 'Random generator seed: %d\n', model.seed);
    fprintf(1, '\nBeginning simulation...\n\n');
end

%% Parallel Business
try
    parpool; % create a pool of workers for parallel computation
             % (as many workers as cores available)
catch
    fprintf(1, '\n[!] unable to start parallel pool. Multiprocessing disabled.\n');
    fprintf(1, '    Consider compiling the SQGMU software (see /compile) if Distributed Computing licenses are unavailable.\n\n');
end

%% Loop over time
% gisela_daymin = 14.; % [TMP] hack for gisela
% gisela_daymax = 16.;
t = 1; % [TODO] deal with restart
while t<=N_t   
    % write the inital state
    if t==1
        output()
    end
    
%     % [TMP] frequency hack--------
%     tmp_days =  t*model.advection.dt_adv/(3600*24);
%     if (tmp_days>=gisela_daymin) && (tmp_days<gisela_daymax)
%         model.output.frequency = 100;
%     else
%         model.output.frequency = output_frequency;
%     end;
%     fprintf(1, '[!!!!!!!] day %.3f - setting freq to %d\n', tmp_days, model.output.frequency);
%     %-------------------------------
    
    % set the current time interval
    begin_t = t;
    end_t = min((t+model.output.frequency-1), N_t);
    fprintf(1 ,'%s - beginning iterations: [%d, %d]\n', datestr(datetime), begin_t, end_t);
    
    %% for each sample (s_)
    parfor sampl=1:N_ech
        % retrieve the current buoyancy
        s_fft_b = fft_buoy_part(:,:,sampl);
        % Note: defined these 2 temporary variables here just to suppress some warning
        s_sigma_dBt_over_dt = 0.;
        s_A = 0.;
        % iterate over the time interval
        for tt=begin_t:end_t
            
            %% Time-correlated velocity
            tmp_fft_w = SQG_large_UQ(model, s_fft_b);
            s_w = real(ifft2(tmp_fft_w)) + model.advection.forcing.F; %include optional forcing
            
            %% Time-uncorrelated velocity 
            % either when stochastic, or 1st step of deterministic ensemble
            % with random I.C.
            if (model.is_stochastic) || (tt==1 && ~model.is_stochastic && N_ech>1)
                switch model.sigma.type_noise
                    % Spectral noise (isotropic and homogeneous in space)
                    case 'Spectrum'
                        % Fourier transform of white noise
                        tmp_dBt_C_on_sq_dt = fft2( randn( model.grid.MX ));
                        % Multiplication by the Fourier transform of the kernel \tilde \sigma
                        tmp_fft_sigma_dBt = bsxfun(@times, ...
                                                   model.sigma.sigma_on_sq_dt, ...
                                                   tmp_dBt_C_on_sq_dt);
                        % Homogeneous velocity field
                        s_sigma_dBt_over_dt = real(ifft2(tmp_fft_sigma_dBt));
                        % Stochastic diffusion (scalar here)
                        s_A = model.sigma.a0;
                    % Pseudo-observation noise (anisotropic, inhomogeneous in space)
                    case 'SVD'
                        % Build basis for noise
                        [tmp_U, tmp_S, tmp_A, ~] = svd_noise_basis(s_w(:,:,1), ...
                                                                   s_w(:,:,2), ...
                                                                   model.sigma.P, ...
                                                                   model.sigma.N_obs, ...
                                                                   model.sigma.scaling);
                        % At that point A is the *variance* of small scale
                        % velocity, from which we compute characteristic time.
                        tmp_tau = characteristic_time(model, tmp_A);
                        % Scale A:
                        % we multiply it by tau
                        tmp_A = bsxfun(@times, tmp_A, tmp_tau);
                        % Then we interpolate it (I know...) to match
                        % simulation resolution, and we're done with A
                        s_A = imresize(tmp_A, model.grid.MX, ...
                                           model.sigma.interp_method);
                        % Scale tmp_S:
                        %   i) 1/sqrt(dt) scaling of tmp_S is due to implementation
                        % specificities, see also further down section "Transport of tracer".
                        % It should not be applied to tmp_A, only to tmp_S,
                        % therefore it is not included in svd_noise_basis()
                        % scaling factor.
                        % - the BM has a sqrt(dt) factor;
                        % - the advection terms further down have dt in factor,
                        %   so that sigma_dBt neds a (1/dt) factor to
                        %   compensate: we compute sigma_dBt_over_dt;
                        % The combined factor is then sqrt(dt)/dt = 1/sqrt(dt).
                        %   ii) there is also the sqrt(tau) factor. We already
                        % applied a tau factor to A, this is the corresponding one.
                        %   So, the overall factor is sqrt(tau/dt), which we 
                        % apply to S to minimize computations.
                        tmp_S = tmp_S.*sqrt(tmp_tau/model.advection.dt_adv);
                        % Draw random coefficients
                        % Note: N_obs is S dimension, P^2 the number of points
                        % in each patch.
                        tmp_gamma = randn(model.sigma.N_obs, model.sigma.P^2);
                        % Multiply by the normalized eigenvalues
                        tmp_S = bsxfun(@times, tmp_gamma, tmp_S);
                        % Rebuild the noise sigma.dBt.(1/dt) 
                        tmp_sigma_dBt_over_dt = tmp_U*tmp_S;
                        % Reshape for the final value
                        s_sigma_dBt_over_dt = reshape(model.sigma.R*tmp_sigma_dBt_over_dt(:), ...
                                                      model.sigma.shape_field);
                    case 'SVDfull'
                        [tmp_U, tmp_S, tmp_A] = svd_noise_basis_full(s_w(:,:,1), ...
                                                                     s_w(:,:,2), ...
                                                                     model.sigma.idxpatch, ...
                                                                     model.sigma.N_obs, ...
                                                                     model.sigma.scaling);
                        % At that point A is the *variance* of small scale
                        % velocity, from which we compute characteristic time.
                        tmp_tau = characteristic_time(model, tmp_A);
                        % Scale A:
                        % we multiply it by tau, and we're done with A.
                        s_A = bsxfun(@times, tmp_A, tmp_tau);
                        % Scale tmp_S:
                        %   i) 1/sqrt(dt) scaling of tmp_S is due to implementation
                        % specificities, see also further down section "Transport of tracer".
                        % It should not be applied to tmp_A, only to tmp_S,
                        % therefore it is not included in svd_noise_basis()
                        % scaling factor. Indeed:
                        % - the BM has a sqrt(dt) factor;
                        % - the advection terms further down have dt in factor,
                        %   so that sigma_dBt neds a (1/dt) factor to
                        %   compensate: we compute sigma_dBt_over_dt;
                        % The combined factor is then sqrt(dt)/dt = 1/sqrt(dt).
                        %   ii) there is also the sqrt(tau) factor. We already
                        % applied a tau factor to A, this is the corresponding one for S.
                        %   So, the overall factor is sqrt(tau/dt), which we 
                        % apply to S (instead of sigma_dBt_over_dt) to minimize computations.
                        tmp_S = tmp_S.*sqrt(tmp_tau/model.advection.dt_adv);
                        % Draw random coefficients
                        tmp_gamma = randn(numel(tmp_S), 1);
                        % Multiply by the normalized eigenvalues
                        tmp_S = bsxfun(@times, tmp_gamma, tmp_S);
                        % Rebuild the noise sigma.dBt.(1/dt) 
                        s_sigma_dBt_over_dt = reshape(tmp_U*tmp_S, [model.grid.MX, 2]);
                    % Noise structure learnt from high-res data
                    case 'SVDLearning'
                        % [TODO] incompatible with parfor?
                        % preload all data for a given time-step? => super
                        % heavy I guess
                    % We shouldn't go there, but just in case...
                    case 'None'
                        s_sigma_dBt_over_dt = 0.;
                    % Error                        
                    otherwise
                        error('SQGMU:fft_advection_sto:ValueError', ...
                              '"%s" is not a valid noise model', ...
                              model.sigma.type_noise);
                end
            else % Deterministic case
                s_sigma_dBt_over_dt = 0.;
                % [DEV] observe the variance
                % [TODO] incompatible with parfor?
            end
            
            %% Adding time-correlated and time decorrelated velocity
            s_wf = s_w + s_sigma_dBt_over_dt; % wf is the "full" velocity
            
            %% Transport of tracer
            if model.is_stochastic
                % Euler scheme
                % Note: dt in factor requires appropriate compensation on
                % sigma_dBt above.
                s_fft_b = s_fft_b + deriv_fft_advection( ...
                    model, s_fft_b, s_wf, s_A)*model.advection.dt_adv;
            else
                % Runge-Kutta 4 scheme
                s_fft_b = RK4_fft_advection(model, s_fft_b, s_wf);
            end
            
            %% Last step? save variables for this sample
            if tt==end_t
                % common variables
                fft_buoy_part(:,:,sampl) = s_fft_b; % buoyancy
                w(:,:,:,sampl) = s_w; % [WIP] large scale (time-correlated) velocity
                % model-dependent variables
                if model.is_stochastic
                    sigma_dBt_over_dt(:,:,:,sampl) = s_sigma_dBt_over_dt; % noise
                    switch model.sigma.type_noise
                        case {'SVD', 'SVDfull'}
                            A(:,:,:,sampl) = s_A; % diffusion tensor
                    end
                end
            end
        end
    end
    % set current global time (end of the interval)
    t = end_t;
    
    %% join: output data
    output()
    %output_gisela() % [TMP] specific hook for gisela
    
    %% move on
    t = t + 1;
end

    %% nested output function to it can access variables right away.
    function output()
        %% function output()
        %
        % Created by P. DERIAN 2016-10-25
        
        % Time information
        ndays = t*model.advection.dt_adv/(3600*24); % Number of days (float)
        % [TMP] hack
        if 1 %double(floor(ndays))==ndays
            day_str = num2str(floor(ndays)); % integer number here
            if model.verbose
                fprintf(1, '%s - iter #%d output: %.3f day(s) of advection\n', ...
                        datestr(datetime), t, ndays);
            end
            % Plots
            if model.output.plot_results
                fct_plot(model,fft_buoy_part,day_str);
            end
            % Save file
            output_file = fullfile(model.output.folder_simu, 'files', day_str);
            % first generic variables
            rstate = rng; % retrieve the random generator data
            save(output_file, 'model','t','fft_buoy_part','w','sigma_dBt_over_dt', 'rstate');
            % then model-specific variables
            if model.is_stochastic
                switch model.sigma.type_noise
                    case {'SVD', 'SVDfull', 'Learning_SVD'}
                        save(output_file, 'A', '-append');
                end
            end
            % Print message
            if model.verbose
                fprintf(1, '\tsaved mat archive: %s\n', output_file);
            end
        end
    end

    %% [TMP] special output for gisela
    function output_gisela()
        ndays = t*model.advection.dt_adv/(3600*24); % Number of days (float)
        if (ndays >= gisela_daymin) && (ndays < gisela_daymax)
            % variables
            velx = w(:,:,1,:);
            vely = w(:,:,2,:);
            x = model.grid.x;
            y = model.grid.y;
            dt = model.advection.dt_adv;
            % file name
            file_name = sprintf('step_%05d', t);
            output_file = fullfile(model.output.folder_simu, 'files', file_name);
            % save
            save(output_file, 'velx', 'vely', 'x', 'y', 't', 'dt');
            fprintf(1, '\tsaved mat archive: %s\n', output_file);
        end
    end

end

