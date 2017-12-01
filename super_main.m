
%%%%%%%%%%%%%%%%%%%%
%%% Super main
%%%%%%%%%%%%%%%%%%%%
init;

%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
% dynamics = '2D';

% Type of initial condtions
type_data ='Constantin_case2' ;
% 'Vortices' : 2 large anticyclones and 2 large cyclones
%   (used in "Geophysical flow under location uncertainty", Resseguier V.,
%    Memin E., Chapron B.)
% 'Vortices2' : same as 'Vortices' but properly periodized (by Pierre Derian).
% 'Perturbed_vortices' : Same flow with slight small-scale modifications
%   (used in "Chaotic transitions and location uncertainty in geophysical
%    fluid flows", Resseguier V., Memin E., Chapron B.)
% 'Spectrum' : Gaussian random field with a spectrum slope deined by
%   the variable slop_b_ini (default value  = -5/3)
% 'Zero' : Field equal to zero everywhere

% Resolution
%resolution = 128;
%resolution = 512;
resolution = 128;
%resolution = 1024;
%resolution = 2048;

% The number of grid point is resolution^2
% It has to be an even integer

% Forcing

% Forcing or not
forcing = false;
% If yes, there is a forcing
% F = ampli_forcing * odg_b * 1/T_caract * sin( 2 freq_f pi y/L_y)
% % If yes, there is an additionnal velocity V = (0 Vy)
% % with Vy = ampli_forcing * odg_b *  sin( 2 freq_f pi y/L_y)

%% Deterministic Smag model
% Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
Smag.bool = false;

%% Stochastic terms

% Deterministic or random model
stochastic_simulation = true;
sigma.sto = stochastic_simulation;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

if sigma.sto
    % Type of spectrum for sigma dBt
    % type_spectrum = 'Band_Pass_w_Slope'; % as in GAFD part II
    %type_spectrum = 'Low_Pass_w_Slope';
    % Spectrum cst for k<km ans slope for k>km
    % type_spectrum = 'Low_Pass_streamFct_w_Slope';
    % Matern covariance for the streamfunction
    % spectrum = cst. * k2 .* ( 1 + (k/km)^2 )^slope )
    % ~ k2 for k<km ans slope for k>km
    % type_spectrum = 'BB';
    % type_spectrum = 'Bidouille';
    type_spectrum = 'SelfSim_from_LS';
    %  Sigma computed from self similarities from the large scales
    sigma.type_spectrum = type_spectrum;
    
    % Homogeneous dissipation associated with the spectrum slope
    sigma.assoc_diff = false;
    
    % Smagorinsky-like control of dissipation
    sigma.Smag.bool = true;
    
    %     % Sigma computed from self similarities from the large scales
    %     sigma.SelfSim_from_LS.bool = true;
    
    %     if sigma.SelfSim_from_LS.bool
    %         % Sigma computed from a energy of absolute diffusivity spectrum
    %         % sigma.SelfSim_from_LS.spectrum = 'energy';
    %         sigma.SelfSim_from_LS.spectrum = 'abs_diff';
    %     end
    
    % if strcmp(type_spectrum,'SelfSim_from_LS')
    % Heterrogeenosu energy flux epsilon
    sigma.hetero_energy_flux = false;
    
    % Modulation by local V L (estimated from the velocity and from
    % thegradient of the velocity)
    sigma.hetero_modulation = false;
    
    % Modulation by local V^2
    sigma.hetero_modulation_V2 = false;
    
    %     %if strcmp(type_spectrum,'SelfSim_from_LS')
    %     if sigma.hetero_modulation & strcmp(type_spectrum,'SelfSim_from_LS')
    if sigma.hetero_modulation | sigma.hetero_energy_flux ...
            | sigma.hetero_modulation_V2
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1/8;
        % Smag.dealias_ratio_mask_LS = 1/4;
        
    end
    % end
    
    % Force sigma to be diveregence free
    sigma.proj_free_div = true;
    
    if ( (sigma.Smag.bool + sigma.hetero_modulation + ...
            sigma.hetero_energy_flux + sigma.hetero_modulation_V2 ) > 1 ) ...
            || ( (sigma.Smag.bool + sigma.assoc_diff ) > 1 )
        error('These parametrizations cannot be combined');
    end
    
    if sigma.Smag.bool || sigma.assoc_diff
        % Rate between the smallest wave number of the spatially-unresolved
        % (not simulated) component of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappaMinUnresolved_on_kappaShanon = 1;
        
        % Rate between the largest wave number of the spatially-unresolved
        % (not simulated) component of sigma dBt and the largest wave
        % number of the simulation
        sigma.kappaMaxUnresolved_on_kappaShanon = 8;
        
    end
    
    %         % Factor in front of the additional constant dissipation
    %         % Set to 0 for no additional constant dissipation
    %         sigma.Smag.weight_cst_dissip = 0;
    
    % Heterogeneity of the noise
    sigma.Smag.SS_vel_homo = false;
    
    % Desactivate the noise
    sigma.no_noise = false;
    if sigma.no_noise
        warning('There isno noise here');
    end
    
    %%
    
    % Rate between the largest wave number of sigma dBt and the largest wave
    % number of the simulation
    sigma.kappamax_on_kappaShanon = 1;
    
    % Rate between the smallest and the largest wave number of sigma dBt
    if strcmp(sigma.type_spectrum , 'SelfSim_from_LS')
        if sigma.Smag.bool | ...
                Lap_visco.bool | ( HV.bool & (HV.order<=2) )
            % sigma.kappamin_on_kappamax = 1/2;
            sigma.kappamin_on_kappamax = 1/4;
            % sigma.kappamin_on_kappamax = 1/8;
            
%             v_kappamin_on_kappamax = 1 ./ [ 2 4 ];
%             
%             sigma_ref = sigma;
%             Smag_ref = Smag;
%             for r=1:length(v_kappamin_on_kappamax)
%                 Smag(r) = catstruct(Smag_ref,Smag(1));
%                 sigma(r) = catstruct(sigma_ref, sigma(1));
%                 sigma(r).kappamin_on_kappamax = v_kappamin_on_kappamax(r);
%             end
        elseif ( HV.bool & (HV.order==4) )
            sigma.kappamin_on_kappamax = 1/2;
        else
            warning('kappamin_on_kappamax may be inapropriate');
            sigma.kappamin_on_kappamax = 1/2;
            % sigma.kappamin_on_kappamax = 1/4;
            % sigma.kappamin_on_kappamax = 1/8;
        end
        
        sigma.kappaLS_on_kappamax = 1/8;
    else
        %kappamin_on_kappamax = 1/32;
        sigma.kappamin_on_kappamax = 1/2;
        % sigma.kappamin_on_kappamax = 1/128;
        %         sigma.slope_sigma = - 5;
        % warning('THIS PARAMETER NEEDS TO BE CHANGED -- TEST');
        
        sigma.kappaLS_on_kappamax = 1/8;
    end
    
    %%
    
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    if sigma.Smag.bool
        
        % Use a spatial derivation scheme for the herogeneous
        % disspation
        Smag.spatial_scheme = false;
        
        % Smagorinsky energy budget (dissipation epsilon)
        % without taking into account the noise intake
        sigma.Smag.epsi_without_noise = false;
        
        %         % Ratio between the Shanon resolution and filtering frequency used to
        %         % filter the heterogenous diffusion coefficient
        %         % Smag.dealias_ratio_mask_LS = 1;
        %         % Smag.dealias_ratio_mask_LS = 1/8;
        %         % Smag.dealias_ratio_mask_LS = 1/4;
        %         %Smag.dealias_ratio_mask_LS = 1/2;
        %         Smag.dealias_ratio_mask_LS = 1;
        %
        %         %         % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        %         %         % and the targeted diffusion scale
        %         %         % %        sigma.Smag.kappamax_on_kappad = 2;
        %         %         % sigma.Smag.kappamax_on_kappad = 1;
        %
        %         % sigma.Smag.kappamax_on_kappad = 0.5; % (better(?))
        %         % sigma.Smag.kappamax_on_kappad = 1 / 4;
        %         sigma.Smag.kappamax_on_kappad = 1 / ...
        %             sigma.kappaMaxUnresolved_on_kappaShanon;
        %%
        
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        v_dealias_ratio_mask_LS = 1./ [1 2 4 8]';
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        % and the targeted diffusion scale
        v_kappamax_on_kappad = 1./ [1 2 4 8]' ;
        
        % Rate between the smallest and the largest wave number of sigma dBt
        v_kappamin_on_kappamax = 1 ./ [ 2 4 ];
        
%         sigma_ref = sigma;
%         Smag_ref = Smag;
%         for r=1:length(v_kappamin_on_kappamax)
%             Smag(r) = catstruct(Smag_ref,Smag(1));
%             sigma(r) = catstruct(sigma_ref, sigma(1));
%             sigma(r).kappamin_on_kappamax = v_kappamin_on_kappamax(r);
%         end
        
        
        sigma_ref = sigma;
        Smag_ref = Smag;
        for p=1:length(v_dealias_ratio_mask_LS)
            for q=1:length(v_kappamax_on_kappad)
                for r=1:length(v_kappamin_on_kappamax)
                    Smag(p,q,r) = catstruct(Smag_ref,Smag(1,1,1));
                    %Smag(p,q) = Smag_ref;
                    Smag(p,q,r).dealias_ratio_mask_LS = v_dealias_ratio_mask_LS(p);
                    sigma(p,q,r) = catstruct(sigma_ref, sigma(1,1,1));
                    %sigma(p,q) = sigma_ref;
                    sigma(p,q,r).Smag.kappamax_on_kappad = v_kappamax_on_kappad(q);
                    sigma(p,q,r).kappamin_on_kappamax = v_kappamin_on_kappamax(r);
                end
            end
        end
%         sigma_ref = sigma;
%         Smag_ref = Smag;
%         for p=1:length(v_dealias_ratio_mask_LS)
%             for q=1:length(v_kappamax_on_kappad)
%                 Smag(p,q) = catstruct(Smag_ref,Smag(1,1));
%                 %Smag(p,q) = Smag_ref;
%                 Smag(p,q).dealias_ratio_mask_LS = v_dealias_ratio_mask_LS(p);
%                 sigma(p,q) = catstruct(sigma_ref, sigma(1,1));
%                 %sigma(p,q) = sigma_ref;
%                 sigma(p,q).Smag.kappamax_on_kappad = v_kappamax_on_kappad(q);
%             end
%         end
        %%
        
        
    end
    
end


%% Deterministic subgrid tensor



% Viscosity
Lap_visco.bool = false;

% % Smagorinsky-like viscosity
% Smag.bool = false;
% % HV.bool = false;

% Hyper-viscosity
HV.bool = false;


if Smag(1,1).bool
    if Lap_visco.bool
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        v_dealias_ratio_mask_LS = 1./ [1 2 4 8]';
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        % and the targeted diffusion scale
        v_kappamax_on_kappad = 1./ [1 2 4 8]' ;
        
        for p=1:length(v_dealias_ratio_mask_LS)
            for q=1:length(v_kappamax_on_kappad)
                Smag(p,q).bool = true;
                Smag(p,q).dealias_ratio_mask_LS = v_dealias_ratio_mask_LS(p);
                Smag(p,q).kappamax_on_kappad = v_kappamax_on_kappad(q);
                Smag(p,q).weight_cst_dissip = 0;
                % Use a spatial derivation scheme for the herogeneous
                % disspation
                Smag(p,q).spatial_scheme = false;
            end
        end
    elseif HV.bool
        % Ratio between the Shanon resolution and filtering frequency used to
        % filter the heterogenous diffusion coefficient
        Smag.dealias_ratio_mask_LS = 1
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
        % and the targeted diffusion scale
        Smag.kappamax_on_kappad = 1.1;% still small oscillations or just pixels?
        
        % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
        % Factor in front of the additional constant dissipation
        % Set to 0 for no additional constant dissipation
        % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
        %     % HV.weight_cst_dissip = 1/3; % still a bit of (stable) aliasing
        %     HV.weight_cst_dissip = 1/3;
        %     % HV.weight_cst_dissip = 0;
        
        Smag.weight_cst_dissip = 1/1;
        % % %     HV.weight_cst_dissip = 1/10;% bit of (stable) aliasing
    else
        Smag.kappamax_on_kappad = 0;
    end
    %     % Smag.kappamax_on_kappad = 1.1;
    %     % % Smag.kappad_on_kappamax = 1/2;
    %     if Smag.bool
    %         warning('This value needed to be tuned?')
    %     end
end

%% Main
s_Smag = size(Smag);
Smag = Smag(:);
sigma = sigma(:);
%Smag = Smag(1:2);
ll=length(Smag);
if ~ sigma(1).sto
    for j=1:ll
        sigma(j)=sigma(1);
    end
end
% for j=1:ll
parfor j=1:ll
    main(stochastic_simulation,type_data,resolution,forcing, ...
        sigma(j),Lap_visco,HV,Smag(j));
end

%% PLots

for j=1:ll
%parfor j=1:ll
    plot_post_process_2(stochastic_simulation,type_data,resolution,forcing, ...
        sigma(j),Lap_visco,HV,Smag(j));
end

%% Compared to reference
nb_days =30
resolution_HR = 1024;

error_vs_t = nan([nb_days,2,ll]);
for j=1:ll
%parfor j=1:ll
    error_vs_t(:,:,j) = post_process_error_grid(...
        stochastic_simulation,type_data,resolution,resolution_HR,...
        forcing,sigma(j),Lap_visco,HV,Smag(j));
end
error_vs_t = reshape(error_vs_t,[nb_days 2 s_Smag]);
