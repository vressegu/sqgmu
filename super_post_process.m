
%%%%%%%%%%%%%%%%%%%%
%%% Super main
%%%%%%%%%%%%%%%%%%%%
init;

%% Main parameters to choose
% Type of dynamics
dynamics = 'SQG';
% dynamics = '2D';

% Deterministic or random model
stochastic_simulation = false;
% Usual SQG model (stochastic_simulation=false)
% or SQG_MU model (stochastic_simulation=true)

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

% Viscosity
Lap_visco.bool = true;

% % Smagorinsky-like viscosity
% Smag.bool = false;
% % HV.bool = false;

% Hyper-viscosity
HV.bool = false;

% Smagorinsky-like diffusivity/viscosity or Hyper-viscosity
Smag.bool = true;

if Lap_visco.bool
    % Ratio between the Shanon resolution and filtering frequency used to
    % filter the heterogenous diffusion coefficient
    v_dealias_ratio_mask_LS = 1;
    % For Smagorinsky-like diffusivity/viscosity or Hyper-viscosity,
    % Ratio between the Shanon resolution cut-off ( = pi / sqrt( dx*dy) )
    % and the targeted diffusion scale
    v_kappamax_on_kappad = [ 0.4 0.6 1.2 ] ;
    
    for p=1:length(v_dealias_ratio_mask_LS)
        for q=1:length(v_kappamax_on_kappad)
            Smag(p,q).bool = true;
            Smag(p,q).dealias_ratio_mask_LS = v_dealias_ratio_mask_LS(p);
            Smag(p,q).kappamax_on_kappad = v_kappamax_on_kappad(q);
            Smag(p,q).weight_cst_dissip = 0;
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


%% Main
Smag = Smag(:);
%Smag = Smag(1:2);
ll=length(Smag);
for j=1:ll
%parfor j=1:ll
    plot_post_process_2(stochastic_simulation,type_data,resolution,forcing, ...
        Lap_visco,HV,Smag(j));
end

