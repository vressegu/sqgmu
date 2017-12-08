function fct_create_folder_plots(model,random_IC_large)
% Create folders to save plots and files
%

folder_simu = model.folder.folder_simu;

if model.sigma.sto % Stochastic case
    if ~ (exist([folder_simu '/1st_2nd_order_moments'],'dir')==7)
        mkdir([folder_simu '/1st_2nd_order_moments']);
    end
    if ~ (exist([folder_simu '/3rd_4th_order_moments'],'dir')==7)
        mkdir([folder_simu '/3rd_4th_order_moments']);
    end
end
if ~ (exist([folder_simu '/files'],'dir')==7)
    mkdir([folder_simu '/files']);
end
if ~ (exist([folder_simu '/one_realization'],'dir')==7)
    mkdir([folder_simu '/one_realization']);
end
if ~ (exist([folder_simu '/Spectrum'],'dir')==7)
    mkdir([folder_simu '/Spectrum']);
end
if ~ (exist([folder_simu '/Dissipation'],'dir')==7)
    mkdir([folder_simu '/Dissipation']);
end
if ~ (exist([folder_simu '/dissip_coef'],'dir')==7)
    mkdir([folder_simu '/dissip_coef']);
end
if ~ (exist([folder_simu '/Epsilon_k'],'dir')==7)
    mkdir([folder_simu '/Epsilon_k']);
end
if ~ (exist([folder_simu '/Epsilon_k_meth'],'dir')==7)
    mkdir([folder_simu '/Epsilon_k_meth']);
end
if model.sigma.sto ... % Stochastic case
        & strcmp(model.sigma.type_spectrum,'SelfSim_from_LS')
    if ~ (exist([folder_simu '/AbsDiffByScale_sigma_dB_t'],'dir')==7)
        mkdir([folder_simu '/AbsDiffByScale_sigma_dB_t']);
    end
    if ~ (exist([folder_simu ...
            '/AbsDiffByScale_sigma_dB_t_PostProcess'],'dir')==7)
        mkdir([folder_simu '/AbsDiffByScale_sigma_dB_t_PostProcess']);
    end
end
if nargin > 1
    comp_large_IC_perturb(folder_simu,random_IC_large);
end

    function comp_large_IC_perturb(folder_simu,random_IC_large)
        if  random_IC_large
            folder_simu = [ folder_simu ...
                '/large_IC_perturb' ];
        else
            folder_simu = [ folder_simu ...
                '/small_IC_perturb' ];
        end
        if ~ (exist([folder_simu '/spatial_error'],'dir')==7)
            mkdir([folder_simu '/spatial_error']);
        end
        if ~ (exist([folder_simu '/Estim_spatial_bias'],'dir')==7)
            mkdir([folder_simu '/Estim_spatial_bias']);
        end
        if ~ (exist([folder_simu '/spatial_bias'],'dir')==7)
            mkdir([folder_simu '/spatial_bias']);
        end
        if ~ (exist([folder_simu '/Estim_spectral_error'],'dir')==7)
            mkdir([folder_simu '/Estim_spectral_error']);
        end
    end
end
