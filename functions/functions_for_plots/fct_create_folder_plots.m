function fct_create_folder_plots(model,random_IC_large,plot_random_IC)
% Create folders to save plots and files
%
if nargin < 3
    plot_random_IC = false;
end
if nargin < 2
    random_IC_large = false;
end

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
        
        switch model.sigma.type_spectrum
            case {'SelfSim_from_LS','EOF'}
                if ~ (exist([folder_simu '/AbsDiffByScale_sigma_dB_t'],'dir')==7)
                    mkdir([folder_simu '/AbsDiffByScale_sigma_dB_t']);
                end
                if ~ (exist([folder_simu ...
                        '/AbsDiffByScale_sigma_dB_t_PostProcess'],'dir')==7)
                    mkdir([folder_simu '/AbsDiffByScale_sigma_dB_t_PostProcess']);
                end
                
                if nargin > 1
                    comp_large_IC_perturb(folder_simu,random_IC_large,plot_random_IC,...
                        model.sigma.type_spectrum);
                end
                % elseif strcmp(model.sigma.type_spectrum,'EOF')
        end
end
    function comp_large_IC_perturb(folder_simu,random_IC_large,...
            plot_random_IC,type_spectrum)
        if strcmp(type_spectrum,'EOF')
            folder_simu = [ folder_simu ...
                '/comp_EOF_SelfSim' ];
        elseif plot_random_IC
            if  random_IC_large
                folder_simu = [ folder_simu ...
                    '/large_IC_perturb' ];
            else
                folder_simu = [ folder_simu ...
                    '/small_IC_perturb' ];
            end
        else
            folder_simu = [ folder_simu ...
                '/no_IC_perturb' ];
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
