function fct_create_folder_plots(model)
% Create folders to save plots and files
%

folder_simu = model.folder.folder_simu;

if ~ isinf(model.sigma.k_c) % Stochastic case
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
