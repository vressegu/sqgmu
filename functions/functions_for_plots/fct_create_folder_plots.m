function fct_create_folder_plots(model)
% Create folders to save plots and files
%
% Modified by P. DERIAN 2016-08-19
%   - simpified things;
%   - improved path generation for cross-platform compatibility.

% Root directory
folder_simu = model.output.folder_simu;
% Subdirectories to be created
subdirs = {'files', 'one_realization', 'Spectrum'};
if model.is_stochastic && model.output.plot_moments % Stochastic case, plotting moments
    % Add subdirectories for moments plots.
    subdirs = [subdirs, {'1st_2nd_order_moments', '3rd_4th_order_moments'}];
end

% For each subdirectory
for k=1:length(subdirs)
    output_dir = fullfile(folder_simu, subdirs{k});
    % If does not exist, create.
    if ~ (exist(output_dir,'dir')==7)
        mkdir(output_dir);
    end
end


