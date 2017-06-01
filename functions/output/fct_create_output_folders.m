function fct_create_output_folders(model)
% Create folders to save plots and files
%
% Modified by P. DERIAN 2016-08-19
%   - simpified things;
%   - improved path generation for cross-platform compatibility.
% Modified by P. DERIAN 2016-10-11
%   - made plotting optional.
% Modified by P. DERIAN 2017-05-04
%   - change function name to "fct_create_output_folder"

% Root directory
folder_simu = model.output.folder_simu;
% Subdirectories to be created
subdirs = {'files'}; %default: files only
% If plots are enabled
if model.output.plot_results
    % subdirectories for 1 realization: b and spectrum
    subdirs = [subdirs, {'one_realization', 'spectrum'}];
    % If stochastic case and plotting moments
    if model.output.plot_moments 
        % Add subdirectories for moments plots.
        subdirs = [subdirs, {'1st_2nd_order_moments', '3rd_4th_order_moments'}];
    end
end

% For each subdirectory
for k=1:length(subdirs)
    output_dir = fullfile(folder_simu, subdirs{k});
    % If does not exist, create.
    if ~ (exist(output_dir,'dir')==7)
        mkdir(output_dir);
    end
end


