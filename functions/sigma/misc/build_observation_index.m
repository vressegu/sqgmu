function index = build_observation_index(observation_dir, observation_pattern)
%% function index = build_observation_index(observation_dir, observation_pattern)
% Returns an index of observations file:time 
% 
% Arguments:
%   - observation_dir: the directory where observation files are found;
%   - observation_pattern: the pattern to search for.
%
% Output: a struct array with fields
%   - 'filename': the full filename;
%   - 'time': the corresponding simulation time in [s].
%
% Written by P. DERIAN 2016-09-12.

%% find observation files
observation_list = dir(fullfile(observation_dir, observation_pattern));

%% for each file
N = numel(observation_list);
fnames = cell(N,1); %holds file names
t = cell(N,1); %holds corresponding time 
for i=1:N
    % filename
    f = fullfile(observation_dir, observation_list(i).name);
    fnames{i} = f;
    % load time data
    tmp = load(f, 't', 'dt');
    t{i} = tmp.t*tmp.dt;
end

%% return a struct array
index = struct('filename', fnames, 'time', t);