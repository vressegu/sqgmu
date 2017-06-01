function buoy_scrambl = ens_ScrambleSc(model)
%% buoy_scrambl = ens_ScrambleSc(model)
% Generate an ensemble of low-res realizations from high-res input.
%
% Arguments:
%   - model: the current (low-res) simulation model of grid size [M, N];
% with expected model.init attributes:
%   - ratio: the high_res/low_res resolution ratio;
%   - type_data: the name of the (high-res) initial condition;
%   - circshift: better account for periodic bounday condition.
%   - scaling: scale the variance as scaling^2
%
% Output: a [M, N, 1, N_ech] array.
%
% Note: don't know why it has to be a 4-d array ... to match with velocity?
%       [TODO] ask Valentin.
%
% Written by P. DERIAN 2017-04-10
ratio = model.init.resolution_ratio;

% the high-res data: create a minimal model at appropriate resolution
model_highres = fct_physical_param();
model_highres.is_stochastic = false;
model_highres.resolution = ratio*model.resolution;
model_highres.init.type_data = model.init.type_data;
model_highres.init.type_rand = 'None';
model_highres.advection.N_ech = 1;
[fft_buoy_highres, model_highres] = fct_buoyancy_init_ens(model_highres);
buoy_highres = real(ifft2(fft_buoy_highres));

% the grids
M0 = model_highres.grid.MX(1); %highres gridsize
N0 = model_highres.grid.MX(2); %highres gridsize
MN0 = M0*N0; %highres number of points
M = model.grid.MX(1); %lowres gridsize 
N = model.grid.MX(2);
MN = M*N; %lowres number of points
% intialize indices arrays
iglob0 = reshape(1:MN0, [M0, N0]); %global indices of input grid
if model.init.circshift
    halfratio = ratio/2;
    iglob0 = circshift(iglob0, [halfratio, halfratio]);
end
% and for the output 
ratio_sqr = ratio*ratio;
patchidx = zeros([MN, ratio_sqr], 'uint32'); %holds the global indices of input grid at every point of output

% for each subsampled point,
% find the indices of candidate values
for jn=1:N
    pn = (jn-1)*ratio + 1; %patch first column
    for jm=1:M
        pm = (jm-1)*ratio + 1; %path first row
        itmp = iglob0(pm:(pm+ratio-1), pn:(pn+ratio-1)); %the patch indices
        jglob = (jn-1)*M + jm; %global index in subsampled grid
        patchidx(jglob, :) = itmp(:);
    end
end

% now scramble things up
idxrand = randi(ratio_sqr, [MN, model.advection.N_ech]); % the column of the obs at each point grid
idxrand = bsxfun(@plus, idxrand*MN, ((1-MN):0)'); %equivalent to idxrand(j,:) = self.idxPatch(j, idxrand(j,:));
idxrand = patchidx(idxrand); %these are the global indices of buoy
buoy_scrambl = buoy_highres(idxrand);

% remove the mean accross columns
buoy_mean = mean(buoy_scrambl, 2);
buoy_scrambl = bsxfun(@minus, buoy_scrambl, buoy_mean);
% and get a basis by the SVD
[Phi, Sigma, ~] = svd(buoy_scrambl, 'econ');
% rescale Sigma to get the desired variance
Sigma = diag(Sigma)*(model.init.scaling/sqrt(double(model.advection.N_ech) - 1.));

% draw realizations
gamma = randn(numel(Sigma), model.advection.N_ech);
buoy_scrambl = Phi*bsxfun(@times, gamma, Sigma);

% restore mean
buoy_scrambl = bsxfun(@plus, buoy_scrambl, buoy_mean);

% reshape
buoy_scrambl = reshape(buoy_scrambl, [M, N, 1, model.advection.N_ech]);
end

