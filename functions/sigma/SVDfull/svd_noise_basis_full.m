function [U, S, C] = svd_noise_basis_full(u, v, idxpatch, nobs, varargin)
%% function [U, S, C] = svd_noise_basis_full(u, v, idxpatch, nobs, [varargin,])
%
% Arguments:
%
% * u, v: the velocity components (2D arrays);
% * idxpatch: the indices of patch values -- see path_indices();
% * nobs: the number of pseudo-observations to generate.
% * [optional] scaling: scaling for the singular values V (and covariance
%   C, with scaling^2);
% Outputs:
%
% * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
% * V: a nobs-1 vector associated singular values;
% * C: the one-point (co)variance data;
%
% Written by P. DERIAN - 2016-10-20

%% parse inputs
scaling = parse_inputs(varargin);

%% parameters
% check dimensions match
if ~sameshape(u, v)
    error('noiseBasis:RuntimeError', 'u and v shapes do not match');
end;
if size(idxpatch,1)~=numel(u)
    error('noiseBasis:RuntimeError', 'u,v and idxpatch shapes are not compatible');
end

%% dimensions
[ydim_f, xdim_f] = size(u);
length_f = size(idxpatch, 1); %the number of points in the original image -- should be numel(u)
length_p = size(idxpatch, 2); %the number of points in a single patch

%% pseudo observations
psobs = zeros(2*length_f, nobs);
% draw, for each length_f patch, nobs indices between 1 and length_p,
% i.e. draw one point within the patch for each observation of each patch.
idxrand = randi(length_p, [length_f, nobs]);
% for each row (i.e. each patch)
for j=1:length_f
    % pick u, v values at the random indices
    % note: this means we always use (u, v) values at a same location.
    ir = idxpatch(j, idxrand(j,:)); % here we draw nobs indices amongst the patch indices
    %ir = idxpatch(j, :); % [DEBUG] all patch vlaues, requires nobs=size(idxpatch, 1)
    psobs(j,:) = u(ir); % and extract corresponding values
    psobs(j+length_f,:) = v(ir);
    % note: here we fill as upper half u, lower half v
end
% remove the mean along rows to get fluctuations w.r.t. ensemble mean
psobs = bsxfun(@minus, psobs, mean(psobs,2));

%% computing the SVD
% compute the "economy" version, i.e. the nobs first columns.
[U, S, ~] = svd(psobs, 'econ');
% normalize S to get the proper variance
% Note: the covariance matrix is (F.F^T)/(nobs-1)
% i.e. (U.S.S^T.U^T)/(nobs-1) (with S a diagonal matrix)
% so we spread the 1/(nobs-1) factor over Sigma
S = diag(S).*(scaling/sqrt(double(nobs) - 1.));

%% computing the variance matrices
S2 = (S.^2)';
% upper half of a_xxyy is a_xx, lower half is a_yy
c_xxyy = sum(bsxfun(@times, U.^2, S2), 2);
c_xy = sum(bsxfun(@times, U(1:length_f,:).*U(length_f+1:end,:), S2), 2);
C = reshape([c_xxyy; c_xy], [ydim_f, xdim_f, 3]);
return
end

function s = parse_inputs(v)
% Parse optional inputs for svd_noise_basis().
%
% Written by P. DERIAN 2016-08-23.
if isempty(v)
    s = 1.;
elseif isscalar(v)
    s = v{1};
else
    error('SQGMU:svd_noise_basis_full:InvalidInputs', 'Expecting at most 1 (scaling) optional arguments.');
end
end

function result = sameshape(A, B)
%% function result = sameshape(A, B)
%
% Return true if both A and B "arrays" have the same shape.
%
% Arguments:
%
% * A, B arrays to be checked for same shape.
%
% See https://fr.mathworks.com/matlabcentral/answers/25317-array-dimensions-equivalence-check
% Written by P. DERIAN 2016-08-17.
result = isequal(size(A), size(B)) || (isvector(A) && isvector(B) && numel(A) == numel(B));
end
