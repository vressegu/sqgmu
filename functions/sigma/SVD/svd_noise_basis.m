function [U, S, A, shape] = svd_noise_basis(u, v, sizep, nobs, varargin)
%% function [U, S, A, shape] = svd_noise_basis(u, v, sizep, nobs [, scaling])
%
% Arguments:
%
% * u, v: the velocity components (2D arrays);
% * psize: the size of a patch (psize x psize in 2D);
% * nobs: the number of pseudo-observations to generate.
% * [optional] scaling=1.: if provided, scales the singular values S by "scaling"
%       and therefore the (co)variance A by "scaling^2". 
%
% Outputs:
%
% * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
% * S: a nobs-1 vector associated singular values;
% * A: the one-point (co)variance data;
% * shape: the shape of the subsampled space.
%
% Written by P. DERIAN - 2016-08-18

%% parse inputs
scaling = parse_inputs(varargin);

%% dimensions
% check dimensions match
if ~sameshape(u, v)
    ME = MException('noiseBasis:RuntimeError', 'u and v shapes do not match');
    throw(ME)
end;
[ydim, xdim] = size(u);
% compute patch size and subsampled dimensions
length_p = sizep*sizep; %the number of points in a single patch
ydim_s = ydim/sizep;
xdim_s = xdim/sizep;
length_s = ydim_s*xdim_s; %the number of points in the subsampled image

%% the observations
obs = zeros(2*length_s, length_p);
% note: proceed along columns (oppa Matlab style)
for n=1:xdim_s
    % column indices of the patch
    nmin = (n-1)*sizep + 1;
    nmax = nmin + sizep - 1;
    for m=1:ydim_s
        % row indices of the patch
        mmin = (m-1)*sizep + 1;
        mmax = mmin + sizep - 1;
        % the row of observations
        j = (n-1)*ydim_s + m;
        % store U, V values
        tmp = u(mmin:mmax, nmin:nmax);
        obs(j,:) = tmp(:);
        tmp = v(mmin:mmax, nmin:nmax);
        obs(j+length_s,:) = tmp(:);
   end
end
% At this point "obs" contains, in each row, all the values observed within
% a patch, upper half of columns are from U and the lower half V.
% Patches are listed in a column-wise order.

%% generating pseudo observations
psobs = zeros(2*length_s, nobs);
% for each row (i.e. each patch)
for j=1:length_s
    % draw nobs indices between 1 and length_p, i.e. draw one point within 
    % the patch for each observation of each patch.
    itmp = randi(length_p, 1, nobs);
    % pick u, v values at these indices
    % note: this means we always use (u, v) values at a same location.
    psobs(j,:) = obs(j, itmp);
    psobs(j+length_s,:) = obs(j+length_s, itmp);
    % note: here we fill as upper half u, lower half v
end
% remove the mean along the rows to get fluctuations
psobs = bsxfun(@minus, psobs, mean(psobs,2));

%% computing the SVD
% compute the "economy" version, i.e. the nobs first columns.
[U, S, ~] = svd(psobs, 'econ');
% S is of rank nobs-1 since we removed the mean
% so we leave out the last value/column
% [TODO] check if it slows down, otherwise we might as well let it be.
S = diag(S);
S = S(1:end-1);
U = U(:,1:end-1);
% Normalize S to get the proper variance
% Note: the covariance matrix is (F.F^T)/nobs
% i.e. (U.S.S^T.U^T)/nobs (with S a diagonal matrix)
% so we spread the 1/nobs factor over Sigma.
% Here we also add the optional scaling. 
S = S.*(scaling/sqrt(double(nobs)));
% finally write down the shape of the subsampled space
shape = [ydim_s, xdim_s, 2];

%% computing the variance matrices
S2 = (S.^2)';
% upper half of a_xxyy is a_xx, lower half is a_yy
a_xxyy = sum(bsxfun(@times, U.^2, S2), 2);
a_xy = sum(bsxfun(@times, U(1:length_s,:).*U(length_s+1:end,:), S2), 2);
A = reshape([a_xxyy; a_xy], [ydim_s, xdim_s, 3]);

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
    error('SQGMU:svd_noise_basis:InvalidInputs', 'Optional scaling parameter should be a scalar.');
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
