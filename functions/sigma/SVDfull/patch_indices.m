function Ip = patch_indices(gridSize, patchDim, boundaryCondition)
%% Ip = patch_indices(gridSize, patchDim, boundaryCondition)
%
% This function considers a [patchDim, patchDim] (square) patch sliding
% over an image of size gridSize, with appropriate boundary conditions.
% It returns a [prod(gridSize), patchDim^2] array where each row corresponds
% to a grid point (patch center) and contains the global indices
% (in a gridSize matrix) of the points in that patch. 
%
% Arguments:
%   * gridSize:  [M, N] the image dimension;
%   * patchDim: the ODD patch dimension;
%   * boundayCondition: 'replicate', 'symmetric' or 'circular' -- see
%         help padarray
%
% Output: Ip, the [prod(gridSize), patchDim^2] array of global indices.
%
% Written by P. DERIAN 2016-10-25

% Check parameters
if ~mod(patchDim, 2)
    error('patchIndices:ValueError', 'Expecting ODD patchDim (currently %d)', patchDim);
end
if ~strcmp(boundaryCondition, {'replicate', 'symmetric', 'circular'})
    error('patchIndices:ValueError', 'Expecting boundaryCondition to be "replicate", "symmetric" or "circular"');
end

% Dimensions
patchLength = patchDim*patchDim; % number of point sin a patch
gridLength = prod(gridSize); % number of point in the entire grid
% Create an array of global indices for a grid of size gridSize
idx = 1:gridLength;
% Reshape and padd with the correct boundary condition
halfPatchDim = floor(patchDim/2);
idxPadded = padarray(reshape(idx, gridSize), [halfPatchDim,halfPatchDim], ...
                     boundaryCondition, 'both');

% The array of patch indices
[im, in] = ind2sub(gridSize, idx);
Ip = zeros(gridLength, patchLength, 'uint32');
for j=idx
    % patch indices j <=> (m, n)
    m = im(j);
    n = in(j);
    % pixel index in padded array is 
    %     np = n+halfPatchDim;
    %     mp = m+halfPatchDim;
    % so the patch bounds are:
    %     np-halfPatchDim:np+halfPatchDim = n:n+patchDim-1
    %     mp-halfPatchDim:mp+halfPatchDim = m:m+patchDim-1
    % the indices of points within the patche are:
    idxPatch = idxPadded(m:m+patchDim-1, n:n+patchDim-1);
    Ip(j,:) = idxPatch(:);
end
end
