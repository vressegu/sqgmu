classdef SVDnoise
%% class SVDnoise
%
% Reference implementation of the small-scale noise which relies on
% self-similarity ideas and SVD analysis.
%
% Requires: the Image Processing Toolbox for padarray() and imfilter().
%
% Written by P. DERIAN 2017-01-26

    properties(SetAccess=private)
        % main properties
        gridSize = [128, 128]           % the size of the grid
        patchDim = 3                    % the dimension of the (square) patch
        boundaryCondition = 'circular'  % the boundary condition        
        % derived properties (defined by the main ones)
        numelGrid                       % the number of elements in the grid
        numelPatch                      % the number of elements in the patch
        avgFilter                       % the filter to generate the "large-scale" flow
        idxPatch                        % the patch indices
    end
    
    methods
        function obj = SVDnoise(gridSize, patchDim, boundaryCondition)
        %% function obj = SVDnoise(gridSize, patchDim, boundaryCondition)
        %
        % The class constructor.
        %
        % Arguments:
        %   * gridSize:  [M, N] the grid dimension;
        %   * patchDim: the ODD patch dimension;
        %   * boundayCondition: 'replicate', 'symmetric' or 'circular' -- see
        %     help padarray
        %
        % Output: a SVDnoise class instance.
        %
        % Written by P. DERIAN 2017-01-26

            if nargin==3
                % user parameters
                obj.gridSize = gridSize;
                obj.patchDim = patchDim;
                obj.boundaryCondition = boundaryCondition;
            elseif nargin~=0
                error('SVDnoise: 0 (default) or 3 (gridSize, patchDim, boundaryCondition) parameters expected by the contructor.');                
            end
            obj.numelGrid = prod(obj.gridSize);
            obj.numelPatch = obj.patchDim^2;
            obj.avgFilter = ones(obj.patchDim)/obj.numelPatch;
            obj.idxPatch = obj.patch_indices(obj.gridSize, obj.patchDim, ...
                                             obj.boundaryCondition); 
        end

        function [U, S, C, psobs] = noise_basis(self, w, nObs, varargin)
        %% function [U, S, C] = noise_basis(w, nObs, [scaling,])
        %
        % Arguments:
        %
        % * w: the velocity components ([M, N, 2] 3D array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * V: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Note: this is the newer version, a.k.a. "v3".
        %
        % Written by P. DERIAN - 2017-01-26
        % Modified by P. DERIAN - 2017-03-10: removed loop for speed-up.
            % parse input
            scaling = parse_inputs(varargin);
            % compute flutuations
            wf = w - imfilter(w, self.avgFilter, self.boundaryCondition);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:noise_basis:debugMode', 'noise_basis() is being used in debug mode (with nObs=0), randomness disabled.\n');
                nObs = self.numelPatch;
                idxrand = repmat(1:self.numelPatch, [self.numelGrid, 1]);
            else
                % draw, for each numelGrid patch, nObs indices between 1 and numelPatch,
                % i.e. draw one point within the patch for each pseudo-observation of each patch.
                idxrand = randi(self.numelPatch, [self.numelGrid, nObs]);
            end
            % transform these as global coordinates in self.idxPatch
            % and from these, get the coordinates in the field.
            % Note: this is equivalent to
            %     idxrand(j,:) = self.idxPatch(j, idxrand(j,:));
            idxrand = bsxfun(@plus, idxrand*self.numelGrid, ((1-self.numelGrid):0)');
            idxrand = self.idxPatch(idxrand);
            % so we extract corrersponding values and build pseudo observations
            psobs = zeros(2*self.numelGrid, nObs);
            psobs(1:self.numelGrid,:) = wf(idxrand);
            psobs(self.numelGrid+1:end,:) = wf(idxrand+self.numelGrid);            
            % computing the SVD: the "economy" version, i.e. the nObs first columns.
            [U, S, ~] = svd(psobs, 'econ');
            S = diag(S);
            % normalize S to get the proper variance
            % Note: the covariance matrix is (F.F^T)/(nobs-1)
            % i.e. (U.S.S^T.U^T)/(nobs-1) (with S a diagonal matrix)
            % so we spread the 1/(nobs-1) factor over Sigma
            S = S.*(scaling/sqrt(double(nObs) - 1.));
            % computing the variance matrices
            S2 = (S.^2)';
            % upper half of a_xxyy is a_xx, lower half is a_yy
            c_xxyy = sum(bsxfun(@times, U.^2, S2), 2);
            c_xy = sum(bsxfun(@times, U(1:self.numelGrid,:).*U(self.numelGrid+1:end,:), S2), 2);
            C = reshape([c_xxyy; c_xy], [self.gridSize, 3]);
            return
            
            function s = parse_inputs(v)
            % Parse optional inputs for noise_basis().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                elseif isscalar(v)
                    s = v{1};
                else
                    error('SQGMU:SVDnoise:noise_basis:InvalidInputs', 'Expecting at most 1 (scaling) optional arguments.');
                end
            end
        end

        function [U, S, C, psobs] = noise_basis_old(self, w, nObs, varargin)
        %% function [U, S, C] = noise_basis_old(w, nObs, [scaling,])
        %
        % Arguments:
        %
        % * w: the velocity components ([M, N, 2] 3D array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * V: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Note: this is the older version, a.ka. "SVDfull" or "v2". 
        % The pseudo-ensemble is generated from full-field observations
        % (not fluctuations). The noise amplitude is higher, with more
        %  medium scales. It tends to accelerate events.
        %
        % Written by P. DERIAN - 2017-03-06
        % Modified by P. DERIAN - 2017-03-10: removed loop for speed-up.
        
            scaling = parse_inputs(varargin);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:noise_basis_old:debugMode', 'noise_basis_old() is being used in debug mode (with nObs=0), randomness disabled.\n');
                nObs = self.numelPatch;
                idxrand = repmat(1:self.numelPatch, [self.numelGrid, 1]);
            else
                % draw, for each numelGrid patch, nObs indices between 1 and numelPatch,
                % i.e. draw one point within the patch for each pseudo-observation of each patch.
                idxrand = randi(self.numelPatch, [self.numelGrid, nObs]);
            end
            % transform these as global coordinates in self.idxPatch
            % and from these, get the coordinates in the field.
            % Note: this is equivalent to
            %     idxrand(j,:) = self.idxPatch(j, idxrand(j,:));
            idxrand = bsxfun(@plus, idxrand*self.numelGrid, ((1-self.numelGrid):0)');
            idxrand = self.idxPatch(idxrand);
            % so we extract corrersponding values and build pseudo observations
            psobs = zeros(2*self.numelGrid, nObs);
            psobs(1:self.numelGrid,:) = w(idxrand);
            psobs(self.numelGrid+1:end,:) = w(idxrand+self.numelGrid);          
            % remove the mean along rows to get fluctuations w.r.t. ensemble mean
            psobs = bsxfun(@minus, psobs, mean(psobs,2));
            % computing the SVD: the "economy" version, i.e. the nObs first columns.
            [U, S, ~] = svd(psobs, 'econ');
            S = diag(S);
            % normalize S to get the proper variance
            % Note: the covariance matrix is (F.F^T)/(nobs-1)
            % i.e. (U.S.S^T.U^T)/(nobs-1) (with S a diagonal matrix)
            % so we spread the 1/(nobs-1) factor over Sigma
            S = S.*(scaling/sqrt(double(nObs) - 1.));
            % computing the variance matrices
            S2 = (S.^2)';
            % upper half of a_xxyy is a_xx, lower half is a_yy
            c_xxyy = sum(bsxfun(@times, U.^2, S2), 2);
            c_xy = sum(bsxfun(@times, U(1:self.numelGrid,:).*U(self.numelGrid+1:end,:), S2), 2);
            C = reshape([c_xxyy; c_xy], [self.gridSize, 3]);
            return
            
            function s = parse_inputs(v)
            % Parse optional inputs for noise_basis_old().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                elseif isscalar(v)
                    s = v{1};
                else
                    error('SQGMU:SVDnoise:noise_basis_old:InvalidInputs', 'Expecting at most 1 (scaling) optional arguments.');
                end
            end
        end
        
        function [w, gamma] = draw_noise(self, U, S)
        %% function w = draw_noise(U, S)
        %
        % Draw a realization of a small-scale noise from the given
        % structure.
        %
        % Arguments: U, S the basis vectors and associated singular values
        % as returned by SVDnoise.noise_basis().
        %
        % Output:
        %   - w, a [M, N, 2] array with [M, N] the grid size; 
        %   - gamma, a vector of size(S) of the basis vector coefficients
        %     for this realization (i.e. S.*randn).   
        % Written by P. DERIAN 2017-01-26
            
            % Draw normally-distributed variables
            % and multiply by the singular values
            gamma = bsxfun(@times, randn(numel(S), 1), S); 
            % Rebuild the noise sigma.dBt.(1/dt) 
            w = reshape(U*gamma, [self.gridSize, 2]);
        end
        
        function wl = large_scale(self, w)
        %% function wl = get_large_scale(w)
        % Returns the "large-scale" version of field w.
        %
        % Written by P. DERIAN 2017-02-05.
            wl = imfilter(w, self.avgFilter, self.boundaryCondition);
        end;
        
        function [wf, wl] = fine_scale(self, w)
        %% function wf = get_fine_scale(w)
        % Returns the "fine-scale" version of field w.
        %
        % Written by P. DERIAN 2017-02-05.
            wl = imfilter(w, self.avgFilter, self.boundaryCondition);
            wf = w - wl;
        end;
        
        function scale = default_amplitude_scaling(self)
        %% scale = default_amplitude_scaling()
        % Return the default amplitude scaling for the velocity, function
        % of the patch dimension, to be used e.g. in noise_basis().
        %
        % Written by P. DERIAN 2017-03-06.
            scale = sqrt(double(self.patchDim)^(-2./3.));
        end;
    end
    
    methods (Static)
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
        patchLength = patchDim*patchDim; % number of points in a patch
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
            % the indices of points within the patch are:
            idxPatch = idxPadded(m:m+patchDim-1, n:n+patchDim-1);
            Ip(j,:) = idxPatch(:);
        end
        end
    end
    
end