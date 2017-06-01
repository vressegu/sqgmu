classdef SVDnoise3D
%% class SVDnoise
%
% Reference implementation of the small-scale noise which relies on
% self-similarity ideas and SVD analysis.
%
% Requires: the Image Processing Toolbox for padarray() and imfilter().
%
% Written by P. DERIAN 2017-01-26
% Updated by P. DERIAN 2017-03-28: added support of 3d, 3-component fields
% Updated by P. DERIAN 2017-03-30: added mask support.
% Updated by P. DERIAN 2017-05-03: added support of 2d or 3d scalar fields.
% [TODO] mask multiplication faster than find?
% [TODO] n-d default_amplitude_scaling()
    properties(SetAccess=private)
        gridSize = [128, 128]           % the size of the grid
        patchDim = 3                    % the dimension of the (square, cubic) patch
        boundaryCondition = 'circular'  % the boundary condition        
        % these above are the main properties (and default values)
        
        nDimensions                     % the number of dimensions
        patchSize                       % the size of the patch
        numelGrid                       % the number of elements in the grid
        numelPatch                      % the number of elements in the patch
        avgFilter                       % the filter to generate the "large-scale" flow
        idxPatch                        % the patch indices
        % and these are derived properties (defined from the ones above)
    end
    
    methods
        function obj = SVDnoise3D(gridSize, patchDim, boundaryCondition)
        %% function obj = SVDnoise(gridSize, patchDim, boundaryCondition)
        %
        % The class constructor.
        %
        % Arguments:
        %   * gridSize:  [M, N] or [M, N, P] the grid dimension;
        %   * patchDim: the ODD patch dimension - typically 3;
        %   * boundayCondition: 'replicate', 'symmetric' or 'circular' -- see
        %     help padarray
        %
        % Output: a SVDnoise class instance.
        %
        % Written by P. DERIAN 2017-01-26
        % Updated by P. DERIAN 2017-03-28: added 3d grid support.

            if nargin==3
                % user parameters
                obj.gridSize = gridSize;
                obj.patchDim = patchDim;
                obj.boundaryCondition = boundaryCondition;
                % remove singleton dimensions (see squeeze)
                obj.gridSize(obj.gridSize==1) = [];
                if ~isequal(obj.gridSize, gridSize)
                    warning('SVDnoise: singleton dimensions removed in gridSize');
                end
            elseif nargin~=0
                error('SVDnoise: 0 (default) or 3 (gridSize, patchDim, boundaryCondition) parameters expected by the contructor.');                
            end
            % the grid
            obj.numelGrid = prod(obj.gridSize);
            obj.nDimensions = numel(obj.gridSize); 
            % the patch/filter
            obj.patchSize = repmat(obj.patchDim, [1, obj.nDimensions]); % compute the patch size
            obj.numelPatch = obj.patchDim^obj.nDimensions;
            obj.avgFilter = ones(obj.patchSize)/obj.numelPatch;
            obj.idxPatch = obj.patch_indices(obj.gridSize, obj.patchDim, ...
                                             obj.boundaryCondition); 
        end
        
        function [U, S, C, psobs] = vector_noise_basis(self, w, nObs, varargin)
        %% function [U, S, C] = vector_noise_basis(w, nObs, [scaling,])
        % Wrapper for 2d/3d vector_noise_basis() functions.
        % 
        % Written by P. DERIAN 2017-03-28.
            if 2==self.nDimensions
                [U, S, C, psobs] = vector_noise_basis_2d(self, w, nObs, varargin{:});
            elseif 3==self.nDimensions
                [U, S, C, psobs] = vector_noise_basis_3d(self, w, nObs, varargin{:});
            else
                error('SVDnoise:vector_noise_basis:NotYetImplemented', ...
                      'vector_noise_basis() has not been implemented for dimension %d, blame the developer.', ...
                      self.nDimensions);
            end
        end
        
        function [U, S, C, psobs] = vector_noise_basis_2d(self, w, nObs, varargin)
        %% function [U, S, C] = vector_noise_basis_2d(w, nObs, [scaling, mask])
        %
        % Arguments:
        %
        % * w: the velocity components ([M, N, 2] 3d array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        % * [optional] mask: [M, N] logical array representing regions to 
        %   exclude (true where excluded). Ignored when empty, i.e. mask=[].
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * S: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data as a [M, N, 3] array - details below;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Note: C(:,:,j) with j=1:3 contains the following data:
        % * j=1: C_xx, that is Cov(u', u');
        % * j=2: C_yy, that is Cov(v', v');
        % * j=3: C_xy, that is Cov(u', v')=Cov(v', u');
        %
        % Note: this is the older version, a.ka. "SVDfull" or "v2". 
        % The pseudo-ensemble is generated from full-field observations
        % (not fluctuations). The noise amplitude is higher, with more
        %  medium scales. It tends to accelerate events.
        %
        % Written by P. DERIAN - 2017-03-06
        % Modified by P. DERIAN - 2017-03-10: removed loop for speed-up.
        % Modified by P. DERIAN - 2017-03-30: added mask support.
        
            % check grid size
            if 2~=self.nDimensions
               error('SVDnoise:vector_noise_basis_2d:wrongDimensions', ...
                     'vector_noise_basis_2d() can only be used within a 2d context (nDimensions=%d).\n', ...
                     self.nDimensions);
            end
            % parse inputs
            [scaling, mask] = parse_inputs(varargin);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:vector_noise_basis_2d:debugMode', 'vector_noise_basis_2d() is being used in debug mode (with nObs=0), randomness disabled.\n');
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
            % so we extract corresponding values and build pseudo observations
            psobs = zeros(2*self.numelGrid, nObs);
            psobs(1:self.numelGrid,:) = w(idxrand);
            psobs(self.numelGrid+1:end,:) = w(idxrand+self.numelGrid);          
            % remove the mean along rows to get fluctuations w.r.t. ensemble mean
            psobs = bsxfun(@minus, psobs, mean(psobs,2));
            % apply the mask, if any
            % Note: if not empty, assuming dimensions match...
            if ~isempty(mask)
                % the linear index of masked points
                is_masked = find(mask);
                % the mask can be  the same for all components (self.nDimensions)
                % or defined per component (self.nDimensions + 1)
                if ndims(mask)==self.nDimensions
                    % duplicate for all 2 components
                    is_masked = [is_masked; ...
                                 is_masked+self.numelGrid];
                elseif ndims(mask)~=(self.nDimensions+1)
                    error('SVDnoise:vector_noise_basis_2d:invalidData', ...
                          'the mask number of dimensions should be self.nDimensions (%d) or self.nDimensions+1 (%d), currently %d.', ...
                          self.nDimensions, self.nDimensions+1, ndims(mask));
                end
                % and mask out fluctuations
                psobs(is_masked,:) = 0.;
            end 
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
            % upper half of c_xx_yy is c_xx, lower half is c_yy
            c_xx_yy = sum(bsxfun(@times, U.^2, S2), 2);
            c_xy = sum(bsxfun(@times, U(1:self.numelGrid,:).*U(self.numelGrid+1:end,:), S2), 2);
            C = reshape([c_xx_yy; c_xy], [self.gridSize, 3]);
            return
            
            function [s,m] = parse_inputs(v)
            % Parse optional inputs for noise_basis_2d().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                    m = [];
                elseif isscalar(v)
                    s = v{1};
                    m = [];
                elseif 2==numel(v)
                    s = v{1};
                    m = v{2};
                else
                    error('SQGMU:SVDnoise:vector_noise_basis_2d:InvalidInputs', 'Expecting at most 2 (scaling, mask) optional arguments.');
                end
            end
        end
        
        function [U, S, C, psobs] = vector_noise_basis_3d(self, w, nObs, varargin)
        %% function [U, S, C] = vector_noise_basis_3d(w, nObs, [scaling, mask])
        %
        % Arguments:
        %
        % * w: the velocity components ([M, N, P, 3] 4d array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        % * [optional] mask: [M, N, P] logical array representing regions
        %   to exclude (true where excluded). Ignored when empty, i.e. mask=[].
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * S: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data as a [M, N, P, 6] matrix - details below;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Note: C(:,:,:,j) with j=1:6 contains the following data:
        % * j=1: C_xx, that is Cov(u', u');
        % * j=2: C_yy, that is Cov(v', v');
        % * j=3: C_zz, that is Cov(w', w');
        % * j=4: C_xy, that is Cov(u', v')=Cov(v', u');
        % * j=5: C_yz, that is Cov(v', w')=Cov(w', v');
        % * j=6: C_xz, that is Cov(u', w')=Cov(w', u');
        %
        % Note: this is corresponds to the second version, a.ka. "v2". 
        % The pseudo-ensemble is generated from full-field observations
        % (not fluctuations). The noise amplitude is higher, with more
        % medium scales. It tends to accelerate events.
        %
        % Written by P. DERIAN - 2017-03-28
        % Modified by P. DERIAN - 2017-03-30: added mask support.
        
            % check grid size
            if 3~=self.nDimensions
               error('SVDnoise:vector_noise_basis_3d:wrongDimensions', ...
                     'vector_noise_basis_3d() can only be used within a 3d context (nDimensions=%d).\n', ...
                     self.nDimensions);
            end
            % parse inputs
            [scaling, mask] = parse_inputs(varargin);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:vector_noise_basis_3d:debugMode', 'vector_noise_basis_3d() is being used in debug mode (with nObs=0), randomness disabled.\n');
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
            % so we extract corresponding values and build pseudo observations
            psobs = zeros(3*self.numelGrid, nObs);
            psobs(1:self.numelGrid,:) = w(idxrand);
            psobs((self.numelGrid+1):(2*self.numelGrid),:) = w(idxrand+self.numelGrid);
            psobs((2*self.numelGrid+1):end,:) = w(idxrand+2*self.numelGrid);
            % remove the mean along rows to get fluctuations w.r.t. ensemble mean
            psobs = bsxfun(@minus, psobs, mean(psobs,2));
            % apply the mask, if any
            % Note: if not empty, assuming dimensions match...
            if ~isempty(mask)
                % the linear index of masked points
                is_masked = find(mask);
                % the mask can be  the same for all components (self.nDimensions)
                % or defined per component (self.nDimensions + 1)
                if ndims(mask)==self.nDimensions
                    % duplicate for all 3 components
                    is_masked = [is_masked; ...
                                 is_masked+self.numelGrid; ...
                                 is_masked+2*self.numelGrid];
                elseif ndims(mask)~=(self.nDimensions+1)
                    error('SVDnoise:vector_noise_basis_3d:invalidData', ...
                          'the mask number of dimensions should be self.nDimensions (%d) or self.nDimensions+1 (%d), currently %d.', ...
                          self.nDimensions, self.nDimensions+1, ndims(mask));
                end
                % and mask out fluctuations
                psobs(is_masked,:) = 0.;
            end            
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
            % upper third of c_xx_yy_zz is c_xx, middle third is c_yy, 
            % lower third is c_zz
            c_xx_yy_zz = sum(bsxfun(@times, U.^2, S2), 2);
            % upper half of c_xy_yz is c_xy, lower half is c_yz
            c_xy_yz = sum(bsxfun(@times, U(1:(2*self.numelGrid),:).*U((self.numelGrid+1):end,:), S2), 2);
            c_xz = sum(bsxfun(@times, U(1:self.numelGrid,:).*U((2*self.numelGrid+1):end,:), S2), 2);
            % concatenate and reshape
            C = reshape([c_xx_yy_zz; c_xy_yz; c_xz], [self.gridSize, 6]);
            return
            
            function [s, m] = parse_inputs(v)
            % Parse optional inputs for vector_noise_basis_3d().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                    m = [];
                elseif isscalar(v)
                    s = v{1};
                    m = [];
                elseif 2==numel(v)
                    s = v{1};
                    m = v{2};
                else
                    error('SQGMU:SVDnoise:vector_noise_basis_3d:InvalidInputs', 'Expecting at most 2 (scaling, mask) optional arguments.');
                end
            end
        end

        function [U, S, C, psobs] = vector_noise_basis_fluct(self, w, nObs, varargin)
        %% function [U, S, C] = vector_noise_basis_fluct(w, nObs, [scaling,])
        % Wrapper for 2d/3d vector_noise_basis_fluct() functions.
        % 
        % Written by P. DERIAN 2017-03-30.
            if 2==self.nDimensions
                [U, S, C, psobs] = vector_noise_basis_fluct_2d(self, w, nObs, varargin{:});
            elseif 3==self.nDimensions
                [U, S, C, psobs] = vector_noise_basis_fluct_3d(self, w, nObs, varargin{:});
            else
                error('SVDnoise:vector_noise_basis_fluct:NotYetImplemented', ...
                      'vector_noise_basis_fluct() has not been implemented for dimension %d, blame the developer.', ...
                      self.nDimensions);
            end
        end        
        
        function [U, S, C, psobs] = vector_noise_basis_fluct_2d(self, w, nObs, varargin)
        %% function [U, S, C] = vector_noise_basis_fluct_2d(w, nObs, [scaling, mask])
        %
        % Arguments:
        %
        % * w: the velocity components ([M, N, 2] 3d array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        % * [optional] mask: [M, N] logical array representing regions to 
        %   exclude (true where excluded). Ignored when empty, i.e. mask=[].
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * S: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data as a [M, N, 3] array - details below;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Note: C(:,:,j) with j=1:3 contains the following data:
        % * j=1: C_xx, that is Cov(u', u');
        % * j=2: C_yy, that is Cov(v', v');
        % * j=3: C_xy, that is Cov(u', v')=Cov(v', u');
        %
        % Note: this is the third iteration of the model, a.k.a. "v3". The
        % pseudo-ensemble is built from fluctuations computed w.r.t. some
        % spatial mean.
        %
        % Written by P. DERIAN - 2017-01-26
        % Modified by P. DERIAN - 2017-03-10: removed loop for speed-up.
        % Modified by P. DERIAN - 2017-03-30: added mask support.
        
            % check grid size
            if 2~=self.nDimensions
               error('SVDnoise:vector_noise_basis_fluct_2d:wrongDimensions', ...
                     'vector_noise_basis_fluct_2d() can only be used within a 2d context (nDimensions=%d.\n', ...
                     self.nDimensions);
            end
            % parse input
            [scaling, mask] = parse_inputs(varargin);
            % compute flutuations
            wf = w - imfilter(w, self.avgFilter, self.boundaryCondition);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:vector_noise_basis_fluct_2d:debugMode', 'vector_noise_basis_fluct_2d() is being used in debug mode (with nObs=0), randomness disabled.\n');
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
            % so we extract corresponding values and build pseudo observations
            psobs = zeros(2*self.numelGrid, nObs);
            psobs(1:self.numelGrid,:) = wf(idxrand);
            psobs(self.numelGrid+1:end,:) = wf(idxrand+self.numelGrid);            
            % apply the mask, if any
            % Note: if not empty, assuming dimensions match...
            if ~isempty(mask)
                % the linear index of masked points
                is_masked = find(mask);
                % the mask can be  the same for all components (self.nDimensions)
                % or defined per component (self.nDimensions + 1)
                if ndims(mask)==self.nDimensions
                    % duplicate for all 2 components
                    is_masked = [is_masked; ...
                                 is_masked+self.numelGrid];
                elseif ndims(mask)~=(self.nDimensions+1)
                    error('SVDnoise:vector_noise_basis_fluct_2d:invalidData', ...
                          'the mask number of dimensions should be self.nDimensions (%d) or self.nDimensions+1 (%d), currently %d.', ...
                          self.nDimensions, self.nDimensions+1, ndims(mask));
                end
                % and mask out fluctuations
                psobs(is_masked,:) = 0.;
            end 
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
            % upper half of c_xx_yy is c_xx, lower half is c_yy
            c_xx_yy = sum(bsxfun(@times, U.^2, S2), 2);
            c_xy = sum(bsxfun(@times, U(1:self.numelGrid,:).*U(self.numelGrid+1:end,:), S2), 2);
            C = reshape([c_xx_yy; c_xy], [self.gridSize, 3]);
            return
            
            function [s, m] = parse_inputs(v)
            % Parse optional inputs for vector_noise_basis_fluct_2d().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                    m = [];
                elseif isscalar(v)
                    s = v{1};
                    m = [];
                elseif 2==numel(v)
                    s = v{1};
                    m = v{2};    
                else
                    error('SQGMU:SVDnoise:vector_noise_basis_fluct_2d:InvalidInputs', 'Expecting at most 2 (scaling, mask) optional arguments.');
                end
            end
        end

        function [U, S, C, psobs] = vector_noise_basis_fluct_3d(self, w, nObs, varargin)
        %% function [U, S, C] = vector_noise_basis_fluct_3d(w, nObs, [scaling, mask])
        %
        % Arguments:
        %
        % * w: the velocity components ([M, N, P, 3] 4d array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        % * [optional] mask: [M, N, P] logical array representing regions to 
        %   exclude (true where excluded). Ignored when empty, i.e. mask=[].
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * S: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data as a [M, N, P, 6] array - details below;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Note: C(:,:,j) with j=1:3 contains the following data:
        % * j=1: C_xx, that is Cov(u', u');
        % * j=2: C_yy, that is Cov(v', v');
        % * j=3: C_zz, that is Cov(w', w');
        % * j=4: C_xy, that is Cov(u', v')=Cov(v', u');
        % * j=5: C_yz, that is Cov(v', w')=Cov(w', v');
        % * j=6: C_xz, that is Cov(u', w')=Cov(w', u');
        %
        % Note: this is the third iteration of the model, a.k.a. "v3". The
        % pseudo-ensemble is built from fluctuations computed w.r.t. some
        % spatial mean.
        %
        % Written by P. DERIAN - 2017-30-30
        
            % check grid size
            if 3~=self.nDimensions
               error('SVDnoise:vector_noise_basis_fluct_3d:wrongDimensions', ...
                     'vector_noise_basis_fluct_3d() can only be used within a 3d context (nDimensions=%d.\n', ...
                     self.nDimensions);
            end
            % parse input
            [scaling, mask] = parse_inputs(varargin);
            % compute flutuations
            wf = w - imfilter(w, self.avgFilter, self.boundaryCondition);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:vector_noise_basis_fluct_3d:debugMode', 'vector_noise_basis_fluct_3d() is being used in debug mode (with nObs=0), randomness disabled.\n');
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
            % so we extract corresponding values and build pseudo observations
            psobs = zeros(3*self.numelGrid, nObs);
            psobs(1:self.numelGrid,:) = wf(idxrand);
            psobs((self.numelGrid+1):(2*self.numelGrid),:) = wf(idxrand+self.numelGrid);
            psobs((2*self.numelGrid+1):end,:) = wf(idxrand+2*self.numelGrid);
            % apply the mask, if any
            % Note: if not empty, assuming dimensions match...
            if ~isempty(mask)
                % the linear index of masked points
                is_masked = find(mask);
                % the mask can be  the same for all components (self.nDimensions)
                % or defined per component (self.nDimensions + 1)
                if ndims(mask)==self.nDimensions
                    % duplicate for all 3 components
                    is_masked = [is_masked; ...
                                 is_masked+self.numelGrid; ...
                                 is_masked+2*self.numelGrid];
                elseif ndims(mask)~=(self.nDimensions+1)
                    error('SVDnoise:vector_noise_basis_fluct_3d:invalidData', ...
                          'the mask number of dimensions should be self.nDimensions (%d) or self.nDimensions+1 (%d), currently %d.', ...
                          self.nDimensions, self.nDimensions+1, ndims(mask));
                end
                % and mask out fluctuations
                psobs(is_masked,:) = 0.;
            end 
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
            % upper third of c_xx_yy_zz is c_xx, middle third is c_yy, 
            % lower third is c_zz
            c_xx_yy_zz = sum(bsxfun(@times, U.^2, S2), 2);
            % upper half of c_xy_yz is c_xy, lower half is c_yz
            c_xy_yz = sum(bsxfun(@times, U(1:(2*self.numelGrid),:).*U((self.numelGrid+1):end,:), S2), 2);
            c_xz = sum(bsxfun(@times, U(1:self.numelGrid,:).*U((2*self.numelGrid+1):end,:), S2), 2);
            % concatenate and reshape
            C = reshape([c_xx_yy_zz; c_xy_yz; c_xz], [self.gridSize, 6]);
            return
            
            function [s, m] = parse_inputs(v)
            % Parse optional inputs for vector_noise_basis_fluct_3d().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                    m = [];
                elseif isscalar(v)
                    s = v{1};
                    m = [];
                elseif 2==numel(v)
                    s = v{1};
                    m = v{2};    
                else
                    error('SQGMU:SVDnoise:vector_noise_basis_fluct_3d:InvalidInputs', 'Expecting at most 2 (scaling, mask) optional arguments.');
                end
            end
        end        
        
        function [U, S, C, psobs] = scalar_noise_basis(self, f, nObs, varargin)
        %% function [U, S, C] = scalar_noise_basis(f, nObs, [scaling, mask])
        %
        % Arguments:
        %
        % * f: the scalar field ([M, N] or [M, N, P] array);
        % * nObs: the number of pseudo-observations to be generated;
        % * [optional] scaling: scaling for the singular values V (and covariance
        %   C, with scaling^2);
        % * [optional] mask: [M, N] or [M, N, P] logical array representing regions to 
        %   exclude (true where excluded). Ignored when empty, i.e. mask=[].
        %
        % Outputs:
        %
        % * U: the left eigenvectors of the SVD decomposition (nobs-1 columns);
        % * S: a nobs-1 vector associated singular values;
        % * C: the one-point (co)variance data as a [M, N] or [M, N, P] array;
        % * psobs: the generated pseudo-observation (useful for debug).
        %
        % Written by P. DERIAN - 2017-05-03.
        
            % parse inputs
            [scaling, mask] = parse_inputs(varargin);
            % the random obs indices
            if 0==nObs
                % [DEBUG] all patch values (no randomness)
                warning('SVDnoise:scalar_noise_basis:debugMode', 'scalar_noise_basis() is being used in debug mode (with nObs=0), randomness disabled.\n');
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
            % so we extract corresponding values and build pseudo observations
            psobs = f(idxrand);
            % remove the mean along rows to get fluctuations w.r.t. ensemble mean
            psobs = bsxfun(@minus, psobs, mean(psobs,2));
            % apply the mask, if any
            % Note: if not empty, assuming dimensions match...
            if ~isempty(mask)
                % the mask dimensions must be same as grid
                if ndims(mask)~=self.nDimensions
                    error('SVDnoise:scalar_noise_basis:invalidData', ...
                          'the mask number of dimensions should be self.nDimensions (%d), currently %d.', ...
                          self.nDimensions, ndims(mask));
                end                
                % and mask out fluctuations
                psobs(mask(:),:) = 0.;
            end 
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
            C = reshape(sum(bsxfun(@times, U.^2, S2), 2), self.gridSize);
            return
            
            function [s,m] = parse_inputs(v)
            % Parse optional inputs for scalar_noise_basis().
            %
            % Written by P. DERIAN 2016-08-23.
                if isempty(v)
                    s = 1.;
                    m = [];
                elseif isscalar(v)
                    s = v{1};
                    m = [];
                elseif 2==numel(v)
                    s = v{1};
                    m = v{2};
                else
                    error('SQGMU:SVDnoise:scalar_noise_basis:InvalidInputs', 'Expecting at most 2 (scaling, mask) optional arguments.');
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
        %   - w, a [self.gridSize, P] array with self.gridSize the grid size
        %     and P the number of components (1 if scalar, 2 or 3 if vector); 
        %   - gamma, a vector of size(S) of the basis vector coefficients
        %     for this realization (i.e. S.*randn).   
        % Written by P. DERIAN 2017-01-26
        % Modified by P. DERIAN 2017-05-03: made compatible with scalar
        % fields.
            
            % Draw normally-distributed variables
            % and multiply by the singular values
            gamma = bsxfun(@times, randn(numel(S), 1), S); 
            % Rebuild the noise sigma.dBt.(1/dt) 
            w = reshape(U*gamma, [self.gridSize, size(U, 1)/self.numelGrid]);
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
            % [TODO] n-d version
            scale = sqrt(double(self.patchDim)^(-2./3.));
        end;
    end
    
    methods (Static)
        function Idp = patch_indices(gridSize, patchDim, boundaryCondition)
        %% Ip = patch_indices(gridSize, patchDim, boundaryCondition)
        %
        % In 2d, this function considers a [patchDim, patchDim] (square) patch sliding
        % over an image of size gridSize, with appropriate boundary conditions.
        % It returns a [prod(gridSize), patchDim^2] array where each row corresponds
        % to a grid point (patch center) and contains the global indices
        % (in a gridSize matrix) of the points in that patch. 
        % The 3d version is its natural extension in volume.
        %
        % Note: consider removing singleton dimensions in gridSiz prior to
        % calling this function for better efficiency.
        %
        % Arguments:
        %   * gridSize:  [M, N] or [M, N, P] the (volumic) image dimension;
        %   * patchDim: the ODD patch dimension;
        %   * boundayCondition: 'replicate', 'symmetric' or 'circular' -- see
        %         help padarray
        %
        % Output: Ip, the [prod(gridSize), patchDim^numel(gridSize)] array of global indices.
        %
        % Written by P. DERIAN 2016-10-25
        % Updated by P.DERIAN 2017-03-28: added 3d grid support. 

            % Check parameters
            if ~mod(patchDim, 2)
                error('patchIndices:ValueError', 'Expecting ODD patchDim (currently %d)', patchDim);
            end
            if ~strcmp(boundaryCondition, {'replicate', 'symmetric', 'circular'})
                error('patchIndices:ValueError', 'Expecting boundaryCondition to be "replicate", "symmetric" or "circular"');
            end
            % Dimensions
            patchLength = patchDim.^numel(gridSize); % number of points in a patch
            gridLength = prod(gridSize); % number of point in the entire grid
            % Create an array of global indices for a grid of size gridSize
            idx = 1:gridLength;
            % Reshape and padd with the correct boundary condition
            halfPatchDim = floor(patchDim/2);
            idxPadded = padarray(reshape(idx, gridSize), ...
                                 repmat(halfPatchDim, [1, numel(gridSize)]), ...
                                 boundaryCondition, 'both');
            % The array of patch indices
            Idp = zeros(gridLength, patchLength, 'uint32');
            % 2d case
            if 2==numel(gridSize)
                [im, in] = ind2sub(gridSize, idx);
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
                    Idp(j,:) = idxPatch(:);
                end
            % 3d case    
            elseif 3==numel(gridSize)
                [im, in, ip] = ind2sub(gridSize, idx);
                for j=idx
                    % patch indices j <=> (m, n, p)
                    m = im(j);
                    n = in(j);
                    p = ip(j);
                    % pixel index in padded array is 
                    %     np = n+halfPatchDim;
                    %     mp = m+halfPatchDim;
                    %     pp = p+halfPatchDim;
                    % so the patch bounds are:
                    %     np-halfPatchDim:np+halfPatchDim = n:n+patchDim-1
                    %     mp-halfPatchDim:mp+halfPatchDim = m:m+patchDim-1
                    %     pp-halfPatchDim:pp+halfPatchDim = p:p+patchDim-1
                    % the indices of points within the patch are:
                    idxPatch = idxPadded(m:m+patchDim-1, n:n+patchDim-1, p:p+patchDim-1);
                    Idp(j,:) = idxPatch(:);
                end
            else
                error('SVDnoise:patch_indices:NotYetImplemented', 'Only 2d and 3d grid are currently supported.');
            end
        end
    end
    
end