function [R, shape_F] = reshape_mat(M, N, P)
%% [R, shape_F] = reshape_mat(M, N, P)
%
% Arguments: M, N, P
%
%   * M: the number of rows in the subsampled space;
%   * N: the number of columns in the subsampled space;
%   * P: the patch side size.
%
% Output: R, shape_F
%
%   * R: a sparse matrix;
%   * shape_F: the shape for final reshaping
%
% The results are given by:
%   F1 = R*S(:) in 1D or F = reshape(R*S(:), shape_F) in 2D 
%
% Notes:
%
%   * fairly slow due to sub2ind()... but only called once?
%   * patch layout (see details below) must be in agreement with that of S.
%     This is a column-wise version, Matlab style.
%
% Concept:
% we have a matrix S, of size (2MN)x(PP):
%   [[1a,1b,1c,1d],
%    [2a,2b,2c,2d],
%    [3a,3b,3c,3d],
%    ...]
% which we want to reshape as F, of size (MP)x(NP)x2:
%   [[1a,1c, 5a,5c, ...],
%    [1b,1d, 5b,5d, ...],
%    [2a,bc, 6a,6c, ...],
%    [2b,bd, 6b,6d, ...],
%    ...]
% which, in turn, can instead be seen as its 1d version F1
%   [1a,1b,2a,2b, ..., 1c,1d,2c,2d, ...
%    5a,5b,6a,6b, ..., 5c,5d,6c,6d, ...]
% So we want to compute the matrix R so that F1 = R*S(:).
% R is (2MNPP)x(2MNPP), sparse.
%
% Written by P. DERIAN 2016-08-19

%% dimensions
% output matrix dimension
MN = M*N;
rdim = 2*MN*P*P; %R is rdimxdim
% input matrix shapes
shape_S = [2*MN, P*P];
shape_F = [M*P, N*P, 2];
% allocate indices arrays
% there are also rdim many non-0 elements
ir = zeros(rdim,1); %row indices
jr = zeros(rdim,1); %column indices
sr = ones(rdim,1); %this is all ones.

%% compute indices
k = 1;
for z=1:2
    for n=1:N
        for p2=1:P
            F_col = (n-1)*P + p2; %column index in F
            for m=1:M
                S_row = (z-1)*MN + (n-1)*M + m; %row index in S
                for p1=1:P
                    S_col = (p2-1)*P + p1; %column index in S
                    F_row = (m-1)*P + p1; %row index in F
                    % global indices
                    i = sub2ind(shape_F, F_row, F_col, z);
                    j = sub2ind(shape_S, S_row, S_col);
                    ir(k) = i;
                    jr(k) = j;
                    k = k+1;
                end
            end
        end
    end
end    
% assenmble matrix and return
R = sparse(ir, jr, sr, rdim, rdim);
return
end




