function [D1, D2] = gradient_periodic(M, N, varargin)
% function [D1, D2] = gradient_periodic(M, N, ...)
% Centered difference gradient with periodic BC.
%
% Arguments:
% * M, N: the domain dimensions;
% * dx [optional], dy [optional]: the spatical step. If none provided, using 1.
%   If only one is provided, dy=dx.
%
% Output: D1, D2
%   Sparses matrices such that reshape(D1*X(:), [M,N]) is the gradient of X
%   computed along the 1st dimension (i.e. vertically) with periodic BC,
%   and similarly with D2 and the 2nd dimension.
%
% Written by P. DERIAN 2016-08-22
%% parse inputs
[d1, d2] = parse_inputs(varargin); %these are the spacings along 1st and 2nd dimensions.

%% dimensions and such
shape = [M, N];
N_nz = 2*M*N; % number of non-zeros elements
i_d1 = zeros(N_nz, 1); % row-indices for D1
j_d1 = zeros(N_nz, 1); % column indices for D1
s_d1 = zeros(N_nz, 1); % weights for D1
i_d2 = zeros(N_nz, 1); % row indices for D2
j_d2 = zeros(N_nz, 1); % colum indices for D2
s_d2 = zeros(N_nz, 1); % weights for D2

%% D1
coeff_prev = -1./(2.*d1);
coeff_next = 1./(2.*d1);
k = 1;
% for each column
for n=1:N
    % BC: top border, m=1
    i = sub2ind(shape, 1, n); % column
    i_d1(k) = i;
    j_d1(k) = sub2ind(shape, M, n); % warping with periodic BC
    s_d1(k) = coeff_prev;
    k = k + 1;
    i_d1(k) = i;
    j_d1(k) = sub2ind(shape, 2, n);
    s_d1(k) = coeff_next;
    k = k + 1;
    % middle
    for m=2:(M-1)
        i = sub2ind(shape, m, n); % column
        i_d1(k) = i;
        j_d1(k) = sub2ind(shape, m-1, n);
        s_d1(k) = coeff_prev;
        k = k + 1;
        i_d1(k) = i;
        j_d1(k) = sub2ind(shape, m+1, n);
        s_d1(k) = coeff_next;
        k = k + 1;
    end
    % BC: bottom border, m=M
    i = sub2ind(shape, M, n); % column
    i_d1(k) = i;
    j_d1(k) = sub2ind(shape, M-1, n);
    s_d1(k) = coeff_prev;
    k = k + 1;
    i_d1(k) = i;
    j_d1(k) = sub2ind(shape, 1, n); % warping with periodic BC
    s_d1(k) = coeff_next;
    k = k + 1;
end


%% D2
coeff_prev = -1./(2.*d2);
coeff_next = 1./(2.*d2);
k = 1;
% for each row
for m=1:M
    % BC: left border, n=1
    i = sub2ind(shape, m, 1); % column
    i_d2(k) = i;
    j_d2(k) = sub2ind(shape, m, N); % warping with periodic BC
    s_d2(k) = coeff_prev;
    k = k + 1;
    i_d2(k) = i;
    j_d2(k) = sub2ind(shape, m, 2);
    s_d2(k) = coeff_next;
    k = k + 1;
    % middle
    for n=2:(N-1)
        i = sub2ind(shape, m, n); % column
        i_d2(k) = i;
        j_d2(k) = sub2ind(shape, m, n-1);
        s_d2(k) = coeff_prev;
        k = k + 1;
        i_d2(k) = i;
        j_d2(k) = sub2ind(shape, m, n+1);
        s_d2(k) = coeff_next;
        k = k + 1;
    end
    % BC: right border, n=N
    i = sub2ind(shape, m, N); % column
    i_d2(k) = i;
    j_d2(k) = sub2ind(shape, m, N-1);
    s_d2(k) = coeff_prev;
    k = k + 1;
    i_d2(k) = i;
    j_d2(k) = sub2ind(shape, m, 1); % warping with periodic BC
    s_d2(k) = coeff_next;
    k = k + 1;
end

%% Assemble
D1 = sparse(i_d1, j_d1, s_d1);
D2 = sparse(i_d2, j_d2, s_d2);

end

function [d1, d2] = parse_inputs(v)
    % Parse inputs for gradient_periodic().
    %
    % Written by P. DERIAN 2016-08-22
    if isempty(v) % gradient_periodic(M, N)
        d1 = 1.;
        d2 = 1.;
    elseif isscalar(v) % gradient_periodic(M, N, d)
        d1 = v{1};
        d2 = v{1};
    elseif numel(v)==2 % gradient_periodic(M, N, d1, d2)
        d1 = v{1};
        d2 = v{2};
    else
        error('MATLAB:gradient_periodic:InvalidInputs', 'invalid number of spacings - should be 0, 1 or 2 values.');
    end
end
