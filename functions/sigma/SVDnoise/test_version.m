
%% parameters
gridSize_2d = [128, 128];
gridSize_3d = [30, 20, 30];
% gridSize_3d = [3,4,5];
patchDim = 3;
boundaryCondition = 'circular';

%% check the 2D versions
sn = SVDnoise(gridSize_2d, patchDim, boundaryCondition);
sn2 = SVDnoise3D([gridSize_2d, 1], patchDim, boundaryCondition); % should warn of singleton dimension
if ~isequal(sn.idxPatch, sn2.idxPatch)
    error('patch indices of sn, sn2 should be equal..?');
end
w2_test = randn([gridSize_2d, 2]);
% the v2 basis
fprintf(1, 'testing noise_basis v2\n');
[U, S, C] = sn.noise_basis_old(w2_test, 0); % shouldn warn of debug mode
[U2, S2, C2] = sn2.noise_basis_2d(w2_test, 0); %should warn of debug mode
if ~isequal(S, S2)
    error('singular values S, S2 should be equal..?');
end
fprintf(1, '\tOK\n');
% the v3 basis
fprintf(1, 'testing noise_basis v3\n');
[U, S, C] = sn.noise_basis(w2_test, 0); % shouldn warn of debug mode
[U2, S2, C2] = sn2.noise_basis_fluct_2d(w2_test, 0); %should warn of debug mode
if ~isequal(S, S2)
    error('singular values S, S2 should be equal..?');
end
fprintf(1, '\tOK\n');

%% check the 3D versions
sn3 = SVDnoise3D(gridSize_3d, patchDim, boundaryCondition); % should warn of singleton dimension
w3_test = randn([gridSize_3d, 3]);
[U3, S3, C3] = sn3.noise_basis(w3_test, 50);