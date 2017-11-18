function v_proj = fct_proj_free_div(model, v)
% Project on the space of free divergence function
%

% warning('NEED TO REMOVE THIS -- TEST');
% % v(:,:,2)=2*v(:,:,1);
% v = v.^2;

%% Grid of wave vectors
% proj_free_div = model.grid.k_aa.proj_free_div;
proj_free_div = model.grid.k.proj_free_div;

%% Fourier transform
fft_v = fft2(v);

%% Projection
fft_v = sum( bsxfun( @times, proj_free_div, fft_v ) , 3);
% fft_v = permute(fft_v,[1 2 4 3]);
fft_v = permute(fft_v,[1 2 5 4 3]);

%% Inverse fourier transform
v_proj = real(ifft2(fft_v));

%% Plots for test
% figure(43);quiver(model.grid.x,model.grid.y,v(:,:,1)',v(:,:,2)');
% axis equal; axis xy; 
% figure(44);quiver(model.grid.x,model.grid.y,v_proj(:,:,1)',v_proj(:,:,2)');
% axis equal; axis xy; 
% DIV = divergence(model.grid.x,model.grid.y,v(:,:,1)',v(:,:,2)');
% CURLZ = curl(model.grid.x,model.grid.y,v(:,:,1)',v(:,:,2)');
% DIV_proj = divergence(model.grid.x,model.grid.y,v_proj(:,:,1)',v_proj(:,:,2)');
% CURLZ_proj = curl(model.grid.x,model.grid.y,v_proj(:,:,1)',v_proj(:,:,2)');
% figure(45);
% subplot(2,2,1);imagesc(model.grid.x(2:end-1),model.grid.y(2:end-1),...
%     DIV(2:end-1,2:end-1)');
% axis equal; axis xy;  colorbar;
% subplot(2,2,2);imagesc(model.grid.x(2:end-1),model.grid.y(2:end-1),...
%     DIV_proj(2:end-1,2:end-1)');
% axis equal; axis xy; colorbar;
% subplot(2,2,3);imagesc(model.grid.x(2:end-1),model.grid.y(2:end-1),...
%     CURLZ(2:end-1,2:end-1)');
% axis equal; axis xy;  colorbar;
% subplot(2,2,4);imagesc(model.grid.x(2:end-1),model.grid.y(2:end-1),...
%     CURLZ_proj(2:end-1,2:end-1)');
% axis equal; axis xy; colorbar;
% keyboard;