function [sigma_on_sq_dt,f_sigma,trace_a,....
    slope_w_a_estim,mult_offset_spectrum_a_estim,...
    km_LS, ...
    spectrum_a_sigma] = ...
    fct_sigma_spectrum_abs_diff_mat(model,ft_w,bool_plot,day)
% - sigma_on_sq_dt is the Fourier transform of the kernel \tilde sigma up to a multiplicative
% constant
% - f_sigma is the Fourier transform of the associted streamfunction
% - trace_a_on_dt measures the total energy of the field and is generally used to
% set the muliplicative constant
% - spectrum_a_sigma is the spectrum
%

% Average over realizations
ft_w2=abs(ft_w).^2;
% ft_w2=mean(abs(ft_w).^2,4);

% Fourier transform norm
ft_w2=sum(ft_w2,3);
% if nargin >2
%     ft2=mean(abs(ft2).^2,4);
%     ft2=sum(ft2,3);
% end

% switch model.dynamics
%     case 'SQG'
%         slope_ref = -5/3;
%     case '2D'
%         slope_ref = -3;
%     otherwise
%         error('Unknown type of dynamics');
% end
slope_ref = model.sigma.slope_sigma_ref;
% slope_ref = model.sigma.slope_sigma;
slope_ref_a = (slope_ref-3)/2;

% Get parameters
MX=model.grid.MX;
dX=model.grid.dX;
if any(size(ft_w2)~=[ MX 1 model.advection.N_ech] )
    error('wrong size');
end
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=MX/2;
ft_w2(PX(1),:,:)=0;
ft_w2(:,PX(2),:,:)=0;

%% Wave vector
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
[kx,ky]=ndgrid(kx,ky);
k=sqrt(kx.^2+ky.^2);
k(PX(1)+1,:)=inf;
k(:,PX(2)+1)=inf;
k=k(:);

%% Wave number
M_kappa=min(model.grid.MX);
P_kappa= M_kappa/2;
d_kappa = 2*pi/sqrt(prod(model.grid.MX.* model.grid.dX));
kappa= d_kappa * ( 0:(P_kappa-1) ) ;
% M_kappa=min(MX);
% P_kappa= M_kappa/2;
% %d_kappa = max(1./dX);
% kappa=1/(M_kappa)* (0:(P_kappa-1)) ;
% kappa=2*pi*max(1./dX)*kappa;
% %kappa=2*pi*d_kappa*kappa;

%% Masks associated with the rings of iso wave number
d_kappa = kappa(2) - kappa(1);
kappa_shift = [ kappa(2:end) kappa(end)+d_kappa ];
% kappa_rep = repmat(kappa,[1 1 1 model.advection.N_ech ]);
% kappa_shift_rep = repmat(kappa_shift,[1 1 1 model.advection.N_ech ]);
% k_rep = repmat(k,[1 1 1 model.advection.N_ech ]);
% idx = sparse( bsxfun(@le,kappa_rep, k_rep ) );
% idx = idx & sparse( bsxfun(@lt,k_rep, kappa_shift_rep ) );
idx = sparse( bsxfun(@le,kappa, k ) );
% idx = idx & sparse( bsxfun(@lt,k, [ kappa(2:end) kappa(end)+d_kappa ] ) );
idx = idx & sparse( bsxfun(@lt,k, kappa_shift ) );

kappa = kappa';

%% Spectrum
% Integration over the rings of iso wave number
ft_w2 = reshape(ft_w2 , [prod(model.grid.MX) model.advection.N_ech ]);
spectrum_w = idx' * ft_w2; %  [ P_kappa model.advection.N_ech ] 

% ft_w2 = reshape(ft_w2 , [prod(model.grid.MX) 1 1 model.advection.N_ech ]);
% % idx = repmat (idx , [1 1 1 model.advection.N_ech ]);
% spectrum_w = sum( idx .* ft_w2 , 1);
% % spectrum_w = sum( bsxfun( @times, idx, ft_w2 ) , 1);
% spectrum_w = multitrans(spectrum_w); %  [ P_kappa 1 1 model.advection.N_ech ] 
% % spectrum_w = idx' * ft_w2(:);
% % % if nargin>2
% % %     spectrum_w2 = idx' * ft2(:);
% % % end

% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% velocity averaged (not just integrated) over the space
spectrum_w = 1/prod(model.grid.MX)^2 * spectrum_w;
% if nargin>2
%     spectrum_w2 = 1/prod(model.grid.MX)^2 * spectrum_w2;
% end

% Division by the wave number step
d_kappa = kappa(2)-kappa(1);
spectrum_w = spectrum_w / d_kappa;
% if nargin>2
%     spectrum_w2 = spectrum_w2 / d_kappa;
% end

% % From buoyancy spectrum to velocity spectrum
% spectrum_w = spectrum_w / model.physical_constant.buoyancy_freq_N^2;

% From velocity spectrum to absolute diffusity by scale
spectrum_w_a = zeros(size(spectrum_w));
spectrum_w_a(2:end,:) = bsxfun( @times, kappa(2:end).^(-3/2) , ...
                                     spectrum_w(2:end,:).^(1/2) );
                                 %  [ P_kappa model.advection.N_ech ] 
% spectrum_w_a = zeros(size(kappa));
% spectrum_w_a(2:end) = kappa(2:end).^(-3/2) .* spectrum_w(2:end).^(1/2);

% Time scales
% %tau_k = 1./sqrt( d_kappa * kappa.^2 .* spectrum_w);
% v_k = sqrt( d_kappa * spectrum_w) * sqrt(prod(model.grid.MX));
% tau_k = 1./sqrt( kappa.^3 .* spectrum_w);
% figure(56);loglog(kappa(2:end),v_k(2:end)/3600/24);
% figure(55);loglog(kappa(2:end),tau_k(2:end)/3600/24);

% Compensated absolute diffusity by scale
spectrum_w_a_comp = bsxfun( @times, kappa.^(-slope_ref_a) , ...
                                     spectrum_w_a );
                                 %  [ P_kappa 1 1 model.advection.N_ech ] 
% spectrum_w_a_comp = kappa.^(-slope_ref_a) .* spectrum_w_a;

%% Absosute diffusivity
abs_diff_w = sum(spectrum_w_a,1) * d_kappa;
                                 %  [ 1 model.advection.N_ech ] 

% %% Time integral of the dissipation per scale epsilon(k)
% int_epsilon = - cumsum(spectrum) * d_kappa;

%% Test
% fft_grad_v(:,:,:,1) = fct_grad(model,ft_w(:,:,1,1));
% fft_grad_v(:,:,:,2) = fct_grad(model,ft_w(:,:,2,1));
% fft_grad_v = permute(fft_grad_v,[1 2 4 3]);
% grad_v = real(ifft2(fft_grad_v));
% n_grad_v = sqrt(sum(sum(grad_v.^2,4),3));
% v = real(ifft2(ft_w));
% n_v2 = sum(v.^2,3);
%
% % mask_aa_LS = model.grid.k_aa.mask; %anti-aliasing mask
% mask_aa_LS = model.grid.k_aa_LS.mask; %anti-aliasing mask
%
% % Filtering the diffusivity coefficient at large scales
% n_v2 = fft2(n_v2 );
% n_v2 = real(ifft2( bsxfun(@times, n_v2, mask_aa_LS) ));
%
% % Filtering the diffusivity coefficient at large scales
% n_grad_v = fft2(n_grad_v );
% n_grad_v = real(ifft2( bsxfun(@times, n_grad_v, mask_aa_LS) ));
%
%
% abs_diff_estim_local = n_v2 ./ (n_grad_v/sqrt(2)); %
% abs_diff_w_estim_global1 = mean(abs_diff_estim_local(:))
% abs_diff_w_estim_global2 = mean(n_v2(:))/mean(n_grad_v(:)/sqrt(2))
%
% coef_a = 1/abs_diff_w * abs_diff_estim_local ;
%
% % warning('ft_w2 can be replaced by ft_w');
% % warning('This treatement should be done realisation by realisation');
%
% figure;imagesc(abs_diff_estim_local');
% figure;imagesc(n_v2');
% figure;imagesc(n_grad_v');
% figure;imagesc(coef_a');

%% Estimation of slope of the spectral absolute diffusivity of
% % the large-scale velocity

% Get the large-scale "length scale" km_LS and the spectral window for the
% linear regression
if ~ model.sigma.sto
    model.sigma.kappamin_on_kappamax = 1/2;
    model.sigma.kappaLS_on_kappamax = 1/8;
end
threshold_k = model.sigma.kappamin_on_kappamax;
threshold_k_LS = model.sigma.kappaLS_on_kappamax;
spectrum_w_a_comp_cut = spectrum_w_a_comp( kappa < kappa(end)*threshold_k , :);
[~,i_first] = max(spectrum_w_a_comp_cut( [ false; ...
    ( kappa(2:end) < kappa(end)*threshold_k_LS )] , :) ,[], 1 );
i_first = i_first +1;
mask_iii_k_LS = (1:size(spectrum_w_a_comp_cut,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS = bsxfun( @ge,  mask_iii_k_LS , i_first );
% % iii_k_LS = [ ones(1,i_first-1) i_first:length(spectrum_w_a_comp_cut) ];
% iii_k_LS = i_first:length(spectrum_w_a_comp_cut);
spectrum_w_a_comp_for_estim = mask_iii_k_LS .* spectrum_w_a_comp_cut;
% spectrum_w_a_comp_for_estim = spectrum_w_a_comp_cut(iii_k_LS);
kappa_cut = kappa( kappa<kappa(end)*threshold_k );
kkk = bsxfun( @times, mask_iii_k_LS , kappa_cut );
% kkk=kappa_cut(iii_k_LS);
mask_iii_k_LS_one_value = (1:size(spectrum_w_a_comp_cut,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS_one_value = bsxfun( @eq,  mask_iii_k_LS_one_value , i_first );
km_LS = sum( bsxfun( @times,mask_iii_k_LS_one_value , kappa_cut ) ,1);
% km_LS=kkk(1);
% % offset_w_a = spectrum_w_a(iii_k_LS(1));
offset_w_a_comp = sum( mask_iii_k_LS_one_value .* spectrum_w_a_comp_cut ,1);
% offset_w_a_comp = spectrum_w_a_comp(iii_k_LS(1));

% Linear regression

% %threshold_k = 1/6;
% iii_reliable =~ bsxfun( @or , ...
%     (isinf(abs(spectrum_w_a_comp_for_estim))|...
%     isnan(spectrum_w_a_comp_for_estim)|...
%     (spectrum_w_a_comp_for_estim/max(spectrum_w_a_comp_for_estim(:)) ...
%     <=eps)) , ...
%     (kkk>kappa(end)*threshold_k) );
iii_reliable =~ (isinf(abs(spectrum_w_a_comp_for_estim))|...
    isnan(spectrum_w_a_comp_for_estim)|...
    (spectrum_w_a_comp_for_estim/max(spectrum_w_a_comp_for_estim(:)) ...
    <=eps)|...
    (kkk>kappa(end)*threshold_k));
%     (kkk'>kidx(end)/3));
%     (kkk'>10*km));
% logspectrum_w_a_comp_for_estim_centered = bsxfun( @plus, ...
%     log10(spectrum_w_a_comp_for_estim.* iii_reliable) .* iii_reliable , ...
%     - log10( offset_w_a_comp) );
% % logspectrum_w_a_comp_for_estim_centered = ...
% %     log10(spectrum_w_a_comp_for_estim(iii_reliable))...
% %     - log10(offset_w_a_comp);
% logkkk = bsxfun( @plus, ...
%     log10(kkk.* iii_reliable) .* iii_reliable , ...
%     - log10( km_LS ) );
% % logkkk = log10(kkk(iii_reliable)) - log10(km_LS);

% Centering in the point ( km_LS, spectrum_w_a_comp(km_LS) )
logspectrum_w_a_comp_for_estim_centered = bsxfun( @plus, ...
    log10(spectrum_w_a_comp_for_estim) , - log10( offset_w_a_comp) );
logkkk_centered = bsxfun( @plus, log10(kkk)  , - log10( km_LS ) );

% Set to zero the unwanted points 
siz = size(iii_reliable);
iii_reliable = iii_reliable(:);
logspectrum_w_a_comp_for_estim_centered( ~ iii_reliable ) = 0; % P_kappa_cut*N_ech
logkkk_centered( ~ iii_reliable ) = 0; % P_kappa*N_ech
logspectrum_w_a_comp_for_estim_centered = reshape ( ...
    logspectrum_w_a_comp_for_estim_centered, siz ); % P_kappa_cut x N_ech
logkkk_centered = reshape ( ...
    logkkk_centered, siz ); % P_kappa_cut x N_ech

% Linear regression
beta_num = sum( ...
    logkkk_centered .* logspectrum_w_a_comp_for_estim_centered , 1);
beta_den = sum( ...
    logkkk_centered .* logkkk_centered , 1);
beta = beta_num ./ beta_den; % 1 x N_ech
% beta =min([zeros([1 model.advection.N_ech]) ; beta]); 

% beta = logkkk_centered \ logspectrum_w_a_comp_for_estim_centered;
% beta =min([0 beta]); % Prevent unstable behavior of this parametrisation.
slope_w_a_estim = beta + slope_ref_a;

% Prevent some possible unstable behaviors of this parametrisation:
% The maximum allowed slope corresponds to a velocity white in space
slope_w_a_estim = min( [zeros([1 model.advection.N_ech]) ; ...
    slope_w_a_estim]); 

slope_w_estim = 2 * slope_w_a_estim + 3;
%offset = beta(2);
clear beta


%% 1D Spectrum

% Largest wave number
k_inf = kappa(min(PX));

% Smallest wave number
k0 = model.sigma.kappamin_on_kappamax * k_inf;

% % Slope of the absolute diffusitiy distribution by scale
% model.sigma.slope_absDif_sigma = (model.sigma.slope_sigma-3)/2;

% Absolute diffusivity by scale with estimated slope
reference_spectrum_a_estim = bsxfun(@power, kappa(2:end), ...
                                          slope_w_a_estim ) ;
%reference_spectrum_a_estim = kappa(2:end) .^ slope_w_a_estim ;

% Absolute diffusivity by scale with theoretical slope
if bool_plot
    reference_spectrum_a = kappa(2:end) .^ slope_ref_a ;
    reference_spectrum_a = repmat( reference_spectrum_a , ...
        [ 1 model.advection.N_ech ]);
end
% reference_spectrum_a = bsxfun(@power, kappa(2:end) , slope_ref_a );
% %reference_spectrum_a = kappa(2:end) .^ model.sigma.slope_sigma ;

%% Offsets
mask_iii_k_LS_one_value_long = (1:size(spectrum_w_a,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS_one_value_long = bsxfun( @eq,  ...
    mask_iii_k_LS_one_value_long , i_first );
spectrum_w_a_km_LS = ...
    sum( bsxfun( @times,mask_iii_k_LS_one_value_long , spectrum_w_a ) ,1);

mask_iii_k_LS_one_value_long_m1 = (1:size(reference_spectrum_a_estim,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS_one_value_long_m1 = bsxfun( @eq,  ...
    mask_iii_k_LS_one_value_long_m1 , i_first - 1 );
reference_spectrum_a_estim_km_LS = ...
    sum( bsxfun( @times,mask_iii_k_LS_one_value_long_m1 , ...
    reference_spectrum_a_estim ) ,1);
if bool_plot
    reference_spectrum_a_km_LS = ...
        sum( bsxfun( @times,mask_iii_k_LS_one_value_long_m1 , ...
        reference_spectrum_a ) ,1);
end

% Apply offset
mult_offset_spectrum_a_estim = spectrum_w_a_km_LS ...
    / reference_spectrum_a_estim_km_LS;
% mult_offset_spectrum_a_estim = spectrum_w_a(iii_k_LS(1)) ...
%     / reference_spectrum_a_estim(iii_k_LS(1)-1);
% %10^offset
% % reference_spectrum_a_estim = 10^offset...
% %     * reference_spectrum_a_estim;
reference_spectrum_a_estim = bsxfun(@times, ...
    mult_offset_spectrum_a_estim , reference_spectrum_a_estim );

if bool_plot
    mult_offset_spectrum_a = spectrum_w_a_km_LS ...
        / reference_spectrum_a_km_LS ;
    % mult_offset_spectrum_a = spectrum_w_a_km_LS ...
    %     / reference_spectrum_a(iii_k_LS(1)-1) ;
    reference_spectrum_a = bsxfun(@times, ...
        mult_offset_spectrum_a , reference_spectrum_a );
end

% % Compute the offset to superimpose the slop -5/3 and the spectrum of w
% idx_not_inf=~(isinf(log10(spectrum_w_a(2:end)))| ...
%     spectrum_w_a(2:end)<1e-4*max(spectrum_w_a(2:end)) | isinf(kappa(2:end)'));
% idx_not_inf= [ false; idx_not_inf ];
% idx_not_inf(end)=false;
% % reference_spectrum = reference_spectrum * 10 .^( ...
% %     mean( log10(spectrum_w(idx_not_inf)') ...
% %     - log10( reference_spectrum(idx_not_inf)) ) );
% mult_offset = 10 .^( ...
%     mean( log10(spectrum_w_a(idx_not_inf)') ...
%     - log10( reference_spectrum_a(idx_not_inf)) ) );
% %reference_spectrum = reference_spectrum * mult_offset;

% Residual absolute diffusivity by scale
residual_spectrum_a = ...
    reference_spectrum_a_estim - spectrum_w_a(2:end,:);

% Clean it
siz = size(residual_spectrum_a);
residual_spectrum_a = residual_spectrum_a(:);
residual_spectrum_a(residual_spectrum_a<0)=0;
residual_spectrum_a = reshape( residual_spectrum_a, siz);

mask_iii_k_LS_long_m1 = (1:size(reference_spectrum_a_estim,1))' ...
    * ones(1,model.advection.N_ech);
mask_iii_k_LS_long_m1 = bsxfun( @gt,  ...
    mask_iii_k_LS_long_m1 , i_first - 1 );
residual_spectrum_a = bsxfun( @times,mask_iii_k_LS_long_m1 , ...
                                                residual_spectrum_a);

% residual_spectrum_a(1:(iii_k_LS(1)-1))=0;
residual_spectrum_a(1:2,:)=0;

%% Spectre of random small-scale velocity
f_sigma = residual_spectrum_a;
%f_sigma = f_sigma * d_kappa;

% switch model.sigma.type_spectrum
%     case 'Band_Pass_w_Slope'
%         f_sigma = reference_spectrum;
%
%         % % Parseval ( * prod(model.grid.MX) )
%         % % and integrated over the space ( * prod(model.grid.MX) )
%         % f_sigma = prod(model.grid.MX)^2 * f_sigma;
%
% Band-pass filter
idx1 = (kappa(2:end) <= k0);
idx3 = (kappa(2:end) > k_inf);
idx2 = ~ (idx1 | idx3);

unit_approx = fct_unity_approx_(sum(idx2));
% unit_approx = ones([1,sum(idx2)]);
% warning('Above line modified for debug !!!!!')

f_sigma(idx1 | idx3)=0;
f_sigma(idx2) = bsxfun(@times, unit_approx' , f_sigma(idx2,:) ) ;

%     case 'Low_Pass_w_Slope'
%         f_sigma = (k0^2 + kappa(2:end).^2) .^ (model.sigma.slope_sigma/2) ;
%     case  'Low_Pass_streamFct_w_Slope'
%         f_sigma = kappa(2:end).^2 .* ...
%             (k0^2 + kappa(2:end).^2) .^ (model.sigma.slope_sigma/2-1) ;
%     case  'BB'
%         f_sigma = ones(size(kappa(2:end))) ;
%     case  'Bidouille'
%         f_sigma = 1/10 * ones(size(kappa(2:end))) ;
%     otherwise
%         error('Unknown spectrum type for the small-scale velocity');
% end

% To remove the 2 pi which appear when we integrate the spectrum over the
% wave-vector angles and add the (2 pi)^2 which appear when we go from k to
% 2*pi*k
% And discretisation to go from continuous to discrete Fourier transform
% f_sigma = f_sigma * (2*pi);
% f_sigma = 1/prod(model.grid.dX) * f_sigma;
f_sigma = (2*pi/prod(model.grid.dX)) * f_sigma;

% Division by prod(model.grid.MX) because of the variance white-in-space
% noise
% Multiplication by prod(model.grid.MX) because the spectrum respresent
% an energy by unit of space and not the sum over the space
% f_sigma = 1/prod(model.grid.MX) * f_sigma;
% f_sigma = prod(model.grid.MX) * f_sigma;


% Division by k^2 to get the spectrum of the streamfunction
% And from omnidirectional spectrum to Fourier tranform square modulus
% Division by k in dimension 2
f_sigma = bsxfun(@times, 1 ./ ( kappa(2:end).^3) , f_sigma );

% % Division by k^2 to get the spectrum of the streamfunction
% f_sigma = bsxfun(@times, 1 ./ ( kappa(2:end) .^2) , f_sigma );
% 
% % From omnidirectional spectrum to Fourier tranform square modulus
% % Division by k in dimension 2
% f_sigma = f_sigma ./ ( kappa(2:end) );

% % Influence of discretisation
% f_sigma = f_sigma * d_kappa ;
% f_sigma = f_sigma / ( prod(MX.*dX) /(2*pi) ) ;

% From square modulus to modulus
f_sigma = sqrt( f_sigma );

% From 1D function to 2D function
f_sigma = interp1(kappa,[0; f_sigma],k);

% % Cleaning
% if strcmp(model.sigma.type_spectrum,'Band_Pass_w_Slope')
%     f_sigma(k<=k0)=0;
% end
f_sigma(k>k_inf)=0;
f_sigma=reshape(f_sigma,MX);

% Antialiasing
f_sigma(PX(1)+1,:)=0;
f_sigma(:,PX(2)+1)=0;

% % Influence of the complex brownian variance
% f_sigma = 1/sqrt(prod(MX))*f_sigma;

% % Choice to make the variance tensor explicitely independent of the
% % resolution (i.e. independent of MX and of dX, but dependent on dX.*MX)
% f_sigma = sqrt(prod(MX))*f_sigma;

% Orthogonal gradient (from streamfunction to velocity)
sigma_on_sq_dt(:,:,1)= 1i * ( - ky ) .* f_sigma;
sigma_on_sq_dt(:,:,2)= 1i * ( + kx ) .* f_sigma;


% Compute the bi-directional spectrum of sigma dBt
ft_sigma=abs(sigma_on_sq_dt).^2;
ft_sigma=sum(ft_sigma,3);

% % % Calcul of energy
% % % One has to divid by prod(model.grid.MX) because of the form of Parseval
% % % theorem for discrete Fourier transform
% % %trace_a_on_dt = sum(spectrum_a_sigma)
% % % trace_a_on_dt = 1/prod(model.grid.MX) * sum(spectrum_a_sigma);
% % % or equivalently
% % trace_a_on_dt = 1/prod(model.grid.MX) * sum(ft_sigma(:))
% trace_a = 1/prod(model.grid.MX) * sum(ft_sigma(:))

% Calcul of energy
% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
% One has to multiply by prod(model.grid.MX) because of the variance of
% the white in space noise
trace_a = 1/prod(model.grid.MX) * sum(ft_sigma(:));

if bool_plot
    
    % alpha = ( 3 - model.sigma.slope_absDif_sigma )/2;
    % eval(['norm_tr_a_theo = fct_norm_tr_a_theo_' ...
    %     model.sigma.type_spectrum '(model,k0,k_inf,alpha);']);
    % %norm_tr_a_theo = fct_norm_tr_a_theo(model,k0,k_inf,alpha);
    %
    % norm_tr_a_theo/trace_a_on_dt
    % % d_kappa
    % % prod(model.grid.MX)
    % % sqrt(prod(model.grid.MX))
    % % 1/(sqrt(prod(model.grid.MX))*d_kappa)
    % % d_kappa*prod(model.grid.MX)
    
    % Compute the omnidirectional spectrum of sigma dBt
    spectrum_a_sigma = idx' * ft_sigma(:);
    
    % Influence of the complex brownian variance
    spectrum_a_sigma = prod(MX)*spectrum_a_sigma;
    
    % % Choice to make the variance tensor explicitely independent of the
    % % resolution (i.e. independent of MX and of dX, but dependent on dX.*MX)
    % spectrum_a_sigma = (1/prod(MX))*spectrum_a_sigma;
    
    % Division by prod(model.grid.MX) because of the Parseval theorem for
    % discrete Fourier transform
    % Division by prod(model.grid.MX) again in order to the integration
    % of the spectrum over the wave number yields the energy of the
    % buoyancy averaged (not just integrated) over the space
    spectrum_a_sigma = 1/prod(model.grid.MX)^2 * spectrum_a_sigma;
    
    % % % Division by the wave number step
    % % spectrum_a_sigma = spectrum_a_sigma / d_kappa;
    %%
    
    % % % Division by prod(model.grid.MX) because of the Parseval theorem for
    % % % discrete Fourier transform
    % % Multiplication by prod(model.grid.MX) to normalize by the white-in-space noise
    % % variance
    % spectrum_a_sigma = prod(model.grid.MX) * spectrum_a_sigma;
    % %spectrum_a_sigma = 1/prod(model.grid.MX) * spectrum_a_sigma;
    
    % Division by the wave number step
    spectrum_a_sigma = spectrum_a_sigma / d_kappa;
    
    % % To remove the 2 pi which appear when we integrate k^(3-2-alpha) over the
    % % wave-vector angles and add the (2 pi)^2 which appear when we go from k to
    % % 2*pi*k
    % spectrum_a_sigma = spectrum_a_sigma * (2*pi);
    
    % Absolute diffusivity of the small-scale velocity (for test)
    trace_a_from_spectrum = d_kappa * sum(spectrum_a_sigma(:));
    
    
    %%
    
    % warning('debug here');
    % k2=(kx.^2+ky.^2);
    % reference_spectrum2 = k2 .^(1-alpha);
    % %reference_spectrum2 = k .^(2-2*alpha);
    % reference_spectrum2 = idx' * reference_spectrum2(:);
    % reference_spectrum2 = reference_spectrum2(2:end)';
    %
    % % Discretisation to go from continuous to discrete Fourier transform
    % reference_spectrum2 = 1/prod(model.grid.dX) * reference_spectrum2;
    %
    % % Division by prod(model.grid.MX) because of the Parseval theorem for
    % % discrete Fourier transform
    % reference_spectrum2 = 1/prod(model.grid.MX) * reference_spectrum2;
    %
    % % Division by the wave number step
    % reference_spectrum2 = reference_spectrum2 / d_kappa;
    %
    % % To remove the 2 pi which appear when we integrate k^(3-2-alpha) over the
    % % wave-vector angles
    % reference_spectrum2 = reference_spectrum2 * (2*pi);
    %
    % reference_spectrum(1)/reference_spectrum2(1)
    % reference_spectrum(1)/reference_spectrum2(1)*reference_spectrum2 ./ reference_spectrum
    %
    % figure(10);plot(kappa(2:end),reference_spectrum);hold on
    % plot(kappa(2:end),reference_spectrum2);
    
    %% Plot spectrum
    
    % % Make the plots appear at the same level thant the large-scale velocity spectrum
    % spectrum_a_sigma_plot = spectrum_a_sigma * mult_offset;
    spectrum_a_sigma_plot = spectrum_a_sigma;
    % reference_spectrum = reference_spectrum * mult_offset;
    
    taille_police = 12;
    X0 = [10 20];
    
    figure10=figure(10);
    widthtemp = 12;
    heighttemp = 6;
    set(figure10,'Units','inches', ...
        'Position',[X0 widthtemp heighttemp], ...
        'PaperPositionMode','auto');
    loglog(kappa(2:end),reference_spectrum_a,'k')
    hold on;
    loglog(kappa(2:end),reference_spectrum_a_estim,'k--')
    loglog(kappa(2:end),spectrum_w_a(2:end))
    %     if nargin>2
    %         loglog(kappa(2:end),spectrum_w2_a(2:end),'c')
    %     end
    loglog(kappa(2:end),spectrum_a_sigma_plot(2:end),'r.-')
    %     loglog(kappa(2:end),spectrum_w_a(2:end) + ...
    %         spectrum_a_sigma_plot(2:end),'g')
    ax=axis;
    ax(4)=max([spectrum_w_a; reference_spectrum_a_estim ; ...
        reference_spectrum_a ; spectrum_a_sigma_plot]);
    ax(4)=ax(4)*2;
    % ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_w_a); ...
    %     max(spectrum_a_sigma_plot)]);
    ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_w_a); ...
        max(reference_spectrum_a_estim); max(reference_spectrum_a); ...
        max(spectrum_a_sigma_plot)]);
    ax(3) = min( [ax(3) min(spectrum_a_sigma_plot) ...
        min(reference_spectrum_a) min(reference_spectrum_a_estim)]);
    ax(1:2)=kappa([2 end]);
    if ax(4)>0
        axis(ax)
    end
    ax = axis;
    loglog(km_LS*[1 1],...
        [min(reference_spectrum_a_estim) ax(4)],'k--')
    loglog(model.sigma.kappamin_on_kappamax * kappa(end)*[1 1],...
        [min(reference_spectrum_a_estim) ax(4)],'k-.')
    hold off
    set(gca,'XGrid','on','XTickMode','manual');
    width = 9;
    height = 3;
    set(figure10,'Units','inches', ...
        'Position',[X0 width height], ...
        'PaperPositionMode','auto');
    set(gca,'YGrid','on')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'FontName','Times')
    ylabel('$|\hat{f}(\kappa)|^2 \tau_\kappa$',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',taille_police,...
        'FontName','Times')
    xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',taille_police,...
        'interpreter','latex',...
        'FontName','Times')
    title('Absolute diffusivity by scale for $w$ and $\sigma dB_t$',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    %% Save plot
    drawnow
    
    folder_simu = model.folder.folder_simu;
    eval( ['print -depsc ' folder_simu '/AbsDiffByScale_sigma_dB_t/'...
        day '.eps']);
    
    % slope_w_a_comp_for_estim
    % km_LS
    % abs_diff_w
    % trace_a_on_2 = trace_a/2
end
end

function t = fct_unity_approx_(N_t)
% Approximation of unity
%

sslop=8;
t=ones(1,N_t);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

end