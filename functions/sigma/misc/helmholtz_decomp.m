function [w1, w2, w3] = helmholtz_decomp(model, w, varargin)
%% function [...] = helmholtz_decomp(model, w [, which])
% Computes the Helmholtz decomposition of velocity field w.
%
% Arguments:
%   * model: the sqgmu model, the grid in Fourier (model.grid.k) having been
%       initialized;
%   * w: the velocity field;
%   * [optional] decomp: the required decomposition, either 'all'
%       (default), 'sol' (solenoidal only) or 'irr' (irrotational only).
%
% Output: depends on 'decomp'.
%   * with 'all' (default), [w_sol, w_irr, w_harm]: the solenoidal, 
%       irrotational and harmonic components;
%   * with 'sol', w_sol the solenoidal component only;
%   * with 'irr', w_irr the irrotational component only.
%
% Written by P. DERIAN 2016-08-31.

% Parse input
decomp = parse_inputs(varargin);
% Fourier transform of the velocity
f_u = fft2(w(:,:,1));
f_v = fft2(w(:,:,2));
% Apply decomposition
switch decomp
    case 'all'
        w1 = solenoidal(model, f_u, f_v);
        w2 = irrotational(model, f_u, f_v);
        w3 = w - (w1 + w2);
    case 'sol'
        w1 = solenoidal(model, f_u, f_v);
    case 'irr'
        w1 = irrotational(model, f_u, f_v);
    otherwise
        error('SQGMU:helmholtz_decomp:InvalidInputs', 'decomposition should be either one of "all", "irr", "sol"');
end
end

function wir = irrotational(model, f_u, f_v)
%% function wir = irrotational(model, f_u, f_v)
% Computes the irrotational part.

% Fourier transform of the irrotational velocity
tmp = (model.grid.k.kx.*f_u + model.grid.k.ky.*f_v);
f_uir = model.grid.k.kx_over_ksqr.*tmp;
f_vir = model.grid.k.ky_over_ksqr.*tmp;

% Inverse transform
wir = zeros(size(f_u));
wir(:,:,1) = real(ifft2(f_uir));
wir(:,:,2) = real(ifft2(f_vir));
end

function wsol = solenoidal(model, f_u, f_v)
    % Fourier transform of the solenoidal velocity
    tmp = (-model.grid.k.ky.*f_u + model.grid.k.kx.*f_v);
    f_usol = -(model.grid.k.ky_over_ksqr.*tmp);
    f_vsol = model.grid.k.kx_over_ksqr.*tmp;

    % Inverse transform
    wsol = zeros(size(f_u));
    wsol(:,:,1) = real(ifft2(f_usol));
    wsol(:,:,2) = real(ifft2(f_vsol));
end

function decomp = parse_inputs(v)
%% parse optional arguments
if isempty(v)
    decomp = 'all';
else
    decomp = v{1};
end
end
