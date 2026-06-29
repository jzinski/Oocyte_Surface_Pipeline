function h = fspecial3_mod(type, hsize, sigma)
%FSPECIAL3_MOD  Create a 3-D filter kernel (Gaussian).
%
%   h = fspecial3_mod('gaussian', hsize, sigma) returns a 3-D Gaussian
%   lowpass filter kernel.
%
%   Inputs
%     type  : filter type. Only 'gaussian' is implemented (the only type
%             used by this pipeline, via smooth3D.m).
%     hsize : kernel size in voxels. Scalar, or a 1x3 vector [nx ny nz].
%     sigma : Gaussian standard deviation. Scalar, or a 1x3 vector
%             [sx sy sz]. Defaults to hsize if omitted.
%
%   Output
%     h     : 3-D array (size hsize) holding a separable Gaussian kernel
%             normalized to sum to 1, suitable for use with imfilter.
%
%   Notes
%     This is a small, self-contained helper so the pipeline runs without
%     any function on a user's private MATLAB path. It reproduces the
%     standard normalized separable 3-D Gaussian used by MATLAB's fspecial
%     family (centered grid via (hsize-1)/2; values below eps*max set to 0;
%     normalized to unit sum). It does not require any toolbox.
%
%   Example
%     h = fspecial3_mod('gaussian', [24 24 3], [24 24 3]);

if nargin < 3 || isempty(sigma)
    sigma = hsize;
end

if isscalar(hsize), hsize = [hsize hsize hsize]; end
if isscalar(sigma), sigma = [sigma sigma sigma]; end

hsize = double(hsize(:)).';   % force 1x3 row
sigma = double(sigma(:)).';

if numel(hsize) ~= 3 || numel(sigma) ~= 3
    error('fspecial3_mod:badSize', ...
          'hsize and sigma must each be a scalar or a 1x3 vector.');
end
if any(sigma <= 0)
    error('fspecial3_mod:badSigma', 'sigma values must be positive.');
end

switch lower(type)
    case 'gaussian'
        siz = (hsize - 1) / 2;
        [x, y, z] = ndgrid(-siz(1):siz(1), -siz(2):siz(2), -siz(3):siz(3));
        h = exp( -(x.^2 / (2*sigma(1)^2) ...
                 + y.^2 / (2*sigma(2)^2) ...
                 + z.^2 / (2*sigma(3)^2)) );
        h(h < eps * max(h(:))) = 0;
        s = sum(h(:));
        if s ~= 0
            h = h / s;
        end
    otherwise
        error('fspecial3_mod:unsupportedType', ...
              'Only the ''gaussian'' filter type is implemented.');
end
end
