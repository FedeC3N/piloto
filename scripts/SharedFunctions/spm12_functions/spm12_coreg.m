function transform = spm12_coreg ( mri1, mris, opts )
% Coregistrates images using information theory.
%
% matrix = SPM_COREG ( im_ref, im_dis, options )
%
% Where:
% * im_ref is the reference (static) image.
% * im_dis is the displaced image.
% * options is an structure with fields:
%          sep      - optimisation sampling steps (mm)
%                     default: [4 2]
%          params   - starting estimates (6 elements)
%                     default: [0 0 0  0 0 0]
%          func     - cost function string:
%                       'mi'  - Mutual Information
%                       'nmi' - Normalised Mutual Information
%                       'ecc' - Entropy Correlation Coefficient
%                       'ncc' - Normalised Cross Correlation
%                     default: 'nmi'
%          tol      - tolerences for accuracy of each param
%                     default: [0.02 0.02 0.02 0.001 0.001 0.001]
%          fwhm     - smoothing to apply to 256x256 joint histogram
%                     default: [7 7]
%
% * matrix is the transformation matrix to coregistrate im_dis to im_ref.

% Based on SPM12 functions:
% * spm_coreg by John Ashburner

% The registration method used here is based on the work described in:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
%
% The original interpolation method described in this paper has been
% changed in order to give a smoother cost function.  The images are
% also smoothed slightly, as is the histogram.  This is all in order to
% make the cost function as smooth as possible, to give faster convergence
% and less chance of local minima.
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_coreg.m 6435 2015-05-14 09:59:54Z guillaume $

%--------------------------------------------------------------------------
% References
%==========================================================================
%
% Mutual Information
% -------------------------------------------------------------------------
% Collignon, Maes, Delaere, Vandermeulen, Suetens & Marchal (1995).
% "Automated multi-modality image registration based on information theory".
% In Bizais, Barillot & Di Paola, editors, Proc. Information Processing
% in Medical Imaging, pages 263--274, Dordrecht, The Netherlands, 1995.
% Kluwer Academic Publishers.
%
% Wells III, Viola, Atsumi, Nakajima & Kikinis (1996).
% "Multi-modal volume registration by maximisation of mutual information".
% Medical Image Analysis, 1(1):35-51, 1996. 
%
% Entropy Correlation Coefficient
% -------------------------------------------------------------------------
% Maes, Collignon, Vandermeulen, Marchal & Suetens (1997).
% "Multimodality image registration by maximisation of mutual
% information". IEEE Transactions on Medical Imaging 16(2):187-198
%
% Normalised Mutual Information
% -------------------------------------------------------------------------
% Studholme, Hill & Hawkes (1998).
% "A normalized entropy measure of 3-D medical image alignment".
% in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.             
%
% Optimisation
% -------------------------------------------------------------------------
% Press, Teukolsky, Vetterling & Flannery (1992).
% "Numerical Recipes in C (Second Edition)".
% Published by Cambridge.
%--------------------------------------------------------------------------

if nargin < 3
    opts = struct;
end

% Sets the defaults.
opts.func   = ft_getopt ( opts, 'func',   'nmi' );
opts.sep    = ft_getopt ( opts, 'sep',    [ 4 2 ] );
opts.tol    = ft_getopt ( opts, 'tol',    [ 20 20 20  1  1  1 10 10 10  1  1  1 ] / 1e3 );
opts.fwhm   = ft_getopt ( opts, 'fwhm',   [ 7 7 ] );
opts.params = ft_getopt ( opts, 'params', [ 0 0 0 0 0 0 ] );


% Gets the gaussian kernel FWHM.
vsize       = sqrt(sum(mri1.transform(1:3,1:3).^2));
fwhm        = sqrt(max([1 1 1]*opts.sep(end)^2 - vsize.^2, [0 0 0]))./vsize;

% Expands the MRI to 0-255 and smooths the result.
mri1.uint8  = vol2uint8 ( mri1 );
mri1.uint8  = smoothuint8 ( mri1.uint8, fwhm );


sc = opts.tol (:)'; % Required accuracy
sc = sc ( 1: length ( opts.params ) );
xi = diag ( sc * 20 );

% Reserves memory for the output.
transform   = zeros ( 4,4, numel ( mris ) );

for vindex = 1: numel ( mris )
    
    % Gets the current MRI.
    mri2        = mris ( vindex );
    
    % Gets the gaussian kernel FWHM.
    vsize       = sqrt ( sum ( mri2.transform ( 1: 3, 1: 3 ) .^ 2 ) );
    fwhm        = sqrt ( opts.sep ( end ) ^ 2 - vsize .^ 2 ) ./ vsize;
    fwhm        = real ( fwhm );
    
    % Expands the MRI to 0-255 and smooths the result.
    mri2.uint8  = vol2uint8 ( mri2 );
    mri2.uint8  = smoothuint8 ( mri2.uint8, fwhm );
    
    vector      = opts.params(:);
    for samp=opts.sep(:)'
        vector      = spm12_powell ( vector (:), xi, sc, 'spm12_coreg_optfun', mri1, mri2, samp, opts.func, opts.fwhm );
    end
    transform ( :, :, vindex ) = spm_matrix ( vector (:)' );
end


function vol = vol2uint8 ( mri )

% Gets only the finite values.
vol = double ( mri.anatomy );

% Gets the maximum and minimum values of the image.
mx  = nanmax ( vol (:) );
mn  = nanmin ( vol (:) );

% Redefines the maximum to include 99.99% of the voxels.
nh  = 2048;
x   = round ( ( vol ( isfinite ( vol (:) ) ) + ((mx-mn)/(nh-1)-mn) ) * ((nh-1)/(mx-mn)));
h   = accumarray ( x (:), 1, [ nh 1 ] );
tmp = [find(cumsum(h)/sum(h)>0.9999); nh];
mx  = (mn*nh-mx+tmp(1)*(mx-mn))/(nh-1);

% Expands the image to the range 0-255.
vol = ( vol - mn ) * 255 / ( mx - mn );
vol = max ( vol, 0 );
vol = min ( vol, 255 );
vol = uint8 ( vol );


function vol = smoothuint8 ( vol, fwhm )
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol ( vol, vol, x, y, z, - [ i j k ] );
