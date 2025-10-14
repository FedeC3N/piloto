function mris_r = spm12_reslice ( mri1, mris, opts )
% Reslices the image to fit a new space.
%
% im_res = SPM_RESLICE ( im_ref, im_ori, options )
%
% Where:
% * im_ref is the reference (static) image.
% * im_ori is the image to reslice.
% * options is an structure with fields:
%   - degree - the B-spline degree [default: 1]
%              Non-finite values result in Fourier interpolation. Note
%              that Fourier interpolation only works for purely rigid
%              body transformations. Voxel sizes must all be identical
%              and isotropic.
%
% * im_res is the image im_ori resliced to fit the space of im_ref.

% Based on SPM12 functions:
% * spm_reslice by John Ashburner

% Rigid body reslicing of images
% FORMAT spm_reslice(P,flags)
%
% P      - matrix or cell array of filenames {one string per row}
%          All operations are performed relative to the first image.
%          ie. Coregistration is to the first image, and resampling
%          of images is into the space of the first image.
%
% flags  - a structure containing various options.  The fields are:
%
%         interp - the B-spline interpolation method [default: 1]
%                  Non-finite values result in Fourier interpolation. Note
%                  that Fourier interpolation only works for purely rigid
%                  body transformations. Voxel sizes must all be identical
%                  and isotropic.
%
%__________________________________________________________________________
%
% The spatially realigned images are written to the original subdirectory
% with the same (prefixed) filename. They are all aligned with the first.
%
% Inputs:
% A series of images conforming to SPM data format (see 'Data Format'). The
% relative displacement of the images is stored in their header.
%
% Outputs:
% The routine uses information in their headers and writes the realigned 
% image files to the same subdirectory with a prefix.
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_reslice.m 5929 2014-03-27 14:47:40Z guillaume $

%__________________________________________________________________________
%
% The headers of the images contain a 4x4 affine transformation matrix 'M',
% usually affected by the `realignment' and `coregistration' modules.
% What these matrices contain is a mapping from the voxel coordinates
% (x0,y0,z0) (where the first voxel is at coordinate (1,1,1)), to 
% coordinates in millimeters (x1,y1,z1).
%
% x1 = M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4)
% y1 = M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4)
% z1 = M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4)
%
% Assuming that image1 has a transformation matrix M1, and image2 has a
% transformation matrix M2, the mapping from image1 to image2 is: M2\M1
% (ie. from the coordinate system of image1 into millimeters, followed
% by a mapping from millimeters into the space of image2).
%
% Several spatial transformations (realignment, coregistration,
% normalisation) can be combined into a single operation (without the
% necessity of resampling the images several times).
%__________________________________________________________________________
%
% Refs:
%
% Friston KJ, Williams SR, Howard R Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 35:346-355
%
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996) Improved Image
% Registration by Using Fourier Interpolation. Mag. Res. Med. 36(6):923-931
%
% R. W. Cox and A. Jesmanowicz (1999)  Real-Time 3D Image Registration
% for Functional MRI. Mag. Res. Med. 42(6):1014-1018
%__________________________________________________________________________

if nargin < 3
    opts = struct;
end

% Sets the defaults.
opts.degree = ft_getopt ( opts, 'degree', 4 );


% Initializes the realigned MRIs to the original ones.
mris_r = mris;

% Fourier method can only been used for rigid transformations.
if ~isfinite ( opts.degree )
    for i= 1: numel ( mris )
        trans = mri1.transform \ mris ( i ).transform;
        
        if any ( abs ( svd ( trans ( 1: 3, 1: 3 ) ) - 1 ) > 1e-7 )
            fprintf('\n  Zooms  or shears  appear to  be needed');
            fprintf('\n  (probably due to non-isotropic voxels).');
            fprintf('\n  These  can not yet be  done  using  the');
            fprintf('\n  Fourier reslicing method.  Switching to');
            fprintf('\n  7th degree B-spline interpolation instead.\n\n');
            opts.degree = 7;
            break
        end
    end
end


% Defines the degree of the B-spline function.
if isfinite ( opts.degree )
    degree = opts.degree;
else
    degree = 0;
end

% Creates a set of points to evaluate the B-spline function in.
x1 = 1: mri1.dim (1); x1 = permute ( x1 (:), [ 1 2 3 ] ); x1 = repmat ( x1, cat ( 2, 1, mri1.dim (2), mri1.dim (3) ) );
x2 = 1: mri1.dim (2); x2 = permute ( x2 (:), [ 3 1 2 ] ); x2 = repmat ( x2, cat ( 2, mri1.dim (1), 1, mri1.dim (3) ) );
x3 = 1: mri1.dim (3); x3 = permute ( x3 (:), [ 2 3 1 ] ); x3 = repmat ( x3, cat ( 2, mri1.dim (1), mri1.dim (2), 1 ) );

% Reslices each volume.
for i = 1: numel ( mris )
    
    % Calculates the desired tranformation for the current volume.
    trans   = mri1.transform \ mris (i).transform;
    
    % Gets the volume and transforms it to double if needed.
    voltype = class  ( mris (i).anatomy );
    vol     = double ( mris (i).anatomy );
    
    % Gets the value of air.
    air     = nanmin ( vol (:) );
    
    % Generates the B-spline function.
    spline  = spm12_bsplinc ( vol, degree );
    
    if ~isfinite ( opts.degree )
        
        % Thansforms the spline function to the new space.
        vol     = abs ( kspace3d ( spline, trans ) );
    else
        
        % Evaluates the function at the points in the referece MRI.
        [ y1, y2, y3 ] = transform ( inv ( trans ), x1, x2, x3 );
        vol     = spm12_bsplins ( spline, y1, y2, y3, degree );
    end
    
    % Sets the NaN to real zeros (air).
    vol     = nanmax ( vol, air );
    
    % Returns the volume to its orginal type.
    vol     = cast ( vol, voltype );
    
    % Stores the new volume and its metadata.
    mris_r (i).dim       = mri1.dim;
    mris_r (i).anatomy   = vol;
    mris_r (i).transform = mri1.transform;
end


%==========================================================================
%-function v = kspace3d(v,M)
%==========================================================================
function v = kspace3d(v,M)
% 3D rigid body transformation performed as shears in 1D Fourier space
% FORMAT v = kspace3d(v,M)
% v        - image stored as a 3D array
% M        - rigid body transformation matrix
%
% v        - transformed image
%
% References:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI
% Magnetic Resonance in Medicine 42(6):1014-1018
%
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996)
% Improved Image Registration by Using Fourier Interpolation
% Magnetic Resonance in Medicine 36(6):923-931

[S0,S1,S2,S3] = shear_decomp(M);

d  = [size(v) 1 1 1];
g = 2.^ceil(log2(d));
if any(g~=d)
    tmp = v;
    v   = zeros(g);
    v(1:d(1),1:d(2),1:d(3)) = tmp;
    clear tmp;
end

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2)
    t        = reshape( exp((j*S3(3,2) + S3(3,1)*(1:g(1)) + S3(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
    v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end

% XZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(2)-1)/2) 0 (-g(2)/2+1):-1])/g(2);
for k=1:g(3)
    t        = exp( (k*S2(2,3) + S2(2,1)*(1:g(1)) + S2(2,4)).'*tmp1);
    v(:,:,k) = real(ifft(fft(v(:,:,k),[],2).*t,[],2));
end

% YZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(1)-1)/2) 0 (-g(1)/2+1):-1])/g(1);
for k=1:g(3)
    t        = exp( tmp1.'*(k*S1(1,3) + S1(1,2)*(1:g(2)) + S1(1,4)));
    v(:,:,k) = real(ifft(fft(v(:,:,k),[],1).*t,[],1));
end

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2)
    t        = reshape( exp( (j*S0(3,2) + S0(3,1)*(1:g(1)) + S0(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
    v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end

if any(g~=d), v = v(1:d(1),1:d(2),1:d(3)); end


%==========================================================================
%-function [S0,S1,S2,S3] = shear_decomp(A)
%==========================================================================
function [S0,S1,S2,S3] = shear_decomp(A)
% Decompose rotation and translation matrix A into shears S0, S1, S2 and
% S3, such that A = S0*S1*S2*S3. The original procedure is documented in:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI
% Magnetic Resonance in Medicine 42(6):1014-1018

A0 = A(1:3,1:3);
if any(abs(svd(A0)-1)>1e-7), error('Can''t decompose matrix'); end

t  = A0(2,3); if t==0, t=eps; end
a0 = pinv(A0([1 2],[2 3])')*[(A0(3,2)-(A0(2,2)-1)/t) (A0(3,3)-1)]';
S0 = [1 0 0; 0 1 0; a0(1) a0(2) 1];
A1 = S0\A0;  a1 = pinv(A1([2 3],[2 3])')*A1(1,[2 3])';  S1 = [1 a1(1) a1(2); 0 1 0; 0 0 1];
A2 = S1\A1;  a2 = pinv(A2([1 3],[1 3])')*A2(2,[1 3])';  S2 = [1 0 0; a2(1) 1 a2(2); 0 0 1];
A3 = S2\A2;  a3 = pinv(A3([1 2],[1 2])')*A3(3,[1 2])';  S3 = [1 0 0; 0 1 0; a3(1) a3(2) 1];

s3 = A(3,4)-a0(1)*A(1,4)-a0(2)*A(2,4);
s1 = A(1,4)-a1(1)*A(2,4);
s2 = A(2,4);
S0 = [[S0 [0  0 s3]'];[0 0 0 1]];
S1 = [[S1 [s1 0  0]'];[0 0 0 1]];
S2 = [[S2 [0 s2  0]'];[0 0 0 1]];
S3 = [[S3 [0  0  0]'];[0 0 0 1]];


%==========================================================================
%-function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
%==========================================================================
function [ y1, y2, y3 ] = transform ( M, x1, x2, x3 )
y1   = M ( 1, 1 ) * x1 + M ( 1, 2 ) * x2 + ( M ( 1, 3 ) * x3 + M ( 1, 4 ) );
y2   = M ( 2, 1 ) * x1 + M ( 2, 2 ) * x2 + ( M ( 2, 3 ) * x3 + M ( 2, 4 ) );
y3   = M ( 3, 1 ) * x1 + M ( 3, 2 ) * x2 + ( M ( 3, 3 ) * x3 + M ( 3, 4 ) );
