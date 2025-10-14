function output = spm12_preproc_write ( res, wtc, wbf, wdf, mrf, cleanup, pbb, pvx )
% Write out VBM preprocessed data
% FORMAT [cls,M1] = spm12_preproc_write8(res,tc,bf,df,mrf,cleanup,bb,vx)
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm12_preproc_write8.m 6137 2014-08-19 12:43:11Z john $

% Prior adjustment factor.
% This is a fudge factor to weaken the effects of the tissue priors.  The
% idea here is that the bias from the tissue priors probably needs to be
% reduced because of the spatial smoothing typically used in VBM studies.
% Having the optimal bias/variance tradeoff for each voxel is not the same
% as having the optimal tradeoff for weighted averages over several voxels.


mri = res.mri;

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end


% Initializes the output.
output    = [];

% Defines the output files base name.
basename  = 'spm12mri';
getwarp   = [ 0 1 ];


% Extracts data and metadata from the MRI.
mridim    = mri.dim ( 1: 3 );
mritrans  = mri.transform;
vsize     = sqrt ( sum ( mritrans ( 1: 3, 1: 3 ) .^2 ) );
volumes   = mri.anatomy;
volumes   = single ( volumes );
numvols   = size ( volumes, 4 );

clear mri;


% Initializes the undefined parameters.
if nargin < 2, wtc = false ( Kb, 4 );      end % native, import, warped, warped-mod
if nargin < 3, wbf = false ( numvols, 2 ); end % field, corrected
if nargin < 4, wdf = false ( 1, 2 );       end % inverse, forward
if nargin < 5, mrf = 1;                    end % MRF parameter
if nargin < 6, cleanup = 1;                end % Run the ad hoc cleanup
if nargin < 7, pbb = NaN ( 2, 3 );         end % Default to TPM bounding box
if nargin < 8, pvx = NaN;                  end % Default to TPM voxel size

% Gets the TPM template dimensions and its transformation matrix.
tpm       = res.tpm;
tpmdim    = size ( tpm.dat {1} );
tpmdim    = tpmdim ( 1: 3 );
tpmtrans  = tpm.M;

afftrans  = res.Affine;

% Extracts the TPM template bounding box and voxel size.
bb        = spm12_get_bbox ( tpm.mri, 'old' );
vx        = sqrt ( sum ( tpmtrans ( 1: 3, 1: 3 ) .^ 2, 1 ) );
% vx        = geomean ( vx );

% Replaces the calculated values with the provided, if any.
bb ( isfinite ( pbb ) ) = pbb ( isfinite ( pbb ) );
vx ( isfinite ( pvx ) ) = pvx ( isfinite ( pvx ) );

% Rounds the bounding box to the nearest voxel.
bb        = bsxfun ( @times, vx, round ( bsxfun ( @rdivide, bb, vx ) ) );

% Gets the bounding box dimensions in voxels.
bbdim     = abs ( round ( diff ( bb ) ./ vx ) ) + 1;


% Calculates TPM if some tissue map or output parameters are requested.
% Calculates the deformation fields if they or the TPM are requested.
do_cls    = any ( wtc (:) ) || nargout >= 1;
do_defs   = any ( wdf ) || do_cls;



% Calculates the deformation fields.
if do_defs
    
    % Gets the warp deformation matrix.
    Swarp = res.Twarp;
    
    % Gets the transformation between MRI voxels and TPM voxels.
    M = tpmtrans \ afftrans * mritrans;
    
    % Defines the voxel coordiantes.
    xs ( :, 1, 1 ) = 1: mridim (1);
    ys ( 1, :, 1 ) = 1: mridim (2);
    zs ( 1, 1, : ) = 1: mridim (3);
    
    % Transforms the voxels coordinates to sampled space.
    downstrans = res.MT;
    upstrans   = inv ( downstrans );
    xss = xs * upstrans (1,1) + upstrans (1,4);
    yss = ys * upstrans (2,2) + upstrans (2,4);
    zss = zs * upstrans (3,3) + upstrans (3,4);
    [ xsr, ysr, zsr ] = ndgrid ( xss, yss, zss );
    
    % Expands the warp from the subsampled to the original spaces.
    prm     = [ 3 3 3 0 0 0 ];
    BSwarp ( :, :, :, 1 ) = spm_bsplinc ( Swarp  ( :, :, :, 1 ), prm );
    Twarp  ( :, :, :, 1 ) = spm_bsplins ( BSwarp ( :, :, :, 1 ), xsr, ysr, zsr, prm );
    BSwarp ( :, :, :, 2 ) = spm_bsplinc ( Swarp  ( :, :, :, 2 ), prm );
    Twarp  ( :, :, :, 2 ) = spm_bsplins ( BSwarp ( :, :, :, 2 ), xsr, ysr, zsr, prm );
    BSwarp ( :, :, :, 3 ) = spm_bsplinc ( Swarp  ( :, :, :, 3 ), prm );
    Twarp  ( :, :, :, 3 ) = spm_bsplins ( BSwarp ( :, :, :, 3 ), xsr, ysr, zsr, prm );
    
    % Warps and transforms the voxel coordinates to TPM-space.
    [ xsr, ysr, zsr ] = warpaffine ( xs, ys, zs, M, Twarp );
    
    clear Twarp
    
    % Saves the deformation fields.
    if wdf (1) || getwarp (1)
        [ NDx, NDy, NDz ] = warpaffine ( xsr, ysr, zsr, tpmtrans );
        Ndef = single ( cat ( 4, NDx, NDy, NDz ) );
        
        clear NDx NDy NDz
    end
    
    % Saves the TPM-world coordinates.
    if wdf (2) || any ( any ( wtc ( :, 2: 4 ) ) ) || nargout >= 1
        y = single ( cat ( 4, xsr, ysr, zsr ) );
    end
    
    % Calculates the TPM.
    if do_cls
        
        % Generates the bias field and the bias corrected field.
        bf = cell ( 1, numvols );
        cr = cell ( 1, numvols );
        
        % Goes through each volume.
        for vindex = 1: numvols
            
            % Calculates the DCT of the bias field and the corrected data.
            bf { vindex } = single ( exp ( dct3 ( res.Tbias { vindex }, mridim, xs, ys, zs ) ) );
            cr { vindex } = bf { vindex } .* volumes ( :, :, :, vindex );
        end
        
        clear volumes
        
        % Uses a parametric representation of intensity distributions.
        if isfield ( res, 'mg' )
            
            % Gets the a-priori probability of each voxel.
            b   = spm12_sample_priors8 ( tpm, xsr (:), ysr (:), zsr (:) );
            b   = cat ( 2, b {:} );
            b   = bsxfun ( @times, res.wp, b );
            b   = bsxfun ( @rdivide, b, sum ( b, 2 ) );
            
            clear tpm xsr ysr zsr
            
            % Gets the likelihood for each voxel.
            Qt  = likelihoods ( cr, [], res.mg, res.mn, res.vr );
            
            Q   = zeros ( size ( b ), 'single' );
            for k1 = 1: Kb
                Q ( :, k1 ) = sum ( Qt ( :, lkp == k1 ), 2 );
            end
            
            % Multiplies the a-priori and the likelihood.
            Q = Q .* b;
            
            clear Qt b
            
            % Reshapes the result in matrix form.
            Q = reshape ( Q, cat ( 2, mridim, Kb ) );
            
        % Uses a nonparametric representation of intensity distributions.
        else
            warning ( 'Not yet programmed.' );
%             q   = spm12_sample_priors8(tpm,t1,t2,t3);
%             wp  = res.wp;
%             s   = zeros(size(q{1}));
%             for k1 = 1:Kb,
%                 q{k1} = wp(k1)*q{k1};
%                 s     = s + q{k1};
%             end
%             for k1 = 1:Kb,
%                 q{k1} = q{k1}./s;
%             end
%             q   = cat(3,q{:});
%             
%             for n = 1: numvols
%                 tmp = round(cr{n}*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
%                 tmp = min(max(tmp,1),size(res.intensity(n).lik,1));
%                 for k1=1:Kb,
%                     likelihood = res.intensity(n).lik(:,k1);
%                     q(:,:,k1)  = q(:,:,k1).*likelihood(tmp);
%                 end
%             end
        end
    end
end

% Sanitizes the tissue classes.
if do_cls
    
    % If no MRF level normalizes the tissues so its sum is 1.
    if mrf == 0
        P = bsxfun ( @rdivide, Q, sum ( Q, 4 ) + eps );
        P = uint8 ( round ( P * 255 ) );
        
    % Otherwise uses an iterative MRF procedure.
    else
        
        % Sets the parameters and initializes the output.
        P   = zeros ( size ( Q ), 'uint8' );
        G   = ones ( Kb, 1, 'single' ) * mrf;
        
        % Uses 10 iterations.
        for iter = 1: 10
            spm12_mrf ( P, Q, G, single ( vsize .^ 2 ) );
        end
    end
end

clear Q

% Uses an ad-hoc cleanup procedure, if requested.
if cleanup
    if size ( P, 4 ) > 5
        
        % Cleanup is performed in 8-bit unsigned integer data.
        P = clean_gwc ( P, cleanup );
    else
        warning ( 'Cleanup not done. Not enough tissues' );
    end
end

% Rewrites the tissue classes as single precision data.
P     = single ( P ) / 255;

% Puts the tissue classes into a cell array.
tc1   = num2cell ( P, [ 1 2 3 ] );
tc1   = squeeze ( tc1 );

clear P

% Creates the rest of tissue class images.
tc2   = cell ( size ( tc1 ) );
tc3   = cell ( size ( tc1 ) );
tc4   = cell ( size ( tc1 ) );

% Creates the old output variable.
cls   = cell ( size ( tc1 ) );


% Computes the voxel to 'imported'-world  transformation.
mm  = cat ( 1, ...
    cat ( 2, bb (1,1), bb (1,2), bb (1,3), 1 ), ...
    cat ( 2, bb (2,1), bb (1,2), bb (1,3), 1 ), ...
    cat ( 2, bb (1,1), bb (2,2), bb (1,3), 1 ), ...
    cat ( 2, bb (2,1), bb (2,2), bb (1,3), 1 ), ...
    cat ( 2, bb (1,1), bb (1,2), bb (2,3), 1 ), ...
    cat ( 2, bb (2,1), bb (1,2), bb (2,3), 1 ), ...
    cat ( 2, bb (1,1), bb (2,2), bb (2,3), 1 ), ...
    cat ( 2, bb (2,1), bb (2,2), bb (2,3), 1 ) )';
vx3 = cat ( 1, ...
    cat ( 2,         1,         1,         1,         1 ), ...
    cat ( 2, bbdim (1),         1,         1,         1 ), ...
    cat ( 2,         1, bbdim (2),         1,         1 ), ...
    cat ( 2, bbdim (1), bbdim (2),         1,         1 ), ...
    cat ( 2,         1,         1, bbdim (3),         1 ), ...
    cat ( 2, bbdim (1),         1, bbdim (3),         1 ), ...
    cat ( 2,         1, bbdim (2), bbdim (3),         1 ), ...
    cat ( 2, bbdim (1), bbdim (2), bbdim (3),         1 ) )';
imptrans = mm / vx3;

% Calculates the ''imported'' images, if requested.
if any ( any ( wtc ( :, 2: 4 ) ) ) || nargout >= 1 || wdf (2)
    
    % Generates the ''imported'' tissue class images.
    if any ( wtc ( :, 2 ) )
        warning ( 'Not yet programmed.' );
        wtc ( :, 2 ) = false;
        
%         % Calculates the real-world coordinates of each voxel.
%         xyzsr  = warpaffine ( xs, ys, zs, mritrans );
%         
%         % Calculates the TPM-world coordinates of each deformed voxel.
%         xyztpm = warpaffine ( y ( :, :, :, 1 ), y ( :, :, :, 2 ), y ( :, :, :, 3 ), tpmtrans );
%         
%         % Computes a rigid-body transformation from real world to TPM.
%         [ ~, R ] = spm12_get_closest_affine ( xyzsr, xyztpm, tc1 {1} );
%         
%         % Calculates the transformation TPM voxels and image's world.
%         mat0   = R \ imptrans;
%         
%         % Defines an anti-alias Gaussian filter.
%         fwhm   = max ( vx ./ vsize - 1, 0.01 );
%         
%         for k1 = 1: Kb
%             if wtc ( k1, 2 ),
%                 
%                 % Filters the image low-pass using the Gaussian filter.
%                 tc2 { k1 } = decimate ( tc1 { k1 }, fwhm );
%                 
%                 % Interpolates the data.
%                 for z = 1: bbdim (3),
%                     tc2 { k1 } ( :, :, z ) = spm_slice_vol ( tmp1, mritrans \ mat0 * spm_matrix ( [ 0 0 i ] ), bbdim ( 1: 2 ), [ 1, NaN ] ) / 255;
%                 end
%             end
%         end
    end
    
    % Adjust the transformations to get the desired bounding box.
    if max ( abs ( imptrans (:) - tpmtrans (:) ) ) > 1e-10
        
        % Transformation the warped coordinates to fit the bounding box.
        y = warpaffine ( y ( :, :, :, 1 ), y ( :, :, :, 2 ), y ( :, :, :, 3 ), imptrans \ tpmtrans );
    end
    
    % From now on uses the imported space.
    tpmtrans = imptrans;
    tpmdim   = bbdim;
    
    if any ( any ( wtc ( :, 3: 4 ) ) ) || nargout >= 1 || wdf (2)
        
        % Prepares the data for the warped tissue class.
        if any ( wtc ( :, 3 ) )
            C = zeros ( cat ( 2, tpmdim, Kb ), 'single' );
            
            % Sets the type of boundary.
            spm12_field ( 'boundary', 1 );
            
            % Goes through each tissue class.
            for k1 = 1: Kb
                
                % Corrects the original output variable.
                tc          = tc1 { k1 };
                [ tc, w ]   = spm12_diffeo ( 'push', tc, y, tpmdim ( 1: 3 ) );
                cls { k1 }  = tc;
                
                % Gets the voxel size.
                vx          = sqrt ( sum ( tpmtrans ( 1: 3, 1: 3 ) .^ 2 ) );
                
                % Regularizes the tissue class.
                C ( :, :, :, k1 ) = spm12_field ( w, tc, [ vx 1e-6 1e-4 0 3 2 ] );
            end
            
            % Calculates the correction parameter.
            correction = sum ( max ( C, eps ), 4 );
            
            % Goes through each tissue class.
            for k1= 1: Kb
                if wtc ( k1, 3 )
                    tc3 { k1 } = C ( :, :, :, k1 ) ./ correction;
                end
            end
            
        % If no warped tissue class only corrects the data.    
        else
            
            % Goes through each tissue class.
            for k1= 1: Kb
                
                % Corrects the original output variable.
                tc          = tc1 { k1 };
                tc          = spm12_diffeo ( 'push', tc, y, tpmdim ( 1: 3 ) );
                cls { k1 }  = tc;
            end
        end
        
        % Calculates the corrected warped tissue class.
        if any ( wtc ( :, 4 ) )
            
            % Calculates the correction parameter.
            correction = abs ( det ( mritrans ( 1: 3, 1: 3 ) ) / det ( tpmtrans ( 1: 3, 1: 3 ) ) );
            
            % Goes through each tissue class.
            for k1 = 1: Kb
                if wtc ( k1, 4 )
                    tc4 { k1 } = cls { k1 } * correction;
                end
            end
        end
    end
end

% Calculates the inverse deformation field.
if wdf (2) || getwarp (2)
    
    % Calculates the inverse deformation field for each dimension.
%     y = spm12_diffeo ( 'invdef', y, tpmdim, eye (4), mritrans );
    y1 = y ( :, :, :, 1 );
    y2 = y ( :, :, :, 2 );
    y3 = y ( :, :, :, 3 );
    [ y1, y2, y3 ] = spm_invdef ( y1, y2, y3, tpmdim, eye (4), mritrans );
    
    % Joins the three dimensions.
    y = cat ( 4, y1, y2, y3 );
    
    % Extrapolate the points with NaNs.
    y = spm12_extrapolate_def ( y, tpmtrans );
end

% Sets the output.
output.mri = cat ( 4, cr {:} );
output.tc1 = tc1;

% Includes the linear affine transformation and the non-linear warp.
if getwarp (1)
    output.sub2mni.trans = afftrans;
    output.sub2mni.warp  = Ndef;
end
if getwarp (2)
    output.mni2sub.trans = inv ( afftrans );
    output.mni2sub.warp  = y;
end

% Saves the requested MRI files.

% Writes out the tissue class images, if requested.
for tindex = 1: Kb
    
    % Writes out the original tissue class image.
    if wtc ( tindex, 1 )
        filename = sprintf ( 'tc%.0i_%s.nii', tindex, basename );
        ft_write_mri ( filename, tc1 { tindex }, 'transform', mritrans, 'dataformat', 'nifti' );
    end
    
    % Writes out the 'imported' tissue class image.
    if wtc ( tindex, 2 )
        filename = sprintf ( 'rc%.0i_%s.nii', tindex, basename );
        ft_write_mri ( filename, tc2 { tindex }, 'transform', imptrans, 'dataformat', 'nifti' );
    end
    
    % Writes out the warped tissue class image.
    if wtc ( tindex, 3 )
        filename = sprintf ( 'wc%.0i_%s.nii', tindex, basename );
        ft_write_mri ( filename, tc3 { tindex }, 'transform', tpmtrans, 'dataformat', 'nifti' );
    end
    
    % Writes out the modified warped tissue class image.
    if wtc ( tindex, 4 )
        filename = sprintf ( 'mwc%.0i_%s.nii', tindex, basename );
        ft_write_mri ( filename, tc4 { tindex }, 'transform', tpmtrans, 'dataformat', 'nifti' );
    end
end


% Writes out the bias fields and the corrected volumes, if requested.
for vindex = 1: numvols
    
    % Writes out the bias field, if requested.
    if wbf ( vindex, 1 )
        filename = sprintf ( 'bf_%.0i_%s.nii', vindex, basename );
        ft_write_mri ( filename, bf { vindex }, 'transform', mritrans, 'dataformat', 'nifti' );
    end
    
    % Writes out the bias field-corrected MRI, if requested.
    if wbf ( vindex, 2 )
        filename = sprintf ( 'm_%.0i_%s.nii', vindex, basename );
        ft_write_mri ( filename, cr { vindex }, 'transform', mritrans, 'dataformat', 'nifti' );
    end
end


% Writes out the inverse warp deformation field.
if wdf (1)
    filename = sprintf ( 'iy_%s.nii', basename );
    ft_write_mri ( filename, Ndef, 'transform', mritrans, 'dataformat', 'nifti' );
end

% Writes out the warp-deformated coordinates.
if wdf (2)
    filename = sprintf ( 'y_%s.nii', basename );
    ft_write_mri ( filename, y, 'transform', tpmtrans, 'dataformat', 'nifti' );
end


%==========================================================================
% function t = transf(B1,B2,B3,T)
%==========================================================================
% function t = transf(B1,B2,B3,T)
% 
% % Reserves memory for the DCT.
% t = zeros ( size ( B1, 1 ), size ( B2, 1 ), size ( B3, 1 ) );
% 
% if isempty ( T ), return, end
% 
% d2 = size ( T );
% t1 = reshape ( T, d2 (1) * d2 (2), d2 (3) );
% t1 = t1 * B3';
% t1 = reshape ( t1, d2 (1), d2 (2), [] );
% for z = 1: size ( B3, 1 )
%     t ( :, :, z ) = B1 * t1 ( :, :, z ) * B2';
% end

function t = dct3 ( T, order, x, y, z )

% Reserves memory for the DCT.
t = zeros ( order );

if isempty ( T ), return, end

% Gets the dimensions of T.
dim = size ( T );

% Calculates the transformation matrix.
B1 = spm_dctmtx ( order (1), dim (1), x );
B2 = spm_dctmtx ( order (2), dim (2), y );
B3 = spm_dctmtx ( order (3), dim (3), z );

d2 = size ( T );
t1 = reshape ( T, d2 (1) * d2 (2), d2 (3) );
t1 = t1 * B3';
t1 = reshape ( t1, d2 (1), d2 (2), [] );
for z = 1: order (3)
    t ( :, :, z ) = B1 * t1 ( :, :, z ) * B2';
end

%==========================================================================
% function p = likelihoods(f,bf,mg,mn,vr)
%==========================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = bsxfun(@minus,cr,mn(:,k)')/chol(vr(:,:,k));
    p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end

% %==========================================================================
% % function dat = decimate(dat,fwhm)
% %==========================================================================
% function dat = decimate(dat,fwhm)
% % Convolve the volume in memory (fwhm in voxels).
% lim = ceil(2*fwhm);
% x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
% y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
% z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
% i  = (length(x) - 1)/2;
% j  = (length(y) - 1)/2;
% k  = (length(z) - 1)/2;
% spm_conv_vol(dat,dat,x,y,z,-[i j k]);

%==========================================================================
% function [P] = clean_gwc(P,level)
%==========================================================================
function [P] = clean_gwc(P,level)
if nargin<4, level = 1; end

b    = P(:,:,:,2);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = double(P(:,:,i,1));
        wp       = double(P(:,:,i,2));
        bp       = double(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = uint8(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
end

% Also clean up the CSF.
if niter2 > 0,
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = double(P(:,:,i,1));
            wp       = double(P(:,:,i,2));
            cp       = double(P(:,:,i,3));
            bp       = double(c(:,:,i))/255;
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = uint8(round(bp));
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
    end
end

th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4),
        slices{k1} = double(P(:,:,i,k1))/255;
    end
    bp        = double(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0,
        cp        = double(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    slices{5} = slices{5}+1e-4; % Add a little to the soft tissue class
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4),
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4),
        P(:,:,i,k1) = uint8(round(slices{k1}./tot*255));
    end 
end
