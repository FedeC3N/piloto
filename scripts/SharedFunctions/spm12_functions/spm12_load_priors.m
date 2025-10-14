function tpm = spm12_load_priors ( tpmname, tpmvols )
% Load the tissue probability maps for segmentation
% FORMAT tpm = spm12_load_priors8(V)
% V   - structures of image volume information (or filenames)
% tpm - a structure for tissue probabilities
%
% This function is intended to be used in conjunction with spm12_sample_priors.
% V = spm12_vol(P);
% T = spm12_load_priors(V);
% B = spm12_sample_priors(T,X,Y,Z);
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm12_load_priors8.m 5962 2014-04-17 12:47:43Z spm $

% Loads the template tissue probabilities MRI.
tpmmri   = ft_read_mri ( tpmname );

% Gets the volumes.
vols     = tpmmri.anatomy ( :, :, :, tpmvols );

% Stores the metadata.
tpm.name = tpmname;
tpm.vols = tpmvols;

% Gets the MRI transformation matrix.
tpm.M    = tpmmri.transform;

% Extracts the backgrouds (average of the first and last slice).
tpm.bg1  = squeeze ( mean ( mean ( vols ( :, :, 1,   : ) ) ) );
tpm.bg2  = squeeze ( mean ( mean ( vols ( :, :, end, : ) ) ) )';

% Calculates the logarithm of the volumes avoiding zeros.
tpm.tiny = 1e-4;
logvols  = log ( vols + tpm.tiny );

% Sets the degree of B-spline.
deg = 1;
tpm.deg  = deg + 1;

% Goes through each volume.
tpm.dat  = cell ( size ( tpmvols ) );
for vindex = 1: numel ( tpmvols )
    
    % Gets the B-spline coefficents of the current volume.
    tpm.dat { vindex } = spm_bsplinc ( logvols ( :, :, :, vindex ), [ deg deg deg  0 0 0 ] );
end

% Stores the MRI.
tpm.mri  = tpmmri;

% Stores metadata of the templates. To fix.
tpm.V.mat = tpmmri.transform;
tpm.V.dim = tpmmri.dim;
