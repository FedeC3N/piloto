function mri = spm12_segment ( mri )

% Includes SPM8 in the path.
ft_hastoolbox ( 'spm8', 1, 1 );

% The MRI can be a 4D matrix with several channels.
% All 3D matrices must share dimensions and transformations.

% Options for the affine transformation to TPM space.
cfg1          = [];
cfg1.fwhm1    = 16;
cfg1.fwhm2    = 0;
cfg1.samp     = 3;

% Options for the preprocessing (warp transformation to TPM space).
cfg2          = [];
cfg2.fwhm     = 0;
cfg2.samp     = 3;
cfg2.biasreg  = 1e-3;
cfg2.biasfwhm = 60;
cfg2.lkp      = [ 1 2 3 3 4 4 4 5 5 5 5 6 6 ];
cfg2.reg      = [ 0 1e-3 .5 .05 .2 ];

% Options to the tissues saving.
cfg3          = [];
cfg3.tissues  ( :, 1 ) = [ 1 1 1 0 0 0 ];
cfg3.tissues  ( :, 2 ) = [ 0 0 0 0 0 0 ];
cfg3.tissues  ( :, 3 ) = [ 0 0 0 0 0 0 ];
cfg3.tissues  ( :, 4 ) = [ 0 0 0 0 0 0 ];
cfg3.tissues  = true ( 6, 4 );
cfg3.wbf      = true ( 1, 2 );
cfg3.wdf      = true ( 1, 2 );
cfg3.tissues  = false ( 6, 4 );
cfg3.wbf      = false ( 1, 2 );
cfg3.wdf      = false ( 1, 2 );


% Loads the template (tissue probability map) as an a-priori structure.
tpm           = spm12_load_priors ( sprintf ( '%s/TPM.nii', fileparts ( mfilename ( 'fullpath' ) ) ), 1: 6 );

% Makes sure that the MRI is in single precission.
mri.anatomy   = single ( mri.anatomy );

% Performs a two-step affine registration of the MRI image to MNI space.
afftrans      = eye (4);
afftrans      = spm12_maff ( mri, cfg1.samp, cfg1.fwhm1, tpm, afftrans, 'mni' );
afftrans      = spm12_maff ( mri, cfg1.samp, cfg1.fwhm2, tpm, afftrans, 'mni' );

% Segments the MRI image in the requested tissues.
preprocinfo   = spm12_preproc       ( mri, tpm, afftrans, cfg2 );

% Writes out the tissue probability maps.
preproc       = spm12_preproc_write ( preprocinfo, cfg3.tissues, cfg3.wbf, cfg3.wdf );


% Replaces the original MRI volume with the bias-corrected one.
mri.anatomy   = preproc.mri;

% Includes the MNI-transformation data in the MRI.
mri.mni2sub   = preproc.mni2sub;

% Includes the segmentation data in the MRI.
mri.gray      = preproc.tc1 {1};
mri.white     = preproc.tc1 {2};
mri.csf       = preproc.tc1 {3};
mri.bone      = preproc.tc1 {4};
