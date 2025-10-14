function skull = BS_segmentSkull ( mridata, aseg )
% Uses BrainSuite to segment the subject's skull.

% Generates a temporal name for the intermediate files.
basename = tempname;


% Saves the original MRI.
mri   = mridata.anatomy;
mri   = double ( mri );

ft_write_mri ( sprintf ( '%s_mri.nii', basename ), mri, 'dataformat', 'nifti' );


fprintf ( 1, 'Creating the brain mask.\n' );

% Generates a brain mask from the FreeSurfer segmentation.
brain = aseg.anatomy ~= 0;
brain = imclose ( brain, ones ( 5,5,5 ) );
brain = imfill  ( brain, 'holes' );
brain = 255 * double ( brain );

ft_write_mri ( sprintf ( '%s_brainmask.nii', basename ), brain, 'dataformat', 'nifti' );


fprintf ( 1, 'Segmentating the skull.\n' );

% Performs the skull segmentation and filters the resulting mask.
system ( sprintf ( './functions/skullfinder -i %s_mri.nii -m %s_brainmask.nii -o %s_skull.nii', basename, basename, basename ) );
% system ( sprintf ( './functions/scrubmask -i %s_skull.nii -o %s_skull.nii', basename, basename ) );


fprintf ( 1, 'Fixing errors in the skull mask.\n' );

% Loads the skull (and inner tissues) mask.
skull = ft_read_mri ( sprintf ( '%s_skull.nii', basename ) );
skull = skull.anatomy >= 17;

% Fix the segmented skull to include the whole brain and avoid the scalp.
brain = mridata.brain;
scalp = mridata.scalp;

brain = imdilate ( brain, ones ( 5, 5, 5 ) );
scalp = imerode  ( scalp, ones ( 5, 5, 5 ) );

skull = skull | brain;
skull = skull & scalp;

% Performs a morphological closing to fix holes and imprefections.
skull = imclose ( skull, ones ( 15, 15, 15 ) );

skull = skull | brain;
skull = skull & scalp;

% Removes the temporal files.
delete ( sprintf ( '%s*', basename ) );
