function mri = my_read_mri ( filename )

if strcmp ( filename ( end - 2: end ), '.gz' )
    rmfile = true;
    filename = gunzip ( filename, tempdir );
    filename = filename {1};
else
    rmfile = false;
end

mri = ft_read_mri ( filename );

if rmfile
    delete ( filename )
end