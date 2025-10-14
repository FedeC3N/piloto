function myom_write_geometry ( basename, headmodel )

% Based on OpenMEEG functions:
% * om_save_tri by Alexandre Gramfort
% * om_write_geom by Alexandre Gramfort & Paul Czienskowski
% * om_write_cond by Alexandre Gramfort

% Checks the geometrical definition of the head model.
headmodel = myom_check_headmodel ( headmodel );


% Gets the number of meshes.
nmeshes = numel ( headmodel.bnd );

% Writes the surface files.

% Goes through each mesh.
for mindex = 1: nmeshes
    mesh = headmodel.bnd ( mindex );
    npos = size ( mesh.pos, 1 );
    ntri = size ( mesh.tri, 1 );
    
    % Opens the file to write.
    fid = fopen ( sprintf ( '%s_%s.tri', basename, mesh.tissue ), 'w' );
    
    % Writes the vertex and the triangles.
    fprintf ( fid, '- %g\n', npos );
    fprintf ( fid, '%g %g %g %g %g %g\n', [ mesh.pos , mesh.nrm ]' );
    fprintf ( fid, '- %g %g %g\n', ntri, ntri, ntri );
    fprintf ( fid, '%g %g %g\n', mesh.tri' - 1 );
    
    % Closes the file.
    fclose ( fid );
end

% Writes the geometry file, with the information of all meshes.
% om_write_geom ( sprintf ( '%s.geom', basename ), strcat ( basename, '_', headmodel.tissue, '.tri' ), headmodel.tissue );


% Writes the geometry file.

% Opens the file to write.
fid = fopen ( sprintf ( '%s.geom', basename ), 'w' );

% Writes the header.
fprintf ( fid, '# Domain Description 1.0\n' );
fprintf ( fid, '                        \n' );
fprintf ( fid, 'Interfaces %d Mesh      \n', nmeshes );
fprintf ( fid, '                        \n' );

% Lists the surface files.
for mindex = 1: nmeshes
    mesh = headmodel.bnd ( mindex );
    fprintf ( fid, '%s                  \n', strcat ( basename, '_', mesh.tissue, '.tri' ) );
end

fprintf ( fid, '                        \n' );
fprintf ( fid, 'Domains %d              \n', nmeshes + 1 );
fprintf ( fid, '                        \n' );

% Lists the domains and its surrounding surfaces.
mesh = headmodel.bnd (1);
fprintf ( fid, 'Domain %s -%d\n', mesh.tissue, 1 );
for mindex = 2: nmeshes
    mesh = headmodel.bnd ( mindex );
    fprintf ( fid, 'Domain %s %d -%d\n', mesh.tissue, mindex - 1, mindex );
end
fprintf ( fid, 'Domain Air %d           \n', nmeshes );

% Closes the file.
fclose ( fid );


% Writes the conductivity file.

% Opens the file to write.
fid = fopen ( sprintf ( '%s.cond', basename ), 'w' );

% Writes the header.
fprintf ( fid, '# Properties Description 1.0 (Conductivities)\n\n' );
fprintf ( fid, '%-15s %.16f\n', 'Air', 0 );

% Goes through each mesh.
for mindex = 1: nmeshes
    mesh = headmodel.bnd ( mindex );
    fprintf ( fid, '%-15s %.16f\n', mesh.tissue, headmodel.cond ( mindex ) );
end

% Closes the file.
fclose ( fid );
