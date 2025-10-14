function headmodel = myom_sourcematrix ( headmodel, grid )

% Based on FiedTrip functions:
% * ft_leadfield_openmeeg by Daniel D.E. Wong, Sarang S. Dalal
%
% Based on the OpenMEEG functions:
% * openmeeg_dsm by Alexandre Gramfort
% * openmeeg_megm by Emmanuel Olivi


% Adds OpenMEEG to the path.
ft_hastoolbox ( 'openmeeg', 1, 1 );


% If there is an 'inside' field in the grid, takes only those grid points.
if isfield ( grid, 'inside' )
    dipoles = grid.pos ( grid.inside, : );
else
    dipoles = grid.pos;
end

% If no original dipole orientation, uses the identity matrix.
if ~isfield ( grid, 'ori' )
    grid.ori = eye (3);
end

% Gets the number of inside dipoles.
ndipoles = size ( dipoles, 1 );


% Gets the number of meshes.
nmeshes = numel ( headmodel.bnd );

% Operates through each mesh.
for mindex = 1: nmeshes
    
    % Saves a triangles file for the current mesh.
    files.tri ( mindex ).name = sprintf ( '%s.tri', tempname );
    om_save_tri ( files.tri ( mindex ).name, headmodel.bnd ( mindex ).pnt, headmodel.bnd ( mindex ).tri );
end

% Writes the geometry file, with the information of all meshes.
files.geom.name = sprintf ( '%s.geom', tempname );
om_write_geom ( files.geom.name, { files.tri.name }, headmodel.tissue );

% Writes the conductivity file.
files.cond.name = sprintf ( '%s.cond', tempname );
om_write_cond ( files.cond.name, headmodel.cond, headmodel.tissue );


% Triplicates the dipole positions, one repetition for each orientation.
dipoles  = cat ( 2, kron ( dipoles, ones ( 3, 1 ) ), kron ( ones ( ndipoles, 1 ), grid.ori ) );

% Writes the dipole file.
files.dip.name = sprintf ( '%s.bin', tempname );
om_save_full ( dipoles, files.dip.name );


% Sets a name for the output file.
files.dsm.name  = sprintf ( '%s.mat', tempname );


% Checks the integrity of the OpenMEEG binaries.
om_checkombin;

% Executes 'om_assemble' to get the dipoles matrix.
status = system ( sprintf ( 'om_assemble -dsm "%s" "%s" "%s" "%s" "%s"\n', files.geom.name, files.cond.name, files.dip.name, files.dsm.name, headmodel.tissue { end } ) );

% Checks for the completion of the execution.
if status == 0
    
    % Recovers the calculated dipoles model matrix.
    headmodel.dsm = importdata ( files.dsm.name );
    
% If the execution stopped rises an error.
else
    fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );

    % Removes all the temporal files and exits.
    cleaner ( files );
    return
end

% If the headmodel matrix is present calculates hm_dsm.
if isfield ( headmodel, 'hm' )
    
    fprintf ( 1, 'Calculating inv ( hm ) * dsm.\n' );
    
    % Transforms the head model matrix to a symmetric matrix.
    if isstruct ( headmodel.hm )
        headmodel.hm = struct2sym ( headmodel.hm );
    end
    
    % Calculates inv ( hm ) * dsm.
    headmodel.hm_dsm = headmodel.hm \ headmodel.dsm;
    
    % Deletes the head model and sources matrices to save memory.
    headmodel = rmfield ( headmodel, { 'hm' 'dsm' } );
end

% Removes all the temporal files.
cleaner ( files );



function cleaner ( files )

% Gets the list of file types.
filetypes = fieldnames ( files );

% Goes through each file type.
for tindex = 1: numel ( filetypes )
    
    filetype = filetypes { tindex };
    
    % Removes all the files of the current type.
    for findex = 1: numel ( files.( filetype ) )
        if exist ( files.( filetype ) ( findex ).name, 'file' )
            delete ( files.( filetype ) ( findex ).name )
        end
    end
end


function matrix = struct2sym ( structure )

% Generates a matrix of the given size.
matrix = zeros ( structure.size );

% Copies the upper diagonal from the structure data to the matrix.
matrix ( triu ( true ( structure.size ) ) ) = structure.data;

% Fills the lower diagonal transposing and removing the diagonal.
matrix = matrix +  matrix' - diag ( diag ( matrix ) );
