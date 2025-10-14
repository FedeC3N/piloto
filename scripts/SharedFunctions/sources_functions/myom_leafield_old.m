function grid = myom_leafield ( headmodel, grid, sens )

% Based on FiedTrip functions:
% * ft_leadfield_openmeeg by Daniel D.E. Wong, Sarang S. Dalal
%
% Based on the OpenMEEG functions:
% * openmeeg_dsm by Alexandre Gramfort
% * openmeeg_megm by Emmanuel Olivi


% Adds OpenMEEG to the path.
ft_hastoolbox ( 'openmeeg', 1, 1 );


% Gets the sensors type.
ismeg = isfield ( sens, 'coilpos' ) &&  isfield ( sens, 'coilori' );
iseeg = isfield ( sens, 'elecpos' );

% Gets sure that the sensors are correctly identified.
if ~xor ( ismeg, iseeg ), error ( 'The sensor type could not be identified as EEG or MEG. Aborting.\n' ); end


% If no headmodel matrix, calculates it using myom_headmodel.
if ~isfield ( headmodel, 'hm' )
    fprintf ( 'Head model matrix not present in the data. Calculating it.\n' );
    headmodel = myom_headmodel ( headmodel );
end

% Transforms the stored head model structure to a symmetric matrix.
headmodel.hm = struct2sym ( headmodel.hm );


% If there is an 'inside' field in the grid, takes only those grid points.
if isfield ( grid, 'inside' )
    dipoles = grid.pos ( grid.inside, : );
else
    dipoles = grid.pos;
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
dipoles  = cat ( 2, kron ( dipoles, ones ( 3, 1 ) ), kron ( ones ( ndipoles, 1 ), eye (3) ) );

% Writes the dipole file.
files.dip.name = sprintf ( '%s.bin', tempname );
om_save_full ( dipoles, files.dip.name );


% Gets the number of sensors.
nsens = numel ( sens.label );

% Writes the sensors position.
if ismeg
    
    % Writes the gradiometers file.
    files.sens.name = sprintf ( '%s.txt', tempname );
    om_save_full ( cat ( 2, sens.coilpos, sens.coilori ), files.sens.name, 'ascii' );
end
if iseeg
    
    % Writes the electrodes file.
    files.sens.name = sprintf ( '%s.txt', tempname );
    om_save_full ( sens.elecpos, files.sens.name, 'ascii' );
end


% Sets a name for the output files.
files.dsm.name  = sprintf ( '%s.mat', tempname );
files.s2mm.name = sprintf ( '%s.mat', tempname );
files.h2mm.name = sprintf ( '%s.mat', tempname );
files.h2em.name = sprintf ( '%s.mat', tempname );


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

if ismeg
    
    % Executes 'om_assemble' to get the dipoles to MEG matrix.
    status = system ( sprintf ( 'om_assemble -ds2mm "%s" "%s" "%s"\n', files.dip.name, files.sens.name, files.s2mm.name ) );
    
    % Checks for the completion of the execution.
    if status == 0
        
        % Recovers the calculated dipoles model matrix.
        headmodel.s2mm = importdata ( files.s2mm.name );
        
        % If the execution stopped rises an error.
    else
        fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );
        
        % Removes all the temporal files and exits.
        cleaner ( files );
        return
    end
    
    % Executes 'om_assemble' to get the head surface to MEG matrix.
    status = system ( sprintf ( 'om_assemble -h2mm "%s" "%s" "%s" "%s"\n', files.geom.name, files.cond.name, files.sens.name, files.h2mm.name ) );
    
    % Checks for the completion of the execution.
    if status == 0
        
        % Recovers the calculated dipoles model matrix.
        headmodel.h2mm = importdata ( files.h2mm.name );
        
    % If the execution stopped rises an error.
    else
        fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );
        
        % Removes all the temporal files and exits.
        cleaner ( files );
        return
    end
    
    % Calculates the leadfield.
    fprintf ( 1, 'Building the leadfield matrix.\n' );
    leadfield = headmodel.s2mm + headmodel.h2mm * ( headmodel.hm \ headmodel.dsm );
end

if iseeg
    
    % Executes 'om_assemble' to get the head surface to EEG matrix.
    status = system ( sprintf ( 'om_assemble -h2em "%s" "%s" "%s" "%s"\n', files.geom.name, files.cond.name, files.sens.name, files.h2em.name ) );
    
    % Checks for the completion of the execution.
    if status == 0
        
        % Recovers the calculated dipoles model matrix.
        headmodel.h2em = importdata ( files.h2em.name );
        
    % If the execution stopped rises an error.
    else
        fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );
        
        % Removes all the temporal files and exits.
        cleaner ( files );
        return
    end
    
    % Calculates the leadfield.
    fprintf ( 1, 'Building the leadfield matrix.\n' );
    leadfield = headmodel.h2em * ( headmodel.hm \ headmodel.dsm );
end

% Applies the coil/electrode to channel transformation.
if isfield ( sens, 'tra' )
    leadfield = sens.tra * leadfield;
end

% Puts the leadfield in FieldTrip form (cell array of Nx3 matrices).
leadfield = mat2cell ( leadfield, nsens, 3 * ones ( 1, ndipoles ) );


% Stores the leadfield in the grid.
if isfield ( grid, 'inside' )
    grid.leadfield = cell ( 1, size ( grid.pos, 1 ) );
    grid.leadfield ( grid.inside ) = leadfield;
else
    grid.leadfield = leadfield;
end

% Adds the channel labels.
grid.label = sens.label;

% Adds the dimension order in FieldTrip format.
grid.leadfielddimord = '{pos}_chan_ori';

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
