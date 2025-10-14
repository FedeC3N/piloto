function headmodel = myom_headmodel ( geometry )
% Based on FiedTrip functions:
% * ft_headmodel_openmeeg


% Adds OpenMEEG to the path.
ft_hastoolbox ( 'openmeeg', 1, 1 );


% Gets the number of meshes.
nmeshes = numel ( geometry.bnd );

% Sorts the meshes from outter to inner.
order = flipud ( ordermesh ( geometry.bnd ) );

if isfield ( geometry, 'bnd' ),    geometry.bnd    = geometry.bnd ( order );    end
if isfield ( geometry, 'tissue' ), geometry.tissue = geometry.tissue ( order ); end
if isfield ( geometry, 'cond' ),   geometry.cond   = geometry.cond ( order );   end


% Initializes the headmodel.
headmodel = [];

% Copies the mesh definition.
headmodel.bnd = geometry.bnd;

% Copies the tissue labels, if provided.
if isfield ( geometry, 'tissue' )
    headmodel.tissue = geometry.tissue;
else
    headmodel.tissue = cellfun ( @(x) sprintf ( 'Tissue %i', x ), num2cell ( 1: nmeshes ), 'UniformOutput', false );
end

% Copies the conductivities, if provided.
if isfield ( geometry, 'cond' )
    headmodel.cond = geometry.cond;
else
    if nmeshes == 1
        headmodel.cond = 1;
    elseif nmeshes == 3
        headmodel.cond = [1 1/80 1] * 0.33;
    else
        error ( 'Conductivity values are required for 2 shells. More than 3 shells not allowed.' )
    end
    warning ( 'No conductivity is declared. Assuming standard values.\n' )
end

% Marks the skin surfaces and the source bounding surface.
headmodel.skin_surface = 1;
headmodel.source = nmeshes;


% Remembers the type of volume conduction model.
headmodel.type = 'openmeeg';

% Operates through each mesh.
for mindex = 1: nmeshes
    
    % For OpenMEEG the normals of the meshes mush be inward-oriented.
    if ~checkomnormals ( headmodel.bnd ( mindex ) )
        
        % Flips the triangle definition.
        fprintf ( 1, 'Flipping the normals in mesh %i to fit OpenMEEG convention.\n', mindex );
        headmodel.bnd ( mindex ).tri = fliplr ( headmodel.bnd ( mindex ).tri );
    end
    
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

% Sets a name for the head model file.
files.hm.name = sprintf ( '%s.mat', tempname );


% Checks the integrity of the OpenMEEG binaries.
om_checkombin;

% Executes 'om_assemble' to get the head model.
status = system ( sprintf ( 'om_assemble -hm "%s" "%s" "%s"\n', files.geom.name, files.cond.name, files.hm.name ) );

% Checks for the completion of the execution.
if status == 0
    
    % Recovers the calculated head model matrix.
    headmodel.hm = importdata ( files.hm.name );
    
% If the execution stopped rises an error.
else
    fprintf ( 2, 'OpenMEEG program ''om_assemble'' exited with error code %i.\n', status );
    
    % Removes all the temporal files and exits.
    cleaner ( files );
    return
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


function correct = checkomnormals ( bnd )

% Gets the vertex and triangles definition.
pnt = bnd.pnt;
tri = bnd.tri;

% Centers the mesh in the origin.
pnt = bsxfun ( @minus, pnt, mean ( pnt ) );

% Gets the solid angle of the mesh.
sa  = sum ( solid_angle ( pnt, tri ) );

% If the solid angle differs from 4pi the mesh is irregular or not closed.
if ( abs ( sa ) - 4 * pi ) > 1000 * eps
    error ( 'The mesh is irregular or not closed.' );
end

% If the solid angle is possitive the normals are correct.
if sa > 0
    correct = true;
else
    correct = false;
end


function order = ordermesh ( bnd )
% Orders the meshes from inner to outter.

% Based on FieldTrip functions:
% * surface_nesting
% * bounding_mesh by Robert Oostenveld
% * solid_angle by Robert Oostenveld


% Defines the variables.
npnt    = 10;
maxiter = 5;

% Gets the metadata.
nbnd    = numel ( bnd );
nesting = zeros ( nbnd );

% Compares each pair of meshes.
for i = 1: nbnd
    for j = 1: nbnd
        
        % If the meshes are the same continues.
        if i == j, continue, end
            
        % Initializes the solid angle.
        sangs = 0;
        
        % Takes several random points of the mesh i.
        for pindex = 1: npnt
            
            % Gets a random point of the mesh i.
            pnt = randi ( size ( bnd ( i ).pnt, 1 ) );
            pnt = bnd ( i ).pnt ( pnt, : );
            
            % Centers the mesh j in the selected point.
            tmp = bsxfun ( @minus, bnd ( j ).pnt, pnt );
            
            % Calculates the solid angle of the mesh j.
            sang  = nansum ( solid_angle ( tmp, bnd ( j ).tri ) );
            sangs = sangs + abs ( sang );
        end
        
        nesting ( i, j ) = round ( sangs / ( npnt * 4 * pi ) );
    end
end

% Gets the order of the meshes.
[ ~, order ] = sort ( -sum ( nesting, 2 ) );

% Checks for ambiguities.
if any ( any ( nesting & nesting' ) )
    
    % Gets the number of iterations.
    stack = dbstack;
    if sum ( strcmp ( { stack.name } , stack ( end ).name ) ) > maxiter
        error ( 'Order not found after %i tries.\n', maxiter )
    end
    
    % Repeats the process iteratively.
    order = NFT_mfc ( bnd );
end
