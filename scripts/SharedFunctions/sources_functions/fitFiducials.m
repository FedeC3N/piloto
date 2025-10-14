function transform = fitFiducials ( static, moving )

% Checks the inputs.
if nargin ~= 2
    error ( 'This function requires two input arguments.' );
end
if ...
        ~isfield ( static,     'fid' )   || ...
        ~isfield ( static.fid, 'label' ) || ...
        ~isfield ( static.fid, 'pos' )   || ...
        ~isfield ( moving,     'fid' )   || ...
        ~isfield ( moving.fid, 'label' ) || ...
        ~isfield ( moving.fid, 'pos' )
    error ( 'This functions requires fiducial data as input.' );
end
if size ( static.fid.pos, 2 ) ~= 3 || size ( moving.fid.pos, 2 ) ~= 3
    error( 'Fiducial position is not 3D.' );
end

% Gets the list of shared fiducials.
fiducials = intersect ( static.fid.label, moving.fid.label );

% If less than 3 shared fiducials rises an error.
if numel ( fiducials ) < 3
    error ( 'Both inputs must share at least 3 fiducials.' );
end

% Initializes the coulds of points.
x = zeros ( numel ( fiducials ), 3 );
y = zeros ( numel ( fiducials ), 3 );

% Gets the shared fiducials' position.
for findex = 1: numel ( fiducials )
    x ( findex, : ) = static.fid.pos ( strcmp ( static.fid.label, fiducials ( findex ) ), : );
    y ( findex, : ) = moving.fid.pos ( strcmp ( moving.fid.label, fiducials ( findex ) ), : );
end

% Fits both clouds of points.
transform = fitPoints ( x, y );
