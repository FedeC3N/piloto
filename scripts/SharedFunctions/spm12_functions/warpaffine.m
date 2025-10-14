function [ tx, ty, tz ] = warpaffine ( x, y, z, affine, warp )

% If x, y and z are vectors constructs the grid.
if ~isequal ( size ( x ), size ( y ), size ( z ) )
    x  = x ( :, ones ( size ( y, 2 ), 1 ), ones ( size ( z, 3 ), 1 ) );
    y  = y ( ones ( size ( x, 1 ), 1 ), :, ones ( size ( z, 3 ), 1 ) );
    z  = z ( ones ( size ( x, 1 ), 1 ), ones ( size ( y, 2 ), 1 ), : );
end

% If no warp, uses zeros.
if nargin == 5
    
    % Splits the warp per dimension.
    xw = warp ( :, :, :, 1 );
    yw = warp ( :, :, :, 2 );
    zw = warp ( :, :, :, 3 );
    
    % Applies the nonlinear warp.
    x  = x + xw;
    y  = y + yw;
    z  = z + zw;
end

% Applies the affine transformation matrix.
tx = affine (1,1) * x + affine (1,2) * y + affine (1,3) * z + affine (1,4);
ty = affine (2,1) * x + affine (2,2) * y + affine (2,3) * z + affine (2,4);
tz = affine (3,1) * x + affine (3,2) * y + affine (3,3) * z + affine (3,4);

% If only one output argument, concatenates the three 3D matrices.
if nargout == 1
    tx = cat ( 4, tx, ty, tz );
end
