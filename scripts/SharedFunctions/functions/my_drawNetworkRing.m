function drawNetwork ( matrix, cfg )

% drawNetwork ( matrix, cfg )
% 
% 'config' is an struture with fields:
% * handle (axes handler). Handle of the axes where the plot must be drawn.
% * label (nx1 cell). The name of the nodes.
% * width (number or false). Width of the lines.
% * shadow (true or false). To draw the shadow of the lines.

if nargin == 0
    matrix = randn ( 21 );
    matrix ( abs ( matrix ) < 2 ) = 0;
end

if nargin < 2
    cfg = struct;
end

if ~isfield ( cfg, 'label' )
    label = cellfun ( @(x) sprintf ( 'Area %2.0f', x ), num2cell ( 1: length ( matrix ) ), 'UniformOutput', false );
else
    label = cfg.label;
end


% Gets the number of nodes.
nodes  = size ( matrix, 1 );

% Draws the N nodes linearly spaced over a circunference.
angles = linspace ( pi, -pi, nodes + 1 );
angles = angles ( 1: end - 1 ) - pi / 2 - pi / nodes;

if isfield ( cfg, 'handle' )
    ahandle = cfg.handle;
else
    figure
    set ( gcf, 'Position', get ( gcf, 'Position' ) .* [ 1 1 0 0 ] + [ 0 0 500 500 ] );
    ahandle = axes ( gcf );
end
hold ( ahandle, 'on' )
axis ( ahandle, 'equal', 'off' )
xlim ( ahandle, [ -1.1 1.1 ] )
ylim ( ahandle, [ -1.1 1.1 ] )

% Draws the circle.
th = linspace ( -pi, pi, 360 );
x  = cos ( th );
y  = sin ( th );

plot ( ahandle, x, y, 'Color', [ .9 .9 .9 ] );


% Gets the list of active conections.
[ node1, node2 ] = find ( matrix );

% Removes the diagonal and the lower triangular.
remove = node1 >= node2;
node1  = node1 ( ~remove );
node2  = node2 ( ~remove );

% Goes through each conexion.
for lindex = 1: numel ( node1 )
    
    % Defines the extremes of the arc.
    x1 = cos ( angles ( node1 ( lindex ) ) );
    y1 = sin ( angles ( node1 ( lindex ) ) );
    x2 = cos ( angles ( node2 ( lindex ) ) );
    y2 = sin ( angles ( node2 ( lindex ) ) );
    
    % Finds the center and radius of the arc.
    x0 = -( y1 - y2 ) / ( x1 * y2 - y1 * x2 );
    y0 =  ( x1 - x2 ) / ( x1 * y2 - y1 * x2 );
    r  = sqrt ( x0 ^ 2 + y0 ^ 2 - 1 );
    
    % If the radius is greater than 1000 writes a rect line.
    if r > 1e3
        
        % Draws the shadow of the rect, if required.
        if cfg.shadow
            shandle = line ( ahandle, linspace ( x1, x2, 50 ), linspace ( y1, y2, 50 ) );
        end
        
        % Draws the rect.
        lhandle = line ( ahandle, linspace ( x1, x2, 50 ), linspace ( y1, y2, 50 ) );
        
    % Otherwise calculates the arc.
    else
        
        % Defines the extreme angles of the arc.
        th1 = atan2 ( y1 - y0, x1 - x0 );
        th2 = atan2 ( y2 - y0, x2 - x0 );
        
        if x1 >= 0 && x2 >= 0
            % ensure the arc is within the unit disk
            theta = [linspace(max(th1,th2),pi,50),...
                linspace(-pi,min(th1,th2),50)].';
        else
            theta = linspace ( th1, th2 );
        end
        
        % Draws the shadow of the arc, if required.
        if cfg.shadow
            shandle = line ( ahandle, r * cos ( theta ) + x0, r * sin ( theta ) + y0 );
        end
        
        % Draws the arc.
        lhandle = line ( ahandle, r * cos ( theta ) + x0, r * sin ( theta ) + y0 );
    end
    
    % Modifies the width, if required.
    if cfg.width
        set ( lhandle, 'LineWidth', cfg.width );
    end
    
    % Sets the link color.
    if matrix ( node1 ( lindex ), node2 ( lindex ) ) > 0
        
        if cfg.shadow
            set ( shandle, 'Color', [ 1 .7 .7 ], 'LineWidth', 3 );
        end
        
        set ( lhandle, 'Color', [ 1 0 0 ] );
    else
        
        if cfg.shadow
            set ( shandle, 'Color', [ .7 .7 1 ], 'LineWidth', 3 );
        end
        
        set ( lhandle, 'Color', [ 0 0 1 ] );
    end
    
    % If version greater than 8.4 uses transparency.
    if cfg.shadow && ~verLessThan ( 'matlab', '8.4' )
        set ( shandle, 'Color', cat ( 2, get ( lhandle, 'Color' ), .3 ), 'LineWidth', 3 );
    end
end


% Draws the nodes over the link lines.
for nindex = 1: nodes
    
    % Gets the x and y coordinates.
    x = cos ( angles ( nindex ) );
    y = sin ( angles ( nindex ) );
    nhandle = plot ( ahandle, x, y );
    
    set ( nhandle, 'Marker', 'o' );
    set ( nhandle, 'MarkerEdgeColor', [ 0 0 0 ] );
    set ( nhandle, 'MarkerFaceColor', [ 0 0 0 ] );
    
    % Creates the text label.
    lhandle = text ( ahandle, x * 1.1, y * 1.1, label { nindex } );
    
    % Rotates the text label.
    % Clockwise for negative values of x.
    if x < 0
        set ( lhandle, 'HorizontalAlignment', 'right' );
        set ( lhandle, 'Rotation', 180 * angles ( nindex ) / pi - 180 );
        set ( lhandle, 'FontSize', 9 );
        
    % Counterclockwise otherwise.
    else
        set ( lhandle, 'HorizontalAlignment', 'left' );
        set ( lhandle, 'Rotation', 180 * angles ( nindex ) / pi );
        set ( lhandle, 'FontSize', 9 );
    end
end
