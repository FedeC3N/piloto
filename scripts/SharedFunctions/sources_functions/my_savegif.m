function my_savegif ( varargin )

% If no axes returns.
if ~numel ( findall ( 0, 'Type', 'axes' ) ), return, end

if ~nargin
    ax       = gca;
    filename = 'file.gif';
end

if nargin == 1
    if ishandle ( varargin {1} )
        ax       = varargin {1};
        filename = 'file.gif';
    else
        ax       = gca;
        filename = varargin {1};
    end
end

if nargin == 2
    ax       = varargin {1};
    filename = varargin {2};
end

if isempty ( ax )
    return
end


% Gets the figure.
fig = get ( ax, 'Parent' );

% Defines the views.
views = linspace ( 0, 360, 100 );
views = views ( 1: end - 1 );

% Rotates the figure to the first position.
view ( ax, views (1), 0 );

% Defines the map.
frame = getframe ( fig );
im    = frame2im ( frame );
[ ~, cm ] = rgb2ind ( im, 256 );

% Gets the figure bitmap.
frame = getframe ( fig );
im    = frame2im ( frame );
imind = rgb2ind  ( im, cm );

% Writes the first slice.
imwrite ( imind, cm, filename, 'gif', 'DelayTime', .04, 'Loopcount', inf );

% Goes through each slice.
for vindex = 2: numel ( views )
    
    % Rotates the figure.
    view ( ax, views ( vindex ), 0 );
    
    % Gets the figure bitmap.
    frame = getframe ( fig );
    im    = frame2im ( frame );
    imind = rgb2ind  ( im, cm );
    
    % Writes the figure as a GIF slice.
    imwrite ( imind, cm, filename, 'gif', 'DelayTime', .04, 'WriteMode', 'append' );
end
