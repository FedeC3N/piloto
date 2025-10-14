function my_drawNetworkFull ( matrix, config )

hf = figure;
set ( hf, ...
    'Units', 'centimeters', ...
    'Position', [ 0.00 0.00 1.00 0.50 ] * config.width, ...
    'PaperSize', [ 1.00 0.50 ] * config.width );
set ( hf, 'Units', 'normalized' )

% Creates the axes.
ha1 = axes ( hf, 'Position', [ 0.00 0.50 0.25 0.50 ], 'Visible', 'off' );
ha2 = axes ( hf, 'Position', [ 0.00 0.00 0.25 0.50 ], 'Visible', 'off' );
ha3 = axes ( hf, 'Position', [ 0.25 0.50 0.25 0.50 ], 'Visible', 'off' );
ha4 = axes ( hf, 'Position', [ 0.25 0.00 0.25 0.50 ], 'Visible', 'off' );
ha5 = axes ( hf, 'Position', [ 0.50 0.00 0.50 1.00 ], 'Visible', 'off' );


% Draws the four 3D views.
cfg                  = [];
cfg.model            = 'HD';
cfg.ActiveNodeSize   = 10;
cfg.InactiveNodeSize = 5;
cfg.ShowNodesNames   = 'no';
cfg.transparency     = 0.1;
cfg.ConnLineWidth    = 1;
cfg.nodes_positions  = config.nodes;

cfg.handle          = ha1;
Draw_Network ( matrix, cfg );

cfg.handle          = ha2;
Draw_Network ( matrix, cfg );

cfg.handle          = ha3;
Draw_Network ( matrix, cfg );

cfg.handle          = ha4;
Draw_Network ( matrix, cfg );


% Draws the netwrok view.
cfg                 = [];
cfg.groups          = config.groups;
cfg.label           = config.label;
cfg.shadow          = false;
cfg.width           = false;
cfg.handle          = ha5;

my_drawNetworkRing ( matrix, cfg );


% Rotates the 3D views.
view ( ha1,   0,   0 );
view ( ha2, -90,   0 );
view ( ha3,   0,  90 );
view ( ha4,  90,   0 );


% Resizes the network view.
xlim ( ha5, [ -1.7 1.6 ] )
ylim ( ha5, [ -1.6 1.6 ] )
