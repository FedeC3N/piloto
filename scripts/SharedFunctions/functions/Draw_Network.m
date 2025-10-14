function handle = Draw_Network ( matrix, cfg )

% Draw_Network ( matrix, config )
% 
% 'config' is an struture with fields:
% * nodes_positions (nx3 matrix). MNI position of the nodes in mm
%   (required).
% * handle (axes handler). Handle of the axes where the plot must be drawn.
% * model ('HD', 'LD' or none). Cortex mesh to be plotted.
% * title (string). Title to be showed in the canvas.
% * label (1x2 cell string). Labels of the groups compared.
% * namout (string).
% * ShowNodesNames ('yes' or 'no').
% * NodesNames (nx1 cell).
% * FontLabelSize (number).
% * ShowROIsNames ('yes' or 'no').
% * ColorPositiveConn (1x3 RGB vector).
% * ColorNegativeConn (1x3 RGB vector).
% * transparency (number).
% * azimutal (number).
% * elevation (number).
% * luminosity (number).
% * BackgroundColor (1x3 RGB vector).
% * IntraROInodesize (number).
% * ActiveNodeSize (number).
% * InactiveNodeSize (number).
% * ActiveNodesColor (1x3 RGB vector).
% * InactiveNodesColor (1x3 RGB vector).
% * FixedConnLinesWidth ('yes' or 'no'). All links has the same width.
% * ConnLineWidth (number). Fixed with of the link.
% * DrawPopulation ('yes' or 'no'). Uses -log ( matrix ) as power metric.
%   Useful when matrix is the p-value matrix.


% Checks the input parameters.
if nargin < 2, error ( 'This program requires a configuration structure.' ); end
if ~isfield ( cfg, 'nodes_positions' ), error ( 'This program requires the node''s position.' ); end

% Selects the brain model.
model = importdata('ModelBrainMNI.mat');
if isfield ( cfg, 'model' )
    if strcmp(cfg.model,'LD')
        Faces = model.BrainLD.Faces;
        Vertices = model.BrainLD.Vertices*1.15;
    elseif strcmp(cfg.model,'HD')
%         Faces = model.BrainHD.Faces;
%         Vertices = model.BrainHD.Vertices*1.15;
        Faces = model.BrainDavid.Faces;
        Vertices = model.BrainDavid.Vertices;
%         Faces = model.BrainDSmooth.Faces;
%         Vertices = model.BrainDSmooth.Vertices;
%         Vertices = bsxfun ( @minus, Vertices, [ 0 20 0 ] );
%         Vertices = bsxfun ( @plus, Vertices, [ 0 0 15 ] );
    else
        Faces = model.VolCon.Faces;
        Vertices = model.VolCon.Vertices;
    end
else
    Faces = model.VolCon.Faces;
    Vertices = model.VolCon.Vertices;
end


% Cheacks and initializes the configuration.
if ~isfield ( cfg, 'title' ),               cfg.title               = ''; end
if ~isfield ( cfg, 'label' ),               cfg.label               = { 'Group 1' 'Group 2' }; end
if ~isfield ( cfg, 'nameout' ),             cfg.nameout             = 'Test'; end

if ~isfield ( cfg, 'ShowNodesNames' ),      cfg.ShowNodesNames      = 'yes'; end
if ~isfield ( cfg, 'NodesNames' ),          cfg.NodesNames          = cellstr ( num2str ( ( 1: size ( matrix, 1 ) )' ) ); end
if ~isfield ( cfg, 'FontLabelSize' ),       cfg.FontLabelSize       = 10; end

%   Show ROIs names (1 name per ROI useful when many nodes are depicted)
if ~isfield ( cfg, 'ShowROIsNames' ),       cfg.ShowROIsNames       = 'no'; end

if ~isfield ( cfg, 'ColorPositiveConn' ),   cfg.ColorPositiveConn   = [ .9  .0  .0 ]; end
if ~isfield ( cfg, 'ColorNegativeConn' ),   cfg.ColorNegativeConn   = [ .0  .0  .9 ]; end
if ~isfield ( cfg, 'transparency' ),        cfg.transparency        = 0.65; end
if ~isfield ( cfg, 'azimutal' ),            cfg.azimutal            = 0; end
if ~isfield ( cfg, 'elevation' ),           cfg.elevation           = 90; end
if ~isfield ( cfg, 'luminosity' ),          cfg.luminosity          = 0.6; end
if ~isfield ( cfg, 'BackgroundColor' ),     cfg.BackgroundColor     = 'w'; end

if ~isfield ( cfg, 'IntraROInodesize' ),    cfg.IntraROInodesize    = 35; end
if ~isfield ( cfg, 'ActiveNodeSize' ),      cfg.ActiveNodeSize      = 25; end
if ~isfield ( cfg, 'InactiveNodeSize' ),    cfg.InactiveNodeSize    = 10; end
if ~isfield ( cfg, 'ActiveNodesColor' ),    cfg.ActiveNodesColor    = [ .2 .2 .2 ]; end
if ~isfield ( cfg, 'InactiveNodesColor' ),  cfg.InactiveNodesColor  = [ .7 .7 .7 ]; end

if ~isfield ( cfg, 'FixedConnLinesWidth' ), cfg.FixedConnLinesWidth = 'yes'; end
if ~isfield ( cfg, 'ConnLineWidth' ),       cfg.ConnLineWidth       = 2; end

if ~isfield ( cfg, 'DrawPopulation' ),      cfg.DrawPopulation      = 'no'; end


% Creates the figure, if no handle provided.
if ~isfield ( cfg, 'handle' )
    figure ( 'Color', cfg.BackgroundColor )
    handle = gca;
else
    handle = cfg.handle;
end

% Sets the properties of the canvas.
axes ( handle )
hold on
axis equal off

% Plots two points to allow the legend to work.
plot ( nan, 'Color', cfg.ColorPositiveConn, 'DisplayName', sprintf ( '%s > %s', cfg.label {:} ) );
plot ( nan, 'Color', cfg.ColorNegativeConn, 'DisplayName', sprintf ( '%s < %s', cfg.label {:} ) );

% Draws the brain.
trisurf ( Faces, Vertices ( :, 1 ), Vertices ( :, 2 ), Vertices ( :, 3 ), ...
    'SpecularStrength', 0.2, ...
    'DiffuseStrength',  0.8,'AmbientStrength',0.5, ...
    'FaceLighting',     'phong', ...
    'LineStyle',        'none', ...
    'FaceColor',        [ .75 .75 .75], ...
    'FaceAlpha',        cfg.transparency, ...
    'EdgeColor',        'none', ...
    'tag',              'cortex' );

% Hides the brain in the legend.
set ( get ( get ( findall ( gca, 'tag', 'cortex' ), 'Annotation' ), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

% Lights the scene.
light ( 'Position',[  0  0  1 ], 'color', cfg.luminosity * [ 1 1 1 ] );
light ( 'Position',[  0  1  0 ], 'color', cfg.luminosity * [ 1 1 1 ] );
light ( 'Position',[  0 -1  0 ], 'color', cfg.luminosity * [ 1 1 1 ] );
lighting phong
material dull

% Sets the viewing point.
view ( cfg.azimutal, cfg.elevation )


% Forces the matrix to be simetric.
matrix = triu ( matrix ) + triu ( matrix, 1 )' - diag ( diag ( matrix ) );

% Gets the nodes positions.
nodes = cfg.nodes_positions;
nodes ( :, 4 ) = nan;


% Get the list of significant links.
[ lp1, lp2 ] = find ( matrix ~= 0 );

% Gets the list of active nodes.
anodes = any ( matrix );

% Gets the power of the significant links.
power = matrix ( matrix ~= 0 );

% Gets the intra-ROI connectivity value.
ipower = diag ( matrix );


% Craws all the links.
hlinks = line ( cat ( 2, nodes ( lp1, 1 ), nodes ( lp2, 1 ) )',  cat ( 2, nodes ( lp1, 2 ), nodes ( lp2, 2 ) )',  cat ( 2, nodes ( lp1, 3 ), nodes ( lp2, 3 ) )', 'LineStyle', '-', 'LineWidth', cfg.ConnLineWidth );
% hshade = line ( cat ( 2, nodes ( lp1, 1 ), nodes ( lp2, 1 ) )',  cat ( 2, nodes ( lp1, 2 ), nodes ( lp2, 2 ) )',  cat ( 2, nodes ( lp1, 3 ), nodes ( lp2, 3 ) )', 'LineStyle', '-', 'LineWidth', cfg.ConnLineWidth * 2 );

% Sets the link color, according with its sign.
set ( hlinks ( power > 0 ), 'Color', cfg.ColorPositiveConn );
set ( hlinks ( power < 0 ), 'Color', cfg.ColorNegativeConn );

% set ( hshade ( power > 0 ), 'Color', cat ( 2, cfg.ColorPositiveConn, .3 ) );
% set ( hshade ( power < 0 ), 'Color', cat ( 2, cfg.ColorNegativeConn, .3 ) );

% Modifies the link width, if needed.
if strcmp ( cfg.FixedConnLinesWidth, 'no' )
    
    % Log-scalates the power, if needed.
    if strcmp ( cfg.DrawPopulation, 'yes' )
        power = ( -log ( abs ( power ) ) + eps ) .* sign ( power );
    end
    
    % Goes through each link.
    for lindex = 1: numel ( hlinks )
        set ( hlinks ( lindex ), 'LineWidth', power ( lindex ) );
    end
end

% Hides the links in the legend.
for lindex = 1: numel ( hlinks )
    set ( get ( get ( hlinks ( lindex ), 'Annotation' ), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
%     set ( get ( get ( hshade ( lindex ), 'Annotation' ), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
end


% Draws the nodes.
hnodes = plot3 ( nodes ( :, [ 1 4 ] )', nodes ( :, [ 2 4 ] )', nodes ( :, [ 3 4 ] )', ...
    'LineStyle', 'none', ...
    'Marker', '.' );

% Hides the nodes in the legend.
for nindex = 1: numel ( hnodes )
    set ( get ( get ( hnodes ( nindex ), 'Annotation' ), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
end


% Colours the nodes.
set ( hnodes ( ~anodes ), ...
    'MarkerSize',      cfg.InactiveNodeSize, ...
    'MarkerEdgeColor', cfg.InactiveNodesColor, ...
    'MarkerFaceColor', cfg.InactiveNodesColor );
set ( hnodes (  anodes ), ...
    'MarkerSize',      cfg.ActiveNodeSize, ...
    'MarkerEdgeColor', cfg.ActiveNodesColor, ...
    'MarkerFaceColor', cfg.ActiveNodesColor );

% Modifies the nodes size and color depending on the connectivity.
set ( hnodes ( ipower > 0 ), ...
    'MarkerSize',      cfg.IntraROInodesize, ...
    'MarkerEdgeColor', cfg.ColorPositiveConn, ...
    'MarkerFaceColor', cfg.ColorPositiveConn );
set ( hnodes ( ipower < 0 ), ...
    'MarkerSize',      cfg.IntraROInodesize, ...
    'MarkerEdgeColor', cfg.ColorNegativeConn, ...
    'MarkerFaceColor', cfg.ColorNegativeConn );


% Writes a label for each active node.
if strcmp ( cfg.ShowNodesNames, 'yes' )
    
    % Gets the list of active nodes.
    anodes = any ( matrix );
    
    % Writes the labels for the active nodes.
    text ( nodes ( anodes, 1 ), nodes ( anodes, 2 ), nodes ( anodes, 3 ), cfg.NodesNames ( anodes ), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'middle', ...
        'Color',               'k',      ...
        'FontSize',            cfg.FontLabelSize );
end

% Writes out a label for each ROI.
if strcmp ( cfg.ShowROIsNames, 'yes' )
    
    % Writes the label for each ROI.
    text ( cfg.rois_positions ( :, 1, 1, 1 ) * 1.1, cfg.rois_positions ( :, 2, 1, 1 ) * 1.1, cfg.rois_positions ( :, 3, 1, 1 ) * 1.1, cfg.rois_names, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'middle', ...
        'Color',               'k',      ...
        'FontSize',            cfg.FontLabelSize );
end

% Sets the figure title.
title ( cfg.title )
