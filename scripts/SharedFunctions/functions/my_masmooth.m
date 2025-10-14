function sdata = my_masmooth ( data, span )

% Forces the averaging width to be simmetric.
width = span - 1 + mod ( span, 2 );

% Calculates the starting edge.
dbeg  = cumsum ( data ( 1: width - 2, : ) );
dbeg  = bsxfun ( @rdivide, dbeg ( 1: 2: end, : ), ( 1: 2: ( width - 2 ) )' );

% Calculates the ending edge.
dend  = cumsum ( data ( end: -1: end - width + 3, : ) );
dend  = bsxfun ( @rdivide, dend ( end: -2: 1, : ), ( width - 2: -2: 1 )' );

% Calculates the smoothed data.
cdata = cumsum ( data, 1 ) / width;
sdata = ( cdata ( width + 1: end, : ) - cdata ( 1: end - width, : ) );

% Combines the smoothed data with the calculated edges.
sdata = cat ( 1, dbeg, cdata ( width, : ), sdata, dend );
