function cleandata = kalima ( data, EKG )

% Filters the EKG data and gets its envelope.
fir   = fir1 ( 1000, 10 / ( 1000 / 2 ), 'high' );
EKG   = myfiltfilt ( fir, 1, EKG, true );
EKG   = abs ( EKG );

% Gets the areas where the EKG is over the threshold.
thres = std ( EKG, 0, 1 );
thres = repmat ( thres, [ size( EKG, 1 ) 1 1 ] );
trigs = EKG (:) > 2 * thres (:);
trigs = diff ( trigs );
ups   = find ( trigs > 0 );
downs = find ( trigs < 0 );

% Removes the edges.
downs ( downs < ups   (1)     ) = [];
ups   ( ups   > downs ( end ) ) = [];

% Gets the maximum in each area.
trigs = [ ups downs ];
trigs ( diff ( trigs, [], 2 ) < 20, : ) = [];
peaks = zeros ( size ( trigs, 1 ), 1 );

for trig = 1: size ( trigs, 1 )
    segment = EKG ( trigs ( trig, 1 ): trigs ( trig, 2 ) );
    peak    = max ( segment );
    peaks   ( trig ) = trigs ( trig, 1 ) + find ( segment == peak, 1 );
end

% Creates the triggers.
trigs = false ( size ( EKG ) );
trigs ( peaks ) = true;

% Removes triggers in the windows' edges.
trigs ( 1: 500,         :, : ) = 0;
trigs ( end - 500: end, :, : ) = 0;

% Analizes the problem case by case.
while find ( diff ( find ( trigs ) ) < 500, 1 )
    
    % Gets the peaks in conflict.
    peaks = find ( trigs );
    peak1 = peaks ( find ( diff ( peaks ) < 500, 1 ) + 0 );
    peak2 = peaks ( find ( diff ( peaks ) < 500, 1 ) + 1 );
    
    if EKG ( peak1 ) > EKG ( peak2 ), trigs ( peak2 ) = false; end
    if EKG ( peak1 ) < EKG ( peak2 ), trigs ( peak1 ) = false; end
end

% Gets the position of the R wave.
rpos  = find ( trigs );

% Sets the QRS lenght as the minimum distance between beats.
qrsl  = min ( diff ( rpos ) );


% Rewrites the MEG data as a continous signal.
cdata = permute ( data, [ 2 1 3 ] );
cdata = cdata ( :, : );


% Reserves memory for the QRS waveforms.
qrss  = zeros ( size ( cdata, 1 ), qrsl, numel ( rpos ) );

% Extracts the QRS waveforms window from the EKG components.
for c = 1: numel ( rpos )
    
    % Extracts the QRS complex edges.
    qrsb = rpos ( c ) - round ( qrsl / 3 ) + 1;
    qrse = rpos ( c ) + round ( 2 * qrsl / 3 );
    
    qrss ( :, :, c ) = cdata ( :, qrsb: qrse );
end

% Calculates the average QRS waveform.
qrs   = mean ( qrss, 3 );


% Constructs a signal of QRS complexes: The EKG projection.
proj  = zeros ( size ( cdata ) );

for c = 1: numel ( rpos )
    
    % Extracts the QRS complex edges.
    qrsb = rpos ( c ) - round ( qrsl / 3 ) + 1;
    qrse = rpos ( c ) + round ( 2 * qrsl / 3 );
    
    % Goes through each input signal.
    for signal = 1: size ( cdata, 1 )
        
        % Scales the QRS complex to maximize the projection.
        scale  = regress ( detrend ( qrss ( signal, :, c )', 0 ), detrend ( qrs ( signal, : )', 0 ) );
        
        proj   ( signal, qrsb: qrse ) = qrs ( signal, : ) * scale;
    end
end

cleandata = cdata - proj;
cleandata = reshape ( cleandata, size ( permute ( data, [ 2 1 3 ] ) ) );
cleandata = permute ( cleandata, [ 2 1 3 ] );

% if size ( data, 1 ) == 1
%     
%     fir = fir1 ( 1000, .5 / 500, 'high' );
%     
%     pdata      = permute ( myfiltfilt ( fir, 1, permute ( data,      [ 2 1 3 ] ) ), [ 2 1 3 ] );
%     pEKGor     = permute ( myfiltfilt ( fir, 1, permute ( EKGor,     [ 2 1 3 ] ) ), [ 2 1 3 ] );
%     pcleandata = permute ( myfiltfilt ( fir, 1, permute ( cleandata, [ 2 1 3 ] ) ), [ 2 1 3 ] );
%     
%     pdata      = zscore ( pdata      ( :, 1501: end - 1500, : ), [], 2 );
%     pEKGor     = zscore ( pEKGor     ( :, 1501: end - 1500, : ), [], 2 );
%     pcleandata = zscore ( pcleandata ( :, 1501: end - 1500, : ), [], 2 );
%     
%     figure
%     hold on
%     plot ( pEKGor     ( :, : )', 'g' )
%     plot ( pdata      ( :, : )', 'r' )
%     plot ( pcleandata ( :, : )', 'b' )
%     
%     figure
%     hold on
%     plot ( linspace ( 0, 1000, 4000 ), mean ( abs ( fft ( pEKGor,     [], 2 ) ), 3 ), 'g' )
%     plot ( linspace ( 0, 1000, 4000 ), mean ( abs ( fft ( pdata,      [], 2 ) ), 3 ), 'r' )
%     plot ( linspace ( 0, 1000, 4000 ), mean ( abs ( fft ( pcleandata, [], 2 ) ), 3 ), 'b' )
%     xlim ( [ 0 45 ] )
% end
