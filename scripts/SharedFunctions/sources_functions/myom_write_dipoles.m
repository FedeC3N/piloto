function myom_write_dipoles ( basename, grid )

% If there is an 'inside' field in the grid, takes only those grid points.
if isfield ( grid, 'inside' )
    dipole.pos = grid.pos ( grid.inside, : );
else
    dipole.pos = grid.pos;
end

% If no original dipole orientation, uses the identity matrix.
if isfield ( grid, 'ori' )
    dipole.ori = grid.ori;
else
    dipole.ori = eye (3);
end

% Gets the number of dipoles and orientations.
ndipoles = size ( dipole.pos, 1 );
noris    = size ( dipole.ori, 1 );

% Matches the dipoles and their orientation.
if noris == 3 && det ( dipole.ori ) - 1 < 1e-6
    
    % The orientation must be applied to each dipole.
    dipole.pos = kron ( dipole.pos, ones ( 3, 1 ) );
    dipole.ori = kron ( ones ( ndipoles, 1 ), dipole.ori );
    
elseif noris == ndipoles
    
    % Each dipole has its orientation. Nothing to do.
    
elseif noris == ndipoles * 3
    
    % Each dipole has 3 orientations.
    dipole.pos = kron ( dipole.pos, ones ( 3, 1 ) );
    
else
    error ( 'Imposible to match dipoles and orientations.' );
end

% Triplicates the dipole positions, one repetition for each orientation.
dipoles  = cat ( 2, dipole.pos, dipole.ori );

% Writes the dipole file.
om_save_full ( dipoles, sprintf ( '%s_dip.bin', basename ) );
