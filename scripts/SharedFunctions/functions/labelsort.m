function sorted = labelsort ( template, labels )

% Checks that the teplate labels are unique.
if numel ( unique ( template ) ) ~= numel ( template )
    error ( 'The template has repeated labels.' );
end

% Reserves memory for the sorted labels matrix.
sorted = zeros ( size ( labels ) );

% Goes through each input label.
for label = 1: numel ( labels )
    
    % Finds the position of the label in the template.
    sorted ( label ) = find ( strcmp ( labels ( label ), template ) );
end