
function [conn_clus_matrix,clusters] = my_connclusters ( matrix )


% Gets the size of the data.
nnod       = size ( matrix, 1 );

% Makes sure that the matrix is binary.
if ~islogical ( matrix )
    posmatrix  = matrix > 0;
    negmatrix  = matrix < 0;
else
    posmatrix  = matrix;
    negmatrix  = zeros ( nnod );
end

% Makes sure that the matrix is simmetric and the diagonal is active.
posmatrix  = double ( posmatrix | posmatrix' | eye ( nnod ) > 0 );
negmatrix  = double ( negmatrix | negmatrix' | eye ( nnod ) > 0 );

% Closes the clusters.
conn_posmatrix  = posmatrix ^ ( nnod / 2 ) > 0;
conn_negmatrix  = negmatrix ^ ( nnod / 2 ) > 0;

% Removes the clusters with only one member.
conn_posmatrix  = conn_posmatrix ( sum ( conn_posmatrix, 2 ) > 2, : );
conn_negmatrix  = conn_negmatrix ( sum ( conn_negmatrix, 2 ) > 2, : );

% Gets the list of clusters.
poscluster = unique ( conn_posmatrix, 'rows' )';
negcluster = unique ( conn_negmatrix, 'rows' )';

% Sorts the clusters.
[ ~, ord ] = sort ( sum ( poscluster, 1 ), 'descend' );
poscluster = poscluster ( :, ord );
[ ~, ord ] = sort ( sum ( negcluster, 1 ), 'descend' );
negcluster = negcluster ( :, ord );

% Joins the positive and negative clusters.
clusters  = cat ( 2, poscluster, -negcluster );

% Create the positive connectivity structure
conn_clus_matrix_pos = zeros(size(posmatrix,1),size(posmatrix,1),size(poscluster,2));
for iclus = 1:size(poscluster,2)
    
    dummy_index = find(poscluster(:,iclus));
    conn_clus_matrix_pos(dummy_index,:,iclus) = posmatrix(dummy_index,:);
    
end
conn_clus_matrix_pos = conn_clus_matrix_pos | repmat(eye(size(posmatrix,1),size(posmatrix,1)),1,1,size(conn_clus_matrix_pos,3));

% Create the negative connectivity structure
conn_clus_matrix_neg = zeros(size(negmatrix,1),size(negmatrix,1),size(negcluster,2));
for iclus = 1:size(negcluster,2)
    
    dummy_index = find(negcluster(:,iclus));
    conn_clus_matrix_neg(dummy_index,:,iclus) = negmatrix(dummy_index,:);
    
end
conn_clus_matrix_neg = conn_clus_matrix_neg | repmat(eye(size(negmatrix,1),size(negmatrix,1)),1,1,size(conn_clus_matrix_neg,3));
conn_clus_matrix_neg = -1*conn_clus_matrix_neg;

conn_clus_matrix = cat(3, conn_clus_matrix_pos, conn_clus_matrix_neg);
conn_clus_matrix = double(conn_clus_matrix);