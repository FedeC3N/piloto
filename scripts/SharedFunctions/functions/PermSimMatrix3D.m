function matrix_out = PermSimMatrix3D(matrix_in)

%   v1.0 14/09/2016

[nCols,~,nSubjs] = size(matrix_in);

%%   Convert 2 Vector
nData=nCols*(nCols-1)/2;
vector_aux=zeros(nData,nSubjs);

for isuj = 1:nSubjs
    cont=1;
    for nFila=2:nCols,
        for nCol=1:nFila-1,
            vector_aux(cont,isuj)=matrix_in(nFila,nCol,isuj);
            cont=cont+1;
        end
    end
end

%%   Permutations
for isuj = 1:nSubjs
    vector_aux(:,isuj) = vector_aux(randperm(nData),isuj);
end

%%   Convert to matrix
matrix_out = zeros(nCols,nSubjs);
for isuj = 1:nSubjs
    cont=1;
    for nFila=2:nCols,
        for nCol=1:nFila-1,
            matrix_out(nFila,nCol,isuj)=vector_aux(cont,isuj);
            matrix_out(nCol,nFila,isuj)=vector_aux(cont,isuj);
            cont=cont+1;
        end
    end
end

%%  Permutate the main diagonal

for isuj = 1:nSubjs
    aux = matrix_in(:,:,isuj);
    diag = aux(eye(nCols)==1);
    diag = diag(randperm(nCols));
    for icol = 1:nCols
        matrix_out(icol,icol,isuj) = diag(icol);
    end
end

end