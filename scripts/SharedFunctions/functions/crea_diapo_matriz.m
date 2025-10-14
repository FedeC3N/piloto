function [dia_encod, dia_retri, iguales]  = crea_diapo_matriz(txt)
    
txt = txt{1,1};
regular = 'Encoding\w*\d*';
aux = regexp ( txt, regular);
indice = ~cellfun ( @isempty, aux );
ntrials = sum(indice);
lines = find(indice);

dia_encod = cell(ntrials,1);
dia_retri = cell(ntrials,1);

for itrial = 1:ntrials
    
    if (lines(itrial)+32>size(txt,1))
        aux = txt(lines(itrial):end);
    else
        aux = txt(lines(itrial):lines(itrial)+32);
    end
    
    % La siguiente línea a Encoding contiene la diapositiva.
    vect = regexp(aux,'Encoding: (\w*\d*)','tokens');
    indice = ~cellfun ( @isempty, vect );
    dia_encod{itrial} = vect{indice}{1};
    
    vect = regexp(aux,'Retrieval: (\w*\d*)','tokens');
    indice = ~cellfun ( @isempty, vect );  
    dia_retri{itrial} = vect{indice}{1};

end

iguales = cellfun(@strcmp,dia_encod,dia_retri);

end