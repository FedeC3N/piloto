function dataout = read_ePrime_log_nii(txt, tarea_fif)

%   v1.1 12/10/2016
%%  txt es la cell sacada tras leer con
%   txt = textread('NombreLog.txt','%q',50000);


%   Palabra que identificará el principio de cada trial.
txt = txt{1,1};
regular = 'Procedure: Trial\w*\d*';
aux = regexp ( txt, regular);

indice = ~cellfun ( @isempty, aux );
ntrials = sum(indice);
lines = find(indice);

clear indice aux regular

%   En este caso recogemos de cada trial los siguientes campos:
%   TriggerEncodi, TriggerInterruptio, TriggerRecu
%   Recu.OnsetTime, Recu.RTTime, Recu.RESP
%   Respuesta.RTTime, Respuesta.RESP. En total 8 campos:
tags = {'TriggerEncoding','TriggerInterruption','TriggerRetrieval','ImagenRetrievalOnsetTime','RespuestaImagenRetrievalRTTime','RespuestaImagenRetrieval','RespuestaRTTime','Respuesta','Bloque', 'RespuestaInterrupcion'};

info = zeros(ntrials,10);

% Los logs de interrupcion ocupan más. Cogemos unas 50 líneas.
tamano = 40;
if strcmp(tarea_fif, 'int') tamano = 51; end

for itrial = 1:ntrials
    %   Se extrae el bloque entero
    if (lines(itrial)+tamano>size(txt,1))
        aux = txt(lines(itrial):end);
    else
        aux = txt(lines(itrial):lines(itrial)+tamano);
    end
    
    %   TriggerEncodi
    %     info(itrial,1) = str2double(aux{4}); Lo he modificado
    encoding = regexp ( aux, 'TriggerEncodi: (\d*)','tokens');
    indice = ~cellfun ( @isempty, encoding );
    if numel(indice)>1 indice = find(indice,1,'first'); end
    info(itrial,1) = str2double(encoding{indice}{1});
    clear encoding indice
    %   TriggerInterruptio
    %     info(itrial,2) = str2double(aux{10});
    encoding = regexp ( aux, 'TriggerInterruptio: (\d*)','tokens');
    indice = ~cellfun ( @isempty, encoding );
    if numel(indice)>1 indice = find(indice,1,'first'); end
    info(itrial,2) = str2double(encoding{indice}{1});
    clear encoding indice
    %   TriggerRecu
    %     info(itrial,3) = str2double(aux{12});
    encoding = regexp ( aux, 'TriggerRecu: (\d*)','tokens');
    indice = ~cellfun ( @isempty, encoding );
    if numel(indice)>1 indice = find(indice,1,'first'); end
    info(itrial,3) = str2double(encoding{indice}{1});
    clear encoding indice
    
    %   Recu.OnsetTime
    vect = regexp(aux,'Recu[1-4]?.OnsetTime: (\d*)','tokens');
    indice = ~cellfun ( @isempty, vect );
    if numel(indice)>1 indice = find(indice,1,'first'); end
    info(itrial,4) = str2double(vect{indice}{1});
    if isempty(vect{indice}{1}{1})
        info(itrial,4) = 1;
    else
        info(itrial,4) = str2double(vect{indice}{1}{1});
    end
    clear vect indice
    
    %   Recu.RT
    vect = regexp(aux,'Recu[1-4]?.RT: (\d*)','tokens');
    indice = find(~cellfun(@isempty,vect));
    if numel(indice)>1 indice = find(indice,1,'first'); end
    if isempty(indice)
        info(itrial,5) = nan;
    else
        info(itrial,5) = str2double(vect{indice}{1}{1})/1000;
    end
    clear vect indice
    
    %   Recu.RESP
    vect = regexp(aux,'Recu[1-4]?.RESP: (\d*)', 'tokens');
    indice = find(~cellfun(@isempty,vect));
    if numel(indice)>1 indice = find(indice,1,'first'); end
    if isempty(indice)
        RecuRESP = nan;
    else
        RecuRESP = str2double(vect{indice}{1}{1});
    end
    if isnan(RecuRESP) %   Cuando no hay respuesta se salta una fila
        info(itrial,6) = 0;
    else
        info(itrial,6) = RecuRESP;
    end
    clear vect indice RecuRESP
    
    %   Respuesta.RTTime
    vect = regexp(aux,'Respuesta[1-4]?.RT: (\d*)','tokens');
    indice = find(~cellfun(@isempty,vect));
    if numel(indice)>1 indice = find(indice,1,'first'); end
    if isempty(indice)
        info(itrial,7) = nan;
    else
        info(itrial,7) = (str2double(vect{indice}{1}{1})+1000)/1000;
    end
    clear vect indice
    
    %   Respuesta.RESP
    vect = regexp(aux,'Respuesta[1-4]?.RESP: (\d*)', 'tokens');
    indice = find(~cellfun(@isempty,vect));
    if numel(indice)>1 indice = find(indice,1,'first'); end
    if isempty(indice)
        info(itrial,8) = nan;
    else
        info(itrial,8) = str2double(vect{indice}{1}{1});
    end
    clear vect index
    
    % Bloque del trial
    vect = regexp(aux,'Procedure: TrialProc(\d*)','tokens');
    indice = find(~cellfun(@isempty,vect));
    if numel(indice)>1 indice = find(indice,1,'first'); end
    if isempty(indice)
        info(itrial,9) = nan;
    else
        info(itrial,9) = str2double(vect{indice}{1}{1});
        if info(itrial,9) == 0 || info(itrial,9) == 1
            info(itrial,9) = info(itrial,9) + 1;
        end
    end
    clear vect indice
    
    vect = regexp(aux,'Interruptio[1-4]?.RESP: (\d*)','tokens');
    indice = find(~cellfun(@isempty,vect));
    if numel(indice)>1 indice = find(indice,1,'first'); end
    if isempty(indice)
        info(itrial,10) = nan;
    else
        info(itrial,10) = str2double(vect{indice}{1}{1});
    end
    clear vect indice  aux
    
end
dataout.info = info;
dataout.tags = tags;

end