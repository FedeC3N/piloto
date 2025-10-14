function [ trl, artifact ] = manualArtifacts ( data )
% PREPROCESADO DE DATOS MEG RESTING .FIF  
% version 10 agosto 2012


%   Los registros de entrada son registros de resting state de duraci�n
%   y filtrados con tsss-mc. En el caso de que el filtro sea SSS es
%   necesario verificar si se ha registrado HcPI en la senal raw. De ser
%   asi los datos han de haberse filtrado [1,150] antes de realizar el
%   filtro SSS. Despues ya es posible emplear este script.

%   Todo esto basado en la funcion de Naza 'Preproc_Resting.m', en la que
%   se ha cambiado la parte de rechazo de artefactos debido a
%   incompatibilidades entre las diferentes versiones de Fieltrip.

%   Los datos son segmentados en fragmentos de t_frag'' mediante la funcion
%   trialfun_MCI_rest

%   El proceso consiste en:

%   1-  Se coge el .fif original y se segmenta en fragmentos mediante la
%       funcion trialfun_MCI_rest.
%   2-  Se detecta la posicion en los datos de tres tipos de
%       artefactos: Jump, Muscle y EOG. Esto es personalizable mediante
%       unos umbrales para los valores de la transformada z. Luego los
%       trials que contengan artefactos son descartados.
%   3-  Se realiza un filtrado (pasa banda y band stop), un demean,
%       detrean y correcion por linea base. Todo es personalizable en
%       la seccion de parametros.
%   4-  Se montaa la estructura con los trials finales y se guarda.

%   Los par�metros de entrada son:

%   fif_path:       Ruta completa del fichero fif.
%   path_out:       Ruta completa de los datos preprocesados.
%   preprocesado    Estructura con todos los datos necesarios para el
%                   preprocesado. Se configura en el script RLV

% Checks the data options.
if ~isfield ( data, 'eog'         ), data.eog         = 'EEG61'; end % 'EEGO61'
if ~isfield ( data, 'channel'     ), data.channel     = 'MEG';   end % 'MEG' 'MEGMAG' 'MEGGRAD'
if ~isfield ( data, 'interactive' ), data.interactive = 'no';    end % 'yes' 'no'
if ~isfield ( data, 'artifacts'   ), data.artifacts   = [];      end

if ~isfield ( data.artifacts, 'eog'    ), data.artifacts.eog    = true; end
if ~isfield ( data.artifacts, 'jump'   ), data.artifacts.jump   = true; end
if ~isfield ( data.artifacts, 'muscle' ), data.artifacts.muscle = true; end

% Artifact detection parameters.

% Blink detection parameters.
eog.   remove      = data.artifacts.eog;
eog.   interactive = data.interactive;
eog.   channel     = data.eog;
eog.   cut_zvalue  = 5;        %   Valor por defecto 4
eog.   trlpadding  = 0.5;      %   Valor por defecto 0.5
eog.   fltpadding  = 0.25;     %   Valor por defecto 0.1  (Fragmento colaterales a los trials que se incluyen en el filtro)
eog.   artpadding  = 0.1;      %   Valor por defecto 0.1  (fragmento a cada lado del artefacto extra que elimina)

% Jump detection parameters.
jump.  remove      = data.artifacts.jump;
jump.  interactive = data.interactive;
jump.  cut_zvalue  = 20;       %   20 en tutorial
jump.  padding     = 0.1;      %   Valor 0 por defecto.
jump.  fltpadding  = 0.1;      %   (Fragmento colaterales a los trials que se incluyen en el filtro)

% Muscular activity detection parameters.
muscle.remove      = data.artifacts.muscle;
muscle.interactive = data.interactive;
muscle.cut_zvalue  = 6;        %   8--> Naza antes || 4--> tutorial
muscle.trlpadding  = 0.1;      %   Valor por defecto 0.1
muscle.fltpadding  = 0.1;      %   Valor por defecto 0.1
muscle.artpadding  = 0.1;      %   Valor por defecto 0.1

% Filtering parameters.
% opciones.bpf   = 'yes';        %   Se realiza un filtrado pasa banda a los trials ya limpios.
% opciones.fmin  = 1;            %   Frecuencia de corte inferior filtro pasa banda final.
% opciones.fmax  = 45;           %   Frecuencia de corte superior filtro pasa banda final.
% opciones.bsf   = 'no';         %   Se realiza un filtrado para eliminar componentes de la l�nea de tensi�n
% opciones.place = 'EU';         %   Zona donde se haya la MEG que ha generado el fif:
%                                %   ('EU' --> Europe (50 Hz))   Se filtra 50, 100 y 150 Hz.
%                                %   ('USA' --> USA (60 Hz))  Se filtra 60, 120 y 180 Hz.
% 
% opciones.demean  = 'yes';      %   Realiza un demean.
% opciones.detrend = 'yes';      %   Realiza un detrend.


% Trial definition for the artifact selection.
cfg = [];
cfg.principio    = data.start;          %   Tiempo inicial del periodo de Resting
cfg.final        = data.end;            %   Tiempo final del periodo de resting
cfg.trialfun     = data.trialfunartdet; %   Usamos la funcion personalizada
cfg.periodo      = data.segmentartdet;  %   Fija la duracion de cada trial
cfg.dataset      = data.fiffile;
cfg.headerformat = 'neuromag_fif';
cfg.dataformat   = 'neuromag_fif';
cfg.continuous   = 'yes';

% Gets consecutive trials.
% cfg = ft_definetrial(cfg);
cfg = ft_definetrial_nosaltasinohaytrials ( cfg );
trl = cfg.trl;

if isempty ( trl ), return, end


% ARTIFACT REJECTION

% EOG ARTIFACTS

cfg = [];
cfg.trl        = trl;
cfg.datafile   = data.fiffile;
cfg.headerfile = data.fiffile;
cfg.continuous = 'yes';

% Channel selection, cutoff and padding.
cfg.artfctdef.zvalue.channel       = eog.channel;
cfg.artfctdef.zvalue.cutoff        = eog.cut_zvalue;
cfg.artfctdef.zvalue.trlpadding    = eog.trlpadding;
cfg.artfctdef.zvalue.artpadding    = eog.artpadding;
cfg.artfctdef.zvalue.fltpadding    = eog.fltpadding;

% % Algorithmic parameters.
% cfg.artfctdef.zvalue.bpfilter      = 'yes';
% cfg.artfctdef.zvalue.bpfilttype    = 'but';
% cfg.artfctdef.zvalue.bpfreq        = [ 1 15 ];
% cfg.artfctdef.zvalue.bpfreq        = [ 0 15 ];
% cfg.artfctdef.zvalue.bpfiltord     = 4;
% cfg.artfctdef.zvalue.hilbert       = 'yes';

% Algorithmic parameters.
cfg.artfctdef.zvalue.bpfilter      = 'yes';
cfg.artfctdef.zvalue.bpfilttype    = 'fir';
cfg.artfctdef.zvalue.bpfreq        = [ 5 15 ];
cfg.artfctdef.zvalue.bpfiltord     = 400;
cfg.artfctdef.zvalue.hilbert       = 'yes';

% Interactive artifact rejection.
cfg.artfctdef.zvalue.interactive   = eog.interactive;

if eog.   remove
    fprintf ( 1, '  Searching for blinks.\n' );
    cfg             = ft_artifact_zvalue ( cfg );
    drawnow
    artifact.   eog = cfg.artfctdef.zvalue;
end

% JUMP ARTIFACTS
cfg = [];
cfg.trl        = trl;
cfg.padding    = jump.padding;
cfg.datafile   = data.fiffile;
cfg.headerfile = data.fiffile;
cfg.continuous = 'yes';

% Channel selection, cutoff and padding.
cfg.artfctdef.zvalue.channel    = data.channel;
cfg.artfctdef.zvalue.cutoff     = jump.cut_zvalue;
cfg.artfctdef.zvalue.trlpadding = cfg.padding / 2;
cfg.artfctdef.zvalue.artpadding = cfg.padding / 2;
cfg.artfctdef.zvalue.fltpadding = jump.fltpadding;

% Algorithmic parameters.
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';

% Interactive artifact rejection.
cfg.artfctdef.zvalue.interactive   = jump.interactive;

if jump.  remove
    fprintf ( 1, '  Searching for jumps.\n' );
    cfg             = ft_artifact_zvalue ( cfg );
    drawnow
    artifact.  jump = cfg.artfctdef.zvalue;
end

% MUSCLE ARTIFACTS

cfg = [];
cfg.trl        = trl;
cfg.datafile   = data.fiffile;
cfg.headerfile = data.fiffile;
cfg.continuous = 'yes';

% Channel selection, cutoff and padding.
cfg.artfctdef.zvalue.channel       = data.channel;
cfg.artfctdef.zvalue.cutoff        = muscle.cut_zvalue;
cfg.artfctdef.zvalue.trlpadding    = muscle.trlpadding;
cfg.artfctdef.zvalue.fltpadding    = muscle.fltpadding;
cfg.artfctdef.zvalue.artpadding    = muscle.artpadding;

% Algorithmic parameters.
cfg.artfctdef.zvalue.bpfilter      = 'yes';
cfg.artfctdef.zvalue.bpfreq        = [ 110 140 ];   %   Realiza un filtro y despues realiza una transformada H.
cfg.artfctdef.zvalue.bpfiltord     = 8;             %   Valor por defecto 8
cfg.artfctdef.zvalue.bpfilttype    = 'but';
cfg.artfctdef.zvalue.hilbert       = 'yes';
cfg.artfctdef.zvalue.boxcar        = 0.2;           %   Promedio de la amplitud de la transformada H.

% Interactive artifact rejection.
cfg.artfctdef.zvalue.interactive   = muscle.interactive;

if muscle.remove
    fprintf ( 1, '  Searching for muscular artifacts.\n' );
    cfg             = ft_artifact_zvalue ( cfg );
    drawnow
    artifact.muscle = cfg.artfctdef.zvalue;
end

% VISUAL ARTIFACT DETECTION


% Trial definition for the artifact selection.
cfg = [];
cfg.principio    = data.start;           %   Tiempo inicial del periodo de Resting
cfg.final        = data.end;             %   Tiempo final del periodo de resting
cfg.periodo      = data.segmentartview;  %   Fija la duracion de cada trial
cfg.dataset      = data.fiffile;
cfg.trialfun     = data.trialfunartview; %   Usamos la funcion personalizada
cfg.headerformat = 'neuromag_fif';
cfg.dataformat   = 'neuromag_fif';
cfg.continuous   = 'yes';

% Gets overlapping trials.
cfg = ft_definetrial_nosaltasinohaytrials ( cfg );
trl = cfg.trl;

% Expands the trials by 1 second.
trl ( :, 1 ) = trl ( :, 1 ) - 1000;
trl ( :, 2 ) = trl ( :, 2 ) + 1000;
trl ( :, 3 ) = -1000;

% Gets the MEG data for the defined trials.
cfg = [];
cfg.dataset     = data.fiffile;
cfg.trl         = trl;
cfg.channel     = 'MEG';
cfg.precission  = 'single';

trialdata       = ft_preprocessing ( cfg );

% Filters the data.
fir             = fir1 ( 1500, [ 2 45 ] / ( trialdata.fsample / 2 ) );
trialdata       = ft_myfiltfilt ( fir, 1, trialdata );

% Removes the padding.
cfg = [];
cfg.toilim      = [ 0 data.segmentartview - 1 / trialdata.fsample ];

trialdata       = ft_redefinetrial ( cfg, trialdata );

% Displays the MEG data to append or remove artifacts.
cfg = [];
cfg.channel     = 1: 30;
cfg.gradscale   = 0.05;
cfg.ylim        = [ -2 2 ] * 1e-12;
cfg.plotlabels  = 'some';
cfg.viewmode    = 'vertical';
cfg.continous   = 'no';
cfg.colorgroups = 'chantype';
cfg.artfctdef   = artifact;

cfg             = ft_databrowser ( cfg, trialdata );
drawnow
artifact        = cfg.artfctdef;

% Removes the segments with artifact.
% 
% cfg=[];
% cfg.artfctdef.reject          = 'complete';
% cfg.artfctdef.eog.artifact    = artifact.EOG;
% cfg.artfctdef.jump.artifact   = artifact.jump;
% cfg.artfctdef.muscle.artifact = artifact.muscle;
% 
% preprocdata = ft_rejectartifact(cfg,preprocdata);

% % Old artifact suppression method
% segments = cat ( 1, artifact.eog, artifact.jump, artifact.muscle );
% 
% cleantrials = true ( size ( trl, 1 ), 1 );
% for segment = 1: size ( segments, 1 )
%     
%     % Trials where the artifact starts or ends, and those that covers.
%     starts = segments ( segment, 1 ) > trl ( :, 1 ) & segments ( segment, 1 ) < trl ( :, 2 );
%     ends   = segments ( segment, 2 ) > trl ( :, 1 ) & segments ( segment, 2 ) < trl ( :, 2 );
%     covers = segments ( segment, 1 ) < trl ( :, 1 ) & segments ( segment, 2 ) > trl ( :, 2 );
%     
%     % Marks the trials where the artifact is present.
%     cleantrials ( starts | ends | covers ) = false;
% end

% Redefines the trials to avoid the artifacts.
trl = getRestingCleanTrials ( data, artifact );


% % Saves the preproc data.
% PreprocInf = [];
% 
% PreprocInf.trl      = trl;
% PreprocInf.subject  = data.subject;
% PreprocInf.dataset  = data.fif_path;
% PreprocInf.trlinfo  = data.name_out;
% PreprocInf.t_on     = data.t_on;
% PreprocInf.t_off    = data.t_off;
% PreprocInf.t_frag   = data.t_frag;
% PreprocInf.trialfun = data.trialfun;
% PreprocInf.channel  = data.channel;
% 
% PreprocInf.artifact.jump            = jump;
% PreprocInf.artifact.jump.artifact   = artifact.jump;
% 
% PreprocInf.artifact.muscle          = muscle;
% PreprocInf.artifact.muscle.artifact = artifact.muscle;
% 
% PreprocInf.artifact.eog             = eog;
% PreprocInf.artifact.eog.artifact    = artifact.eog;
% 
% save ( data.name_out, 'PreprocInf' );



function [cfg] = ft_definetrial_nosaltasinohaytrials(cfg)

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});

if ~isfield(cfg, 'trl') && (~isfield(cfg, 'trialfun') || isempty(cfg.trialfun))
    % there used to be other system specific trialfuns in previous versions
    % of fieldtrip, but they are deprecated and not included in recent
    % versions any more
    cfg.trialfun = 'trialfun_general';
    warning('no trialfun was specified, using trialfun_general');
end

% create the trial definition for this dataset and condition
if isfield(cfg, 'trl')
    % the trial definition is already part of the configuration
    fprintf('retaining exist trial definition\n');
    trl = cfg.trl;
    if isfield(cfg, 'event')
        fprintf('retaining exist event information\n');
        event = cfg.event;
    else
        event = [];
    end
    
elseif isfield(cfg, 'trialfun')
    
    % provide support for xxx and trialfun_xxx when the user specifies cfg.trialfun=xxx
    if exist(cfg.trialfun, 'file')
        % evaluate this function, this is the default
    elseif exist(['trialfun_' cfg.trialfun], 'file')
        % prepend trialfun to the function name
        cfg.trialfun = ['trialfun_' cfg.trialfun];
    else
        error('cannot locate the specified trialfun (%s)', cfg.trialfun)
    end
    
    % evaluate the user-defined function that gives back the trial definition
    fprintf('evaluating trialfunction ''%s''\n', cfg.trialfun);
    % determine the number of outpout arguments of the user-supplied trial function
    try
        % the nargout function in Matlab 6.5 and older does not work on function handles
        num = nargout(cfg.trialfun);
    catch %#ok<CTCH>
        num = 1;
    end
    if num==1
        % the user-defined function only gives back the trial definition
        trl   = feval(cfg.trialfun, cfg);
        event = [];
    else
        % the user-defined function also gives back detailed information about
        % conditions, reaction time or any other information
        [trl, event] = feval(cfg.trialfun, cfg);
    end
else
    error('no trialfunction specified, see FT_DEFINETRIAL for help');
end

if isfield(cfg, 'trialdef') && isfield(cfg.trialdef, 'eventtype') && isequal(cfg.trialdef.eventtype, '?')
    % give a gentle message instead of an error
    fprintf('no trials have been defined yet, see FT_DEFINETRIAL for further help\n');
elseif size(trl,1)<1
    trl=[]; %error('no trials were defined, see FT_DEFINETRIAL for help');
end

% add the new trials and events to the output configuration
fprintf('found %d events\n', length(event));
cfg.event = event;
fprintf('created %d trials\n', size(trl,1));
cfg.trl = trl;
