function [dataOut,keepID] = selectData_STREAMspatial(cfg,dataIn)
% selectData_STREAM
% 	phase: 1 = encoding; 2 = retrieval; 3 = catch
% 	perceptual: 1 = line; 2 = photo
% 	semantic: 1 = animate; 2 = inanimate
% 	ncorrect: 0, 1 or 2 correct catch questions

% set some defaults
if ~isfield(cfg,'checkInt'); cfg.checkInt = false; end

% find datatype
if isfield(dataIn, 'powspctrm')
    datatype = 'spectrum';
elseif isfield(dataIn, 'timestamp')
    datatype = 'spikes';
    nunits = length(dataIn.label);
else
    datatype = 'LFPs';
end


%% select trials

% eliminate trials to keep based on cfg
keepID = true(1,size(dataIn.trialinfo,2)); % keep all trials

% only use trials specified in config
if isfield(cfg,'trial')
    keepID(setdiff(find(keepID),cfg.trials)) = false;
end

% only use specified session
if isfield(cfg,'session')
    keepID(~ismember([dataIn.trialinfo{1,:}],cfg.session)) = false;
end

% phase: 1 = familiarization; 2 = encoding; 3 = retrieval
if isfield(cfg,'phase')
    if sum(cfg.phase == 1) || sum(strcmp(cfg.phase,{'fam','familiarization'})) % keep familiarization
        cfg.phase = 1;
        keepID(strcmp(dataIn.trialinfo(2,:), 'encoding') | strcmp(dataIn.trialinfo(2,:), 'retrieval')) = false;
    elseif sum(cfg.phase == 2) || sum(strcmp(cfg.phase,{'enc','encoding'})) % keep retrieval or catch
        cfg.phase = 2;
        keepID(strcmp(dataIn.trialinfo(2,:), 'familiarization') | strcmp(dataIn.trialinfo(2,:), 'retrieval')) = false;
    elseif sum(cfg.phase == 3) || sum(strcmp(cfg.phase,{'ret','retrieval'})) % keep retrieval or catch
        cfg.phase = 3;
        keepID(strcmp(dataIn.trialinfo(2,:), 'familiarization') | strcmp(dataIn.trialinfo(2,:), 'encoding')) = false;
    else
        error('Unknown phase selected')
    end
end

% category: 1 = dog; 2 = bird; 3 = car; 4 = plane;
if isfield(cfg,'category')
    if sum(cfg.category == 1) || strcmp(cfg.category,'dogs') % keep dogs
        keepID([dataIn.stimulus{4,:}] ~= 1) = false;
    elseif sum(cfg.category == 2) || strcmp(cfg.category,'birds') % keep birds
        keepID([dataIn.stimulus{4,:}] ~= 2) = false;
    elseif sum(cfg.category == 3) || strcmp(cfg.category,'cars') % keep car
        keepID([dataIn.stimulus{4,:}] ~= 3) = false;
    elseif sum(cfg.category == 4) || strcmp(cfg.category,'planes') % keep planes
        keepID([dataIn.stimulus{4,:}] ~= 4) = false;
    else
        error('Unknown semantic category selected')
    end
end

% perceptual: 1 = line; 2 = photo
if isfield(cfg,'perc1')
    if cfg.perc1 == 1 % keep line
        keepID([dataIn.stimulus{5,:}] == 2) = false;
    elseif cfg.perc1 == 2 % keep photo
        keepID([dataIn.stimulus{5,:}] == 1) = false;
    else
        error('Unknown perceptual category selected')
    end
end
% perceptual: 1 = colour; 2 = black/white
if isfield(cfg,'perc2')
    if cfg.perc2 == 1 % keep colour
        keepID([dataIn.stimulus{6,:}] == 2) = false;
    elseif cfg.perc2 == 2 % keep black/white
        keepID([dataIn.stimulus{6,:}] == 1) = false;
    else
        error('Unknown perceptual category selected')
    end
end
% semantic: 1 = living; 2 = non-living
if isfield(cfg,'sem1')
    if cfg.sem1 == 1 % keep living
        keepID([dataIn.stimulus{7,:}] == 2) = false;
    elseif cfg.sem1 == 2 % keep non-living
        keepID([dataIn.stimulus{7,:}] == 1) = false;
    else
        error('Unknown semantic category selected')
    end
end
% semantic: 1 = flying; 2 = non-flying
if isfield(cfg,'sem2')
    if cfg.sem2 == 1 % keep living
        keepID([dataIn.stimulus{8,:}] == 2) = false;
    elseif cfg.sem2 == 2 % keep non-living
        keepID([dataIn.stimulus{8,:}] == 1) = false;
    else
        error('Unknown semantic category selected')
    end
end

% select specific trial types within a phase
% for example:
% familiarization: select specific catch questions
% encoding: select DaD trials or normal trials
% retrieval: select non-catch questions of a specific catch question
if isfield(cfg,'phase') && isfield(cfg, 'trialtype')
    keepID(~ismember([dataIn.trialinfo{3,:}],cfg.trialtype)) = false;
end

% keep remembered or forgotten trials only (subjective button press)
if isfield(cfg, 'remembered')
    if isfield(cfg, 'phase') && cfg.phase == 3
        keepID([dataIn.response{2,:}] ~= cfg.remembered) = false;
    else
        error('Cannot select remembered trials for non-retrieval trials')
    end
end

% keep trials with correct of incorrect catch questions
if isfield(cfg, 'correct')
    if isfield(cfg, 'phase') && ismember(cfg.phase, [1,3])
        keepID([dataIn.response{4,:}] ~= cfg.correct) = false;
    else
        error('Cannot select correct/incorrect for encoding trials')
    end
end

% remove RTs that are too short or too long
if isfield(cfg, 'minRT')
    if isfield(cfg,'phase') && ismember(cfg.phase, [2,3])
        keepID([dataIn.response{1,:}] < cfg.minRT) = false;
    else
        error('Cannot use minRT for familiarization')
    end
end
if isfield(cfg, 'maxRT')
    if isfield(cfg,'phase') && ismember(cfg.phase, [2,3])
        keepID([dataIn.response{1,:}] > cfg.maxRT) = false;
    else
        error('Cannot use maxRT for familiarization')
    end
end 

if isfield(cfg, 'minCatchRT')
    if isfield(cfg,'phase') && ismember(cfg.phase, [1,3])
        keepID([dataIn.response{5,:}] < cfg.minCatchRT) = false;
    else
        error('Cannot use minCatchRT for encoding')
    end
end
if isfield(cfg, 'maxCatchRT')
    if isfield(cfg,'phase') && ismember(cfg.phase, [1,3])
        keepID([dataIn.response{5,:}] > cfg.maxCatchRT) = false;
    else
        error('Cannot use maxCatchRT for encoding')
    end
end

% create output struct
dataOut = dataIn;
if isfield(dataIn, 'event')
    dataOut.event = dataIn.event(:,keepID);
    dataOut.cue = dataIn.cue(:,keepID);
    dataOut.stimulus = dataIn.stimulus(:,keepID);
    dataOut.response = dataOut.response(:,keepID);
end
if strcmp(datatype, 'spectrum')
    dataOut.powspctrm = dataIn.powspctrm(keepID);
    dataOut.time = dataIn.time(keepID);
elseif strcmp(datatype,'LFPs')
    dataOut.trial = dataIn.trial(keepID);
    dataOut.time = dataIn.time(keepID);
elseif strcmp(datatype, 'spikes')
    trialOrig = 1:length(dataIn.trialinfo);
    trialNew = zeros(size(trialOrig));
    trialNew(keepID) = 1:sum(keepID);
    for u = 1:nunits
        keepSpikes = (ismember(dataIn.trial{u},trialOrig(keepID)));
        dataOut.time{u} = [dataIn.time{u}(keepSpikes)];
        dataOut.trial{u} = trialNew(dataIn.trial{u}(keepSpikes));
        dataOut.timestamp{u} = dataIn.timestamp{u}(keepSpikes);
%         dataOut.waveform{u} = dataIn.waveform{u}(:,:,keepSpikes);
    end
    dataOut.trialtime = dataIn.trialtime(keepID,:);
end
dataOut.trialinfo = dataIn.trialinfo(:,keepID);
if isfield(dataIn, 'sampleinfo')
    dataOut.sampleinfo = dataIn.sampleinfo(keepID,:);
end
if isfield(dataIn.cfg, 'trl')
    dataOut.cfg.trl = dataIn.cfg.trl(keepID,:);
end

%% select events within trial

if strcmp(datatype,'spikes')
    keepIDout = true(1,size(dataOut.trialtime,1));
else
    keepIDout = true(1,length(dataOut.time));
end

if isfield(cfg, 'toi') && ~isfield(cfg, 'event') || ...
        ~isfield(cfg, 'toi') && isfield(cfg, 'event')
    warning('Event and toi need to both be specified to select data')
elseif isfield(cfg, 'toi') && isfield(cfg, 'event')
    
    % find time and index step
    if length(cfg.toi) == 2 && ~isfield(cfg, 'dt')
        cfg.dt = 1/dataOut.fsample; %dataOut.time{1}(1,2) - dataOut.time{1}(1,1);
        cfg.toi = cfg.toi(1):cfg.dt:cfg.toi(2);
    elseif length(cfg.toi) > 2
        cfg.dt = cfg.toi(2)-cfg.toi(1);
    end
    DT = round(cfg.dt*dataOut.fsample); % index stepsize
    
    if strcmp(datatype, 'spikes')
        
        for u = 1:nunits
            if isempty(dataOut.time{u})
                dataOut.time{u} = [];
                dataOut.trial{u} = [];
                dataOut.timestamp{u} = [];
                dataOut.waveform{u} = [];
                continue
            end
            keepSpikes = dataOut.time{u} > ([dataOut.event{cfg.event,dataOut.trial{u}}] + cfg.toi(1)) & ...
                dataOut.time{u} < ([dataOut.event{cfg.event,dataOut.trial{u}}] + cfg.toi(end));
            
            if isempty(keepSpikes) || sum(keepSpikes) == 0
                dataOut.time{u} = [];
                dataOut.trial{u} = [];
                dataOut.timestamp{u} = [];
                dataOut.waveform{u} = [];
                continue
            end
            dataOut.time{u} = dataOut.time{u}(keepSpikes) - [dataOut.event{cfg.event,dataOut.trial{u}(keepSpikes)}];
            dataOut.trial{u} = dataOut.trial{u}(keepSpikes);
            dataOut.timestamp{u} = dataOut.timestamp{u}(keepSpikes);
            dataOut.waveform{u} = dataOut.waveform{u}(:,:,keepSpikes);
        end
        
        dataOut.trialtime = repmat([cfg.toi(1),cfg.toi(end)],[length(keepIDout),1]);
    else
        
        % select relevant data for each trial
        for tr = 1:length(dataOut.time)
            
            % find event ID
            if cfg.event == 0
                eventID = 1;
            else
                [~,eventID] = min(abs(dataOut.time{tr}-dataOut.event{cfg.event,tr}));
            end
            
            % select data
            selectIDs = eventID+round(cfg.toi(1)*dataOut.fsample): DT: eventID+round(cfg.toi(end)*dataOut.fsample);
            
            if sum(selectIDs < 1)
                error('Window starts too early')
            elseif selectIDs(end) > length(dataOut.time{tr})
                error('Window ends too late')
            end
            
            if isfield(dataIn, 'trial')
                dataOut.trial{tr} = dataOut.trial{tr}(:,selectIDs);
            elseif isfield(dataIn, 'powspctrm')
                dataOut.powspctrm{tr} = dataOut.powspctrm{tr}(:,:,selectIDs);
            end
            
            % check trial for NaNs at beginning or end (to allow for interpolation)
            if cfg.checkInt && isfield(dataIn, 'trial')
                if isnan(sum(dataOut.trial{tr}(:,1),1)) || isnan(sum(dataOut.trial{tr}(:,end),1))
                    % delete trial
                    keepIDout(tr) = false;
                    % also delete from original keepID
                    trIDs = find(keepID==1,tr);
                    keepID(trIDs(end)) = false;
                end
            elseif cfg.checkInt && isfield(dataIn, 'powspctrm')
                warning('Cannot check for NaNs in spectra. Setting ignored.')
            end
            
            % set time axis
            dataOut.time{tr} = cfg.toi;
        end
        if isfield(dataIn, 'event')
            dataOut.event = num2cell(zeros(1,length(dataOut.event(1,:)))); % event is now set to 0
        end
        
        % create output struct
        if isfield(dataIn, 'event')
            dataOut.event = dataOut.event(:,keepIDout);
            dataOut.cue = dataOut.cue(:,keepIDout);
            dataOut.stimulus = dataOut.stimulus(:,keepIDout);
            dataOut.response = dataOut.response(:,keepIDout);
        end
        if isfield(dataIn, 'powspctrm')
            dataOut.powspctrm = dataOut.powspctrm(keepIDout);
        end
        if isfield(dataIn, 'trial')
            dataOut.trial = dataOut.trial(keepIDout);
        end
        dataOut.time = dataOut.time(keepIDout);
        dataOut.trialinfo = dataOut.trialinfo(:,keepIDout);
        dataOut.sampleinfo = dataOut.sampleinfo(keepIDout,:);
    end
end

end
