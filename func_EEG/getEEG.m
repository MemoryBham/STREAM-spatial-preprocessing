% Loads the data from the path specified in cfg, preprocesses based on cfg
% and saves everything to a FieldTrip struct.

function [dataEEG] = getEEG(cfg)

% set paths
filenames_data  = dir([cfg.path, '\*', cfg.subject,'*.bdf']);
cfg.dataset = [cfg.path,'\', filenames_data(1).name];

filenames_log  = dir([cfg.logpath, '\*', cfg.subject(1:2),'*.txt']);
cfg.logpath = [cfg.logpath '\', filenames_log(1).name];

% read log
log = getLog(cfg);
cfg.trlsPerSession = arrayfun(@(x) sum([log.values{:,2}]==x), unique([log.values{:,2}]));
cfg.ntrials = sum(cfg.trlsPerSession);

if isfield(cfg, 'lowmemory') && cfg.lowmemory
    
    for s = 1:cfg.numbSessions
        
        for ch = 1:length(cfg.channel)
            % load data and bandpass filter - using FieldTrip
            cfg_tmp = cfg;
            cfg_tmp = rmfield(cfg_tmp,'trialdef');
            cfg_tmp = rmfield(cfg_tmp,'resamplefs');
            cfg_tmp.channel = cfg.channel(ch)';
            data = ft_preprocessing(cfg_tmp);
            
            % Notch filter
            if isfield(cfg, 'notchfilter') && strcmp(cfg.notchfilter,'yes')
                data = NotchFilter(cfg, data);
            end
            
            % downsample
            data_ds = ft_resampledata(cfg, data);
            cfg.fsample = data_ds.fsample;
            cfg.hdr = data_ds.hdr;
            
            % check trigger for this session
            [cfg.trigger,cfg.ignored_trigger] = getEvents(cfg);
            
            
            % define trials
            cfg_tr = defineTrials_STREAMspatial(cfg,s);
            
            % cut trials
            cfg_tmp = [];
            cfg_tmp.trl = cfg_tr.trialdef.trl;
            data_tr = ft_redefinetrial(cfg_tmp,data_ds);
            
            % add trialinfo to struct
            data_tr.cue = cfg_tr.trialdef.cue;
            data_tr.stimulus = cfg_tr.trialdef.stimulus;
            data_tr.response = cfg_tr.trialdef.response;
            data_tr.event = cfg_tr.trialdef.event;
            data_tr.trialinfo = cfg_tr.trialdef.trialinfo;
            
            % merge data across sessions
            if ch == 1
                dataEEG_ch = data_tr;
            else
                dataEEG_ch = appendFTstruct(dataEEG_ch, data_tr, 'channel');
            end
        end
        % merge data across sessions
        if s == 1
            dataEEG = dataEEG_ch;
        else
            dataEEG = appendFTstruct(dataEEG, dataEEG_ch, 'trial');
        end
    end
else
    % load data and bandpass filter - using FieldTrip
    cfg_tmp = cfg;
    cfg_tmp = rmfield(cfg_tmp,'trialdef');
    cfg_tmp = rmfield(cfg_tmp,'resamplefs');
    data = ft_preprocessing(cfg_tmp);
    
    % Notch filter
    if isfield(cfg, 'notchfilter') && strcmp(cfg.notchfilter,'yes')
        data = NotchFilter(cfg, data);
    end
    
    % downsample
    data_ds = ft_resampledata(cfg, data);
    cfg.fsample = data_ds.fsample;
    cfg.hdr = data_ds.hdr;
    
    % check trigger for this session
    [cfg.trigger,cfg.ignored_trigger] = getEvents(cfg);
    
    for s = 1:cfg.numbSessions
        % define trials
        cfg_tr = defineTrials_STREAMspatial(cfg,s);
        
        % cut trials
        cfg_tmp = [];
        cfg_tmp.trl = cfg_tr.trialdef.trl;
        data_tr = ft_redefinetrial(cfg_tmp,data_ds);
        
        % add trialinfo to struct
        data_tr.cue = cfg_tr.trialdef.cue;
        data_tr.stimulus = cfg_tr.trialdef.stimulus;
        data_tr.response = cfg_tr.trialdef.response;
        data_tr.event = cfg_tr.trialdef.event;
        data_tr.trialinfo = cfg_tr.trialdef.trialinfo;
        
        % merge data across sessions
        if s == 1
            dataEEG = data_tr;
        else
            dataEEG = appendFTstruct(dataEEG, data_tr, 'trial');
        end
    end
end