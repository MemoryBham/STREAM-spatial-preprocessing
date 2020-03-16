
%%%%%%%%%%%% defineTrials_STREAMspatail
% DESCRIPTION:
% This functions identifies the trials for the STREAM paradigm based on the trigger
% and outputs this to a variable trl which can be used by the FieldTrip
% function ft_definetrial.

% INPUT:
% - config specifying:
% .dataset and .eventfile: paths or pre-analysed trigger (output of getEvents);
% trialdef.eventvalue: trigger value to look for
% trialdef.prestim and trialdef.poststim: start and end of trial relative
% to trigger. If left empty, poststim will be taken to be the trigger of
% the next trial.
% - FT datastruct, this can either contain LFP data or spikes
% - session: if multiple sessions were used

% OUTPUT
% data: FieldTrip struct with trial cut data & trial information
% The FieldTrip struct has the following fields:
% .hdr: struct with original recording info
% .label: cell with channel labels
% .time: cell with [1 x ntimes] time stamps for all trials
% .trial: cell with [nchan x ntimes] recordings
% .trialinfo: cell with block number, trial type (encoding/retrieval), switch trial
% .fsample: int. sample frequency
% .cfg
% .event: cell with event info (time stamp, event type) for every trial, in sec relative to trial onset
% .cue: cell with verb cue for every trial
% .stimulus: cell with object for encoding trials
% .response: cell with question info and correct/incorrect for every trial


function cfg = defineTrials_STREAMspatial(cfg, session)

%% check inputs and data type

if isempty(cfg);        error('Config invalid');   end
if isempty(session);    session = 1;               end

%% load header, log and trigger

log          = getLog (cfg);

if isfield(cfg, 'trigger')
    trigger_nods = cfg.trigger;
else
    [trigger_nods,~] = getEvents(cfg,session);
end

%% save some trial info

% identify relevant columns in the log
fieldID_trialid = find(strcmp(log.hdr, 'trial_id'));

fieldID_phase = find(strcmp(log.hdr, 'block_state'));
fieldID_session = find(strcmp(log.hdr, 'sessionID'));
fieldID_catch = find(strcmp(log.hdr, 'catch_type'));
fieldID_trial = find(strcmp(log.hdr, 'trial_type'));

fieldID_cue = find(strcmp(log.hdr, 'cue_id'));
fieldID_cuex = find(strcmp(log.hdr, 'cue_xcoord'));
fieldID_cuey = find(strcmp(log.hdr, 'cue_ycoord'));

fieldID_stim = find(strcmp(log.hdr, 'stim_id'));
fieldID_stimlabel = find(strcmp(log.hdr, 'stim_label'));
fieldID_stimcount = find(strcmp(log.hdr, 'stim_counter'));
fieldID_stimcat = find(strcmp(log.hdr, 'stim_cat'));
fieldID_stimperc1 = find(strcmp(log.hdr, 'stim_perc1'));
fieldID_stimperc2 = find(strcmp(log.hdr, 'stim_perc2'));
fieldID_stimsem1 = find(strcmp(log.hdr, 'stim_sem1'));
fieldID_stimsem2 = find(strcmp(log.hdr, 'stim_sem2'));

fieldID_RTenc = find(strcmp(log.hdr, 'RT_encoding'));
fieldID_RTret = find(strcmp(log.hdr, 'RT_ret_reinst'));
fieldID_retans = find(strcmp(log.hdr, 'ret_reinst_resp'));
fieldID_RTcatch = find(strcmp(log.hdr, 'RT_catch'));
fieldID_acc = find(strcmp(log.hdr, 'catch_resp'));

trialinfo           = {};
% .trialinfo: cell with block number, trial type (encoding/retrieval),
% switch trial and a unique identifier (to match encoding to retrieval)
sessionID = find(ismember([log.values{:,fieldID_session}], session));
id1 = sessionID(1);
id2 = sessionID(end);
ntrials = length(sessionID);

trialinfo(1,:)      = log.values(id1:id2,fieldID_session)'; % block number
trialinfo(2,:)      = log.values(id1:id2,fieldID_phase)'; % phase (familiarization/encoding/retrieval)
tmp = [log.values{id1:id2,fieldID_trial}];
tmp2 = mod(tmp,5);
tmp2(tmp2==0 & tmp>0) = 5;
trialinfo(3,:)      = num2cell(tmp2); % trial type (DaD / catch question)
trialinfo(4,:)      = log.values(id1:id2,fieldID_catch)'; % catch type

cue(1,:) = num2cell((session-1)*8+[log.values{id1:id2,fieldID_cue}]'); % cue id
cue(2,:) = log.values(id1:id2,fieldID_cuex)'; % cue x 
cue(3,:) = log.values(id1:id2,fieldID_cuey)'; % cue x 

stimulus(1,:) = log.values(id1:id2,fieldID_stim)'; % object id
stimulus(2,:) = log.values(id1:id2,fieldID_stimlabel)'; % object name
stimulus(3,:) = log.values(id1:id2,fieldID_stimcount)'; % stimulus count
stimulus(4,:) = log.values(id1:id2,fieldID_stimcat)'; % stimulus cat
stimulus(5,:) = log.values(id1:id2,fieldID_stimperc1)'; % object's perceptual info: line = 1; photo = 2;
stimulus(6,:) = log.values(id1:id2,fieldID_stimperc2)'; % object's perceptual info: color = 1; black/white = 2;
stimulus(7,:) = log.values(id1:id2,fieldID_stimsem1)'; % object's semantic info: living = 1; non-living = 2;
stimulus(8,:) = log.values(id1:id2,fieldID_stimsem2)'; % object's semantic info: flying = 1; non-flying = 2;

dum = nansum([[log.values{id1:id2,fieldID_RTenc}];[log.values{id1:id2,fieldID_RTret}]],1);
dum(1,dum(1,:)==0) = NaN;
response(1,:) = num2cell(dum); % RTs of mental image
response(2,:) = log.values(id1:id2,fieldID_retans)'; % remembered or not?
response(3,:) = log.values(id1:id2,fieldID_catch)'; % catch I type
response(4,:) = log.values(id1:id2,fieldID_acc)'; % catch I response
response(5,:) = log.values(id1:id2,fieldID_RTcatch)'; % catch I RT

%% define trial times, cut trials and store details

event = cell(5,ntrials);
trl = zeros(ntrials,3);

% find time stamps of trial start and end
% trl_sid = start of trial - sample ID after downsampling
% trl_eid = end of trial - sample ID after downsampling

% trl: Ntrials x 3 array:
% 1st column: sample index of trial start in raw data
% 2nd column: sample index of trial end in raw data
% 3rd column: trial start time (in s) relative to trigger 

fieldID_trigonset = find(strcmp(log.hdr, 'onset_trigger'));
fieldID_trialonset = find(strcmp(log.hdr, 'onset_trial'));
fieldID_catchresponset = find(strcmp(log.hdr, 'onset_catch_resp'));

fieldID_cueonset = find(strcmp(log.hdr, 'onset_cue'));
fieldID_stimonset = find(strcmp(log.hdr, 'onset_stim'));
fieldID_catchonset = find(strcmp(log.hdr, 'onset_catch'));

fieldID_dadstimonset = find(strcmp(log.hdr, 'onset_DaD_stim'));
fieldID_dadloconset = find(strcmp(log.hdr, 'onset_DaD_loc'));
fieldID_dadRTloc = find(strcmp(log.hdr, 'RT_DaD_loc'));
fieldID_dadRTreinst = find(strcmp(log.hdr, 'RT_DaD_reinst'));

fieldID_retreinst = find(strcmp(log.hdr, 'onset_ret_reinst'));

for ntr = 1:ntrials
    tr = id1 - 1 + ntr;
    
    % align trigger to downsampled signal
    if isfield(cfg,'resamplefs')
        if isstruct(trigger_nods)
            ts_trigger = round(trigger_nods(ntr).sample / (cfg.hdr.Fs/cfg.fsample)); %[~, ts_trigger] = min(abs(data.trial{1}(2,:) -
        else
            ts_trigger = round(trigger_nods(ntr) / (cfg.hdr.Fs/cfg.fsample));
        end
    else
        if isstruct(trigger_nods) % for intracranial
            ts_trigger = trigger_nods(ntr).sample;
        else % for EEG
            ts_trigger = trigger_nods(ntr);
        end
    end
    
    % start of the trial as defined in the log:
    trl_sid = ts_trigger - (log.values{tr, fieldID_trigonset} - log.values{tr, fieldID_trialonset})*cfg.fsample;
    
    % end of trial = beginning of new trial (unless it's the end of the
    % session)
    if tr == size(log.values,1) || log.values{tr+1,fieldID_session} ~= session
        % last trial of the session/phase
        if sum(strcmp(log.values{tr, fieldID_phase}, {'encoding'}))
            if trialinfo{3,ntr} == 1
                % DaD
                trl_eid = ts_trigger + (log.values{tr, fieldID_dadloconset} + log.values{tr,fieldID_dadRTloc} - log.values{tr, fieldID_triggeronset} + cfg.trialdef.postresp)*cfg.fsample;
            else
                % normal
                trl_eid = ts_trigger + (log.values{tr, fieldID_stimonset} + log.values{tr,fieldID_RTenc} - log.values{tr, fieldID_triggeronset} + cfg.trialdef.postresp)*cfg.fsample;
            end
        else
            if trialinfo{3,ntr} == 0
                % no catch question
                trl_eid = ts_trigger + (log.values{tr, fieldID_retreinst} +3 - log.values{tr, fieldID_trigonset} + cfg.trialdef.postresp)*cfg.fsample;
            else
                trl_eid = ts_trigger + (log.values{tr, fieldID_catchresponset} - log.values{tr, fieldID_trigonset} + cfg.trialdef.postresp)*cfg.fsample;
            end
        end
    else
        trl_eid = ts_trigger + (log.values{tr+1, fieldID_trialonset} - log.values{tr, fieldID_trialonset})*cfg.fsample;
    end
    
    % 0 time point in the trial (split per phase)
    % familiarization: object onset
    % encoding: object onset
    % retrieval: cue onset
    if sum(strcmp(log.values{tr, fieldID_phase}, {'familiarization'}))
        trl_t1 = log.values{tr, fieldID_trialonset} - log.values{tr, fieldID_stimonset};
        trl_cue = NaN;
        trl_object = 0;
        trl_response = NaN;
        trl_catchon = log.values{tr, fieldID_catchonset} - log.values{tr, fieldID_stimonset};
        trl_catchresp = log.values{tr, fieldID_catchresponset} - log.values{tr, fieldID_stimonset};
    elseif sum(strcmp(log.values{tr, fieldID_phase}, {'encoding'}))
        if trialinfo{3,ntr} == 1
            % dad trials
            trl_t1 = log.values{tr, fieldID_trialonset} - log.values{tr, fieldID_dadstimonset};
            trl_cue = NaN;
            trl_object = 0;
            trl_response = log.values{tr, fieldID_dadRTreinst};
            trl_catchon = NaN;
            trl_catchresp = NaN;
        else
            % normal encoding trials
            trl_t1 = log.values{tr, fieldID_trialonset} - log.values{tr, fieldID_stimonset};
            trl_cue = log.values{tr, fieldID_cueonset} - log.values{tr, fieldID_stimonset};
            trl_object = 0;
            trl_response = log.values{tr,fieldID_RTenc};
            trl_catchon = NaN;
            trl_catchresp = NaN;
        end
    elseif sum(strcmp(log.values{tr, fieldID_phase}, {'retrieval'}))
        trl_t1 = log.values{tr, fieldID_trialonset} - log.values{tr, fieldID_cueonset};
        trl_cue = 0;
        trl_object = NaN;
        trl_response = log.values{tr,fieldID_RTret};
        trl_catchon = log.values{tr, fieldID_catchonset} - log.values{tr, fieldID_cueonset};
        trl_catchresp = log.values{tr, fieldID_catchresponset} - log.values{tr, fieldID_cueonset};
    end

    trl_sid = round(trl_sid);
    trl_eid = round(trl_eid);
    
    % store trial details
    event{1,ntr} = trl_cue; % verb onset
    event{2,ntr} = trl_object;% object onset
    event{3,ntr} = trl_response; % time of 'clear' image
    event{4,ntr} = trl_catchon;
    event{5,ntr} = trl_catchresp;
    
    trl(ntr,1) = trl_sid;
    trl(ntr,2) = trl_eid;
    trl(ntr,3) = round(trl_t1*cfg.fsample);%trialtime_full(1);

end
    
%% store trial info in FieldTrip struct

% update FieldTrip struct
cfg.trialdef.event = event;
cfg.trialdef.trialinfo = trialinfo;
cfg.trialdef.cue = cue;
cfg.trialdef.stimulus = stimulus;
cfg.trialdef.response = response;
cfg.trialdef.trl = trl;

end