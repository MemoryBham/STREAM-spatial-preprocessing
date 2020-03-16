function log = getLog(cfg)

% LOG contains: 
% subID
% sessionID
% block_state: encoding or retrieval
% trial_id: number of trial in the presentation sequence
% trial_type: 0 = normal; 1 = drag-and-drop
% cue_id: location number
% cue_xcoord: (pixels)
% cue_ycoord: (pixels)
% stim_label: e.g. german shepherd
% stim_id: stimulus number
% stim_filename: xxx.jpg
% stim_counter: counter of repetitions for each stimulus
% stim_cat: e.g. dog
% stim_perc1: 1 = drawing; 2 = photograph;
% stim_perc2: 1 = left-facing; 2 = right-facing;
% stim_sem1: 1 = animate; 2 = inanimate
% stim_sem2: 1 = flying; 2 = non-flying;
% RT_encoding (s)
% DaD_resp: 0 = forgotten; 1 = remembered
% RT_DaD_reinst (s)
% DaD_numb_attempts: number of attempt to reach correct location
% RT_DaD_loc (s)
% ret_reinst_resp: 0 = forgotten; 1 = remembered
% RT_ret_reinst
% catch_type: perc1, perc2, sem1, sem2, exemplar
% RT_catch (s)
% catch_resp: 0 = incorrect; 1 = correct; 2 = forgotten; 3 no answer
% onset_session (s)
% onset_familiarization (s)
% onset_encoding (s)
% onset_retrieval (s)
% onset_trial (s)
% onset_cue (s)
% onset_trigger (s)
% onset_stim (s)
% onset_DaD_stim (s)
% onset_DaD_reinst (s)
% onset_DaD_loc (s)
% onset_ret_reinst (s)
% onset_catch (s)
% onset_catch_resp (s)

%% open log file

path = cfg.logpath;

% load log header and data
log = [];
fext = path(end-2:end);
if ~strcmp(fext,'txt')
    error('Unsupported file type')
elseif strcmp(fext,'txt')
    data = tdfread(path);
    log.hdr = fieldnames(data);
    log.values = cell(size(data.subID,1),size(log.hdr,1));
    for i = 1:length(log.hdr)
        if isnumeric(eval(['data.',log.hdr{i}, '(1,:)']))
             log.values(:,i) = num2cell(eval(['data.',log.hdr{i}, '(1:end,:)']));
        else
            if ismember(log.hdr{i},{'subID','block_state','stim_label',...
                    'stim_filename','stim_cat','catch_type'})
                log.values(:,i) = cellstr(eval(['data.',log.hdr{i}, '(1:end,:)']));
            else
                for j = 1:size(data.subID,1)
                    log.values{j,i} = str2num(eval(['data.',log.hdr{i}, '(j,:)']));
                end
            end
        end
    end
end

