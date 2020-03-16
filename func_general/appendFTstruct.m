
%%%%%%%%%%%% appendFTstruct
% DESCRIPTION:
% This functions identifies appends trial cut FieldTrip structs that include the
% additional fields:
% - cfg.hdr
% - cue
% - stimulus
% - response
% It works on spike, LFP and spectrum data struct (determined automatically) 
% and can extent along one of two dimensions (specified by appendType): 
% number of channels, or number of trials

% INPUT:
% - in_base: input struct to be extended
% - in_new: input struct to be added to in_base
% - appendType: dimension to append along: 'channel' or 'trial'

% OUTPUT:
% appended output struct of the same format as in_base

function out = appendFTstruct(in_base, in_new, appendType)

out = in_base; 
if isfield(in_base, 'powspctrm')
    datatype = 'spectrum';
    out.label = in_base.label(:,1);
elseif isfield(in_base, 'timestamp')
    datatype = 'spikes';
%     out.label = in_base.label(:,1);
else 
    datatype = 'LFPs';
end

switch appendType
    case 'channel' 
        if strcmp(datatype, 'LFPs') || strcmp(datatype,'spectrum')
            out.hdr.label = vertcat(out.hdr.label, in_new.hdr.label);
            out.hdr.nChans = out.hdr.nChans +1;
            out.hdr.chantype = vertcat(out.hdr.chantype, in_new.hdr.chantype);
            out.hdr.chanunit = vertcat(out.hdr.chanunit, in_new.hdr.chanunit);
        end
        if strcmp(datatype, 'spikes')
            out.trial = [out.trial,in_new.trial];
            out.time = [out.time,in_new.time];
            out.timestamp = [out.timestamp,in_new.timestamp];
            out.waveform = [out.waveform,in_new.waveform];
            
            for i = 1:length(in_new.label)
                try
                    dum = sum(ismember(out.label, in_new.label{i}));
                catch
                    dum = sum(ismember([out.label{:}], in_new.label{i}));
                end
                if dum
                    warning('Label duplicate, creating new label')
                    dumv = [out.label,in_new.label];
                    dum1 = regexp(dumv, ...
                        regexptranslate('wildcard',[in_new.label{i}(1:end-1),'*']));
                    dum2 = find(arrayfun(@(x) isempty(dum1{x}), 1:length(dum1))==0);
                    lastname = max(arrayfun(@(x) str2num(dumv{dum2(x)}(end)), 1:length(dum2)));
                    in_new.label{i} = [in_new.label{i}(1:end-1), num2str(lastname+1)];
                end
            end
            out.label = [out.label, in_new.label];
        elseif strcmp(datatype, 'LFPs')
            out.label{end+1,1} = in_new.label{:,1};
            for tr = 1:length(out.trial)
                out.trial{1,tr} = [out.trial{1,tr};in_new.trial{1,tr}];
            end
        elseif strcmp(datatype, 'spectrum')
            out.label{end+1,1} = in_new.label{:,1};
            for tr = 1:length(out.powspctrm)
                out.powspctrm{1,tr} = [out.powspctrm{1,tr};in_new.powspctrm{1,tr}];
            end
        end
        
    case 'trial'
        
        out.trialinfo = horzcat(out.trialinfo,in_new.trialinfo);
        out.sampleinfo = vertcat(out.sampleinfo,in_new.sampleinfo);
        if isfield(out.cfg,'trl')
            trl_new = in_new.cfg.trl;
            trl_new(:,1) = trl_new(:,1) + (out.cfg.trl(end,1) + 1 - trl_new(1,1));
            trl_new(:,2) = trl_new(:,2) + (out.cfg.trl(end,2) + 1 - trl_new(1,2));
            out.cfg.trl = vertcat(out.cfg.trl, trl_new);
        end
        if strcmp(datatype, 'spikes')
             for u = 1:length(in_new.label)
                out.trial{1,u} = horzcat(out.trial{1,u},in_new.trial{1,u});
                out.time{1,u} = horzcat(out.time{1,u},in_new.time{1,u});
                out.timestamp{1,u} = horzcat(out.timestamp{1,u},in_new.timestamp{1,u});
                out.waveform{1,u} = cat(3,out.waveform{1,u},in_new.waveform{1,u});
            end
            out.trialtime = vertcat(out.trialtime,in_new.trialtime);
        elseif strcmp(datatype, 'LFPs')
            out.hdr.nTrials = out.hdr.nTrials +1;
            out.trial = horzcat(out.trial,in_new.trial);
            out.time = horzcat(out.time,in_new.time);
        elseif strcmp(datatype, 'spectrum')
            out.hdr.nTrials = out.hdr.nTrials +1;
            out.powspctrm = horzcat(out.powspctrm,in_new.powspctrm);
            out.time = horzcat(out.time,in_new.time);
        end
        
        if isfield(out, 'event') %those field do not exist for baseline structs
            out.event = horzcat(out.event,in_new.event);
            out.cue = horzcat(out.cue,in_new.cue);
            out.stimulus = horzcat(out.stimulus,in_new.stimulus);
            if size(out.response,1) ~= size(in_new.response,1)
                out.response = horzcat(out.response(1:size(in_new.response,1),:),in_new.response);
            else
                out.response = horzcat(out.response,in_new.response);
            end
        end
        
end
end