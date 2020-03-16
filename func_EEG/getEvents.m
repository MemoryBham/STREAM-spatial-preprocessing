function [events_clean, events_deleted] = getEvents(cfg)

fprintf('\nFinding events\n');

cfg = ft_definetrial(cfg);

%% find events 

% time points with trigger value of interest
idx         = cfg.trl(:,1);%find([events.value] == cfg.trialdef.eventvalue);
% 
idx_cleanstep1 = setdiff(idx, idx([find(diff(idx) < 1.5*0.15*cfg.hdr.Fs)+1; ...
    find(diff(idx) < 1.5*0.15*cfg.hdr.Fs)]));

% compare events with log - delete extra's
if length(idx_cleanstep1) == cfg.ntrials
    fprintf('Number of events matches expectation - starting data epoching\n')
    idx_cleanstep2 = idx_cleanstep1;
elseif length(idx_cleanstep1) < cfg.ntrials
    error('Fewer events found than expected, check the log-file and trigger')
    idx_cleanstep2 = idx_cleanstep1;
else 
    ans1 = input('More events found than expected; Delete events? y/n \n', 's');
    if ismember(ans1,{'y', 'Y', 'yes'})
        prmt = ['Select the trials you want to delete, choose ',...
            num2str(abs(length(idx_cleanstep1)-cfg.ntrials)),...
            ' trial numbers between 1 and ', num2str(length(idx_cleanstep1)), '\n'];
        ans2 = input(prmt);
        % delete trials in ans2
        idx_cleanstep2 = idx_cleanstep1(setdiff(1:length(idx_cleanstep1),ans2));
    else
        error('Could not align log and trigger - preprocessing aborted\n')  
    end
end

events_clean = idx_cleanstep2;
events_deleted = setdiff(idx, idx_cleanstep2);

end