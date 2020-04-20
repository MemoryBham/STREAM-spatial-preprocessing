
clear all

%% settings

datapath = {'W:\Marije\Data_STREAM-spatial_EEG\Logs\'};

savepath = 'W:\Marije\Data_STREAM-spatial_EEG\Preproc\';
savename = 'results_behavior';

catch_types = {'exemplar', 'perc_1', 'perc_2', 'sem_1', 'sem_2'};
catch_names = {'label','perc1','perc2','sem1','sem2'}; % label, photo/drawing, color/black-white, living, flying
ncatch = length(catch_names);

minRT = 0.2; %s
maxRTsd = 3; %SDs

%% identify log files

logfiles = [];
paths = {};
for p = 1:length(datapath)
    % find log files on path
    id1 = length(logfiles)+1;
    logfiles = cat(1,logfiles,dir([datapath{p},'resultfile*.txt']));
    % store log file names
    paths(id1:length(logfiles)) = cellstr(datapath{p});
end

nfiles = length(logfiles);

if isempty(logfiles)
    error('No logfiles found')
end

%% % preallocate variables

subjID = cell(nfiles,1);

avgRT_fam = nan(nfiles,1);
Perf_fam = nan(nfiles,1);
avgRT_fam_catch = zeros(nfiles,ncatch);
RT_fam_catch = cell(ncatch,1);
Perf_fam_catch = nan(nfiles,ncatch);

avgRT_enc = nan(nfiles,1);
Perf_enc = nan(nfiles,1);

avgRT_ret = nan(nfiles,1);
Perf_ret = nan(nfiles,1);
avgRT_ret_catch = nan(nfiles,ncatch);
RT_ret_catch = cell(ncatch,1);
Perf_ret_catch = nan(nfiles,ncatch);
avgRT_ret_catch_Rep1 = nan(nfiles,ncatch);
avgRT_ret_catch_Rep2 = nan(nfiles,ncatch);

%% loop over subjects and analyze behavior

for lf = 1:nfiles
    
    fprintf('Loading subject %i of %i\n', lf, nfiles)
    
    % load log
    cfg = [];
    cfg.logpath = [paths{lf}, logfiles(lf).name];
    log = getLog(cfg);
    
    % store subject ID
    subjID{lf} = logfiles(lf).name(12:end-4);
    
    % find relevant fields of the log file
    fieldID_phase = find(strcmp(log.hdr, 'block_state'));
    fieldID_RT = find(strcmp(log.hdr, 'RT_catch'));
    fieldID_acc = find(strcmp(log.hdr, 'catch_resp'));
    fieldID_catch = find(strcmp(log.hdr, 'catch_type'));
    fieldID_trial = find(strcmp(log.hdr, 'trial_type'));
    fieldID_enc = find(strcmp(log.hdr, 'RT_encoding'));
    fieldID_DaD = find(strcmp(log.hdr, 'DaD_numb_attempts'));
    fieldID_ses = find(strcmp(log.hdr, 'sessionID'));
    fieldID_count = find(strcmp(log.hdr, 'stim_counter'));
    fieldID_stim = find(strcmp(log.hdr, 'stim_id'));
    
    % FAMILIARIZATION
    
    % find rows for this task phase
    phaseID = find(strcmp(log.values(1:end,fieldID_phase), 'familiarization'));
    % find rows with correct answers
    correctID = find([log.values{phaseID,fieldID_acc}]==1);
    
    % find reaction times for correct trials
    RTs = [log.values{phaseID(correctID),fieldID_RT}];
    
    % remove outliers
    RTs(RTs>nanmean(RTs)+maxRTsd*nanstd(RTs)) = NaN;
    RTs(RTs<minRT) = NaN;
    
    % compute mean RT and performance
    avgRT_fam(lf) = nanmean(RTs);
    Perf_fam(lf) = length(correctID)/length(phaseID);
    
    for c = 1:length(catch_types)
        % find catch_types
        catchID = find(strcmp(log.values(phaseID,fieldID_catch),catch_types{c}));
        
        % RTs per question type
        avgRT_fam_catch(lf,c) = nanmean([log.values{phaseID(intersect(correctID,catchID)),fieldID_RT}]);
        RT_fam_catch{c,1} = cat(1,RT_fam_catch{c,1},[log.values{phaseID(intersect(correctID,catchID)),fieldID_RT}]');
        
        % performance per category
        Perf_fam_catch(lf,c) = length(intersect(correctID,catchID))/length(catchID);

    end
    
    % ENCODING
    phaseID = find(strcmp(log.values(1:end,fieldID_phase), 'encoding'));
    if ~isempty(phaseID)
        % find normal learning trials
        trialID = find([log.values{phaseID,fieldID_trial}]==0);
        RTs = [log.values{phaseID(trialID),fieldID_enc}];
        
        % remove outliers
        RTs(RTs> nanmean(RTs)+maxRTsd*nanstd(RTs)) = NaN;
        RTs(RTs<minRT) = NaN;
        
        avgRT_enc(lf) = nanmean(RTs);
        
        % find minitest trials
        trialID = find([log.values{phaseID,fieldID_trial}]==1);
        % correct on the first attempt:
        correctID = find([log.values{phaseID(trialID),fieldID_DaD}]==1);
        
        Perf_enc(lf) = length(correctID)/length(trialID);
    end
    
    % RETRIEVAL
    phaseID = find(strcmp(log.values(1:end,fieldID_phase), 'retrieval'));
    if ~isempty(phaseID)
        correctID = find([log.values{phaseID,fieldID_acc}]==1);
        
        RTs = [log.values{phaseID(correctID),fieldID_RT}];
        
        % remove outliers
        RTs(RTs> nanmean(RTs)+maxRTsd*nanstd(RTs)) = NaN;
        RTs(RTs<minRT) = NaN;
        
        avgRT_ret(lf) = nanmean(RTs);
        Perf_ret(lf) = length(correctID)/length(find([log.values{phaseID,fieldID_acc}]>=0));
        
        for c = 1:length(catch_types)
            % find catch_types
            catchID = find(strcmp(log.values(phaseID,fieldID_catch),catch_types{c}));
            
            % RTs per question type
            avgRT_ret_catch(lf,c) = nanmean([log.values{phaseID(intersect(correctID,catchID)),fieldID_RT}]);
            RT_ret_catch{c,1} = cat(1,RT_ret_catch{c,1},[log.values{phaseID(intersect(correctID,catchID)),fieldID_RT}]');
            
            % performance per category - this should be highly unintersting
            Perf_ret_catch(lf,c) = length(intersect(correctID,catchID))/length(catchID);
            
            % find the order in which the stimuli were shown for this
            % question
            stimID = [log.values{phaseID(catchID),fieldID_stim}];
            
            % find first and second questions
            [~,firsts] = unique(stimID);
            seconds = setdiff(1:length(stimID),firsts);
            
            avgRT_ret_catch_Rep1(lf,c) = nanmean([log.values{phaseID(intersect(correctID,catchID(firsts))),fieldID_RT}]);  
            avgRT_ret_catch_Rep2(lf,c) = nanmean([log.values{phaseID(intersect(correctID,catchID(seconds))),fieldID_RT}]);
        end
    end
    
end

%% store data

fprintf('Saving data...')
xslfields = {'A','B','C', 'D','E'};
for a = 1:5
    xlswrite([savepath,savename, '_fam.xls'],RT_fam_catch{a,1},1,xslfields{a});
    
    if ~isempty(RT_ret_catch{a,1})
        xlswrite([savepath,savename, '_ret.xls'],RT_ret_catch{a,1},1,xslfields{a});
    end
end
fprintf('Done!\n')

%% averages across catch types

% performance
figure
plot([Perf_fam,Perf_enc,Perf_ret]', '-o')
xlim([0,4])
ylim([0.5,1.05])
set(gca,'xtick',1:3,'xticklabel',{'Fam','Enc','Ret'})
ylabel('Performance')
title('Performance per subj')

% RTs
figure
plot([avgRT_fam,avgRT_enc,avgRT_ret]', '-o')
xlim([0,4])
% ylim([0.5,1.05])
set(gca,'xtick',1:3,'xticklabel',{'Fam','Enc','Ret'})
ylabel('Avg reaction time (s)')
title('Avg RT per subj')

%% plot per catch type

% performance per subject
figure;
subplot(121)
plot(Perf_fam_catch', 'o-');
ylim([0.5,1.05]); xlim([0,6])
set(gca,'xtick',1:length(catch_types),'xticklabel',catch_types);
ylabel('Performance')
title('Familiarization')
subplot(122)
plot(Perf_ret_catch', 'o-')
ylim([0.5,1.05]); xlim([0,6])
set(gca,'xtick',1:length(catch_types),'xticklabel',catch_types);
title('Retrieval')

% RT - avg per subject
figure;
subplot(121)
plot(avgRT_fam_catch', 'o-');
xlim([0,6])
% ylim([0.5,1.05]);
set(gca,'xtick',1:length(catch_types),'xticklabel',catch_types);
ylabel('Avg reaction time (s)')
title('Familiarization')
subplot(122)
plot(avgRT_ret_catch', 'o-')
xlim([0,6])
% ylim([0.5,1.05]);
set(gca,'xtick',1:length(catch_types),'xticklabel',catch_types);
title('Retrieval')


% RT scatter
figure;
subplot(121); hold on
for c = 1:5
    scatter(c+0.3*rand(1,length(RT_fam_catch{c})),[RT_fam_catch{c}], 20,[0,0,0]);
    hold on
    plot(c,nanmedian([RT_fam_catch{c}]), 'ro');
end
xlim([0,6])
set(gca,'xtick',1:length(catch_types),'xticklabel',catch_types);
ylabel('Reaction time (s)')
title(['Familiarization', ' (N<=',num2str(length(logfiles)*16),')'])
subplot(122); hold on
for c = 1:5
    scatter(c+0.3*rand(1,length(RT_ret_catch{c})),[RT_ret_catch{c}], 20,[0,0,0]);
    hold on
    plot(c,nanmean([RT_ret_catch{c}]), 'ro');
end
xlim([0,6])
set(gca,'xtick',1:length(catch_types),'xticklabel',catch_types);
title(['Retrieval', ' (N<=',num2str(length(logfiles)*2*16),')'])
