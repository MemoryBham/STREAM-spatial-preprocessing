
% add fieldtrip to path
addpath C:\Users\marij\Documents\MATLAB\Toolboxes\fieldtrip-20180918
% addpath C:\toolbox\fieldtrip-20160601
ft_defaults

addpath(genpath('STREAM-spatial-preproc\func_EEG\'))
addpath(genpath('STREAM-spatial-preproc\func_general\'))

% comment this out when running multiple participants!
path_basic = 'W:\Marije\Data_STREAM-spatial_EEG\';
path_data = [path_basic,'EEG\'];
path_preproc = [path_basic,'Data_preproc\'];
path_log = [path_basic,'\Logs\'];
subjID = 'subYW01';%

%% load data

if exist([path_preproc, 'dataEEG_',subjID,'.mat'], 'file')
    fprintf('Loading data for artifact rejection\n')
    load([path_preproc, 'dataEEG_',subjID,'.mat'])
else
    cfg                         = [];
    cfg.path                    = path_data;
    cfg.logpath                 = path_log;
    cfg.subject                 = subjID;
    cfg.numbSessions            = 2;
    cfg.channel                 = 1:128;
    cfg.lowmemory               = true; % use for systems with low memory
    
    cfg.trialdef.eventtype      = 'STATUS';
    cfg.trialdef.eventvalue     = [65281];
    cfg.trialdef.postresp       = 2; % in seconds
    
    cfg.hpfilter                = 'yes';
    cfg.hpfreq                  = 0.5;
    % cfg.hpfiltord             = 3;
    cfg.hpfilttype              = 'fir';%
    cfg.hpfiltdir               = 'twopass';
    
    cfg.lpfilter                = 'yes';
    cfg.lpfreq                  = 100;
    % cfg.lpfiltord             = 6;
    cfg.lpfilttype              = 'fir';%
    cfg.lpfiltdir               = 'twopass'; 
    
    cfg.detrend                 = 'no';
    cfg.demean                  = 'yes';
    cfg.baselinewindow          = 'all';%
    cfg.feedback                = 'no';
    cfg.resamplefs              = 512;% Hz
    
    cfg.notchfilter             = 'yes';
    cfg.notchfreq               = 50;
    cfg.notchend                = 200;
    cfg.notchbw                 = 0.5;
    
    dataEEG                     = getEEG(cfg);
    
    save([path_preproc, 'dataEEG_', subjID, '.mat'],'dataEEG', '-v7.3')
end

%% split the data

% ---- select based on trial info

% choose a session
cfg = [];
cfg.session = 2; % only session 2
dataEEGs2 = selectData_STREAMspatial(cfg, dataEEG);

% choose a task phase
cfg = [];
cfg.phase = 1; % 1 = familiarization; 2 = encoding; 3 = retrieval
% ___ OR ___
cfg.phase = 'fam'; % 'fam' / 'familiarization' / 'enc' / 'encoding' / 'ret' / 'retrieval'
dataFam = selectData_STREAMspatial(cfg, dataEEG);

% choose specific trial types
cfg = [];
cfg.phase = 2; % encoding trials only
cfg.trialtype = 0; % use only normal encoding trials, i.e. remove mini-test trials
% ___ OR ___
cfg.trialtype = 1; % keep only mini-test trials
dataEnc = selectData_STREAMspatial(cfg, dataEEG);

% choose stimulus characteristics
cfg = [];
cfg.perc1 = 2; % keep drawings only
dataDraw = selectData_STREAMspatial(cfg, dataEEG);

cfg = [];
cfg.category = 'dogs'; % keep dogs only - 'dogs' / 'birds' / 'cars' /  'planes'
dataDogs = selectData_STREAMspatial(cfg, dataEEG);

% select trials with specific catch questions
cfg = [];
cfg.phase = 3; % retrieval only
cfg.trialtype = [4,5]; % keep the semantic questions
dataQSem = selectData_STREAMspatial(cfg, dataEEG);

% keep remembered trials only
cfg = [];
cfg.phase = 'ret';
cfg.remembered = 1; % remembered trials - 0 for forgotten
dataRetRem = selectData_STREAMspatial(cfg, dataEEG);

% keep trials with correct catch questions only
cfg = [];
cfg.phase = 'fam';
cfg.correct = 1; % correct trials - 0 for incorrect
dataFamCorrect = selectData_STREAMspatial(cfg, dataEEG);

% remove button presses that are too fast or too slow
% the numbers below are arbitary!
cfg = [];
cfg.phase = 2; % encoding trials
cfg.minRT = 0.3; % in s; discard fast encoding button presses
cfg.maxRT = 10; % s; discard slow encoding button presses
dataEncSel = selectData_STREAMspatial(cfg, dataEEG);

cfg = [];
cfg.phase = 1; % encoding trials
cfg.minCatchRT = 0.2; % s; discard fast encoding button presses
cfg.maxCatchRT = 5; % s; discard slow encoding button presses
dataFamSel = selectData_STREAMspatial(cfg, dataEEG);

% ---- select time points
cfg = [];
cfg.phase = 3; % retrieval trials
cfg.event = 3; % reinstatement button press
cfg.toi = [0,3]; % 3 second window starting at reinstatement button press
dataRetImag = selectData_STREAMspatial(cfg, dataEEG);

cfg = [];
cfg.phase = 3; % retrieval trials
cfg.minRT = 0.3; % discard fast reinstatement button presses
cfg.event = 3; % reinstatement button press
cfg.toi = [-1,0]; % 1 second window leading up to reinstatement button press
dataRetPreResp = selectData_STREAMspatial(cfg, dataEEG);

% events:
% 1: cue onset
% 2: object onset
% 3: response onset (for encoding, this is the button press to go to the
% next trial, for retrieval to indicate reinstatement)
% 4: catch question onset (fam and ret only)
% 5: catch response (fam and ret only)

