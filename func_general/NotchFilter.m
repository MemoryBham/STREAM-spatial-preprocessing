function data = NotchFilter(cfg,data)

fprintf('Notch filtering\n')
if isfield(cfg, 'bpfreq')
    steps = cfg.notchfreq:cfg.notchfreq:cfg.bpfreq(2); 
elseif isfield(cfg, 'lpfreq')
    steps = cfg.notchfreq:cfg.notchfreq:cfg.lpfreq;
elseif isfield(cfg, 'notchend')
    steps = cfg.notchfreq:cfg.notchfreq:cfg.notchend;
else
    error('No end frequency for Notch filter')
end

a = 1; b = 1;
for i = 1:length(steps)
    [bt,at] = butter(1,2*[steps(i)-cfg.notchbw, steps(i)+cfg.notchbw]/data.fsample, 'stop');
    if i>1
        b = conv(b,bt);
        a = conv(a,at);
    else
        b = bt; 
        a = at;
    end
end
% fvtool(b,a, 'Fs', cfg.resamplefs)

for tr = 1:length(data.trial)
    for ch = 1:size(data.trial{1},1)
    data.trial{1,tr}(ch,:) = filtfilt(b,a,data.trial{1,tr}(ch,:)')';
    end
end

end

