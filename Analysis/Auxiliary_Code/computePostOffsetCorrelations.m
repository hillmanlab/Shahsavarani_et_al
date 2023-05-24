function runningBlocks = computePostOffsetCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
% This function computes the correlation maps over the resting-state windows 
% for each running block. The resting states are initial rest and sustained
% rest. The initial rest is defined as 10 s rest immediately after
% locomotion offset. The sustained rest is defined as 10 s rest starting
% from 40 s after locomotion offset.
% For each mouse, it returns a struct of running blocks with two new fields: 
% correlations_neural and correlations_HbT, which contain correlation maps
% over the resting-state winows for neural and hemo data.
% 
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

% remove running bouts with less than 5 s duration
duration = [runningBlocks.duration];
durationIDX = duration < 100;
runningBlocks(durationIDX) = [];

% remove running bouts with post-running rest time less than 60 s
time2nextrun = [runningBlocks.time2nextrun];
time2nextrunIDX = time2nextrun < 1200;
runningBlocks(time2nextrunIDX) = [];

startFrame1 = 1; % the start frame affter offset for the initial rest (immediately after run ends)
startFrame2 = 800; % the start frame affter offset for the sustained rest (40 s after run ends)

% for FC maps of 10 seconds, set it to 100 for half a window
nframes = 100;

for run = 1:length(runningBlocks)
    sessionName = runningBlocks(run).session;
    
    fprintf('run %d ... \n', run)
    % load the data
    load(strcat(dataDIR,sessionName))
    
    
    jrgeco = data.jrgeco;
    chbt = data.chbo + data.chbr;    

    offset = runningBlocks(run).offset;
    
    a = offset+startFrame1;
    b = offset+startFrame1+2*nframes;
    
    c = offset+startFrame2;
    d = offset+startFrame2+2*nframes;
    
    % correlation map over the initial rest
    rho_jRGECO(1,:,:) = corr(jrgeco(:,a:b)');
    rho_HbT(1,:,:) = corr(detrend(chbt(:,a+delay:b+delay)'));
    
    % correlation map over the sustained rest
    rho_jRGECO(2,:,:) = corr(jrgeco(:,c:d)');
    rho_HbT(2,:,:) = corr(detrend(chbt(:,c+delay:d+delay)'));
    
    runningBlocks(run).correlations_neural = rho_jRGECO;
    runningBlocks(run).correlations_HbT = rho_HbT;
    
    
    clear jrgeco chbt
    clear rho_jRGECO rho_HbT
    clear offset
    
end

%
save(strcat(corrDIR,'post-offset/',mousename),'runningBlocks','-v7.3')

end