function runningBlocks = computeOnsetCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
% This function computes the correlation map over the onset window for each
% running block. The onset window overlaps locomotion and rest with 5 s
% before locomotion and 5 s beginning of locomotion.
% For each mouse, it returns a struct of running blocks with two new fields: 
% correlations_neural and correlations_HbT, which contain correlation maps
% over the onset winow for neural and hemo data.
% 
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

% remove running bouts with less than 10 s duration
duration = [runningBlocks.duration];
durationIDX = duration < 200;
runningBlocks(durationIDX) = [];

% remove running bouts with pre-running rest time less than 60 s
time2previousrun = [runningBlocks.time2previousrun];
time2previousrunIDX = time2previousrun < 1200;
runningBlocks(time2previousrunIDX) = [];

% for FC maps of 10 seconds, set it to 100 for half a window
nframes = 100;

for run = 1:length(runningBlocks)
    sessionName = runningBlocks(run).session;
    
    fprintf('run %d ... \n', run)
    % load the data
    load(strcat(dataDIR,sessionName))
    
    jrgeco = data.jrgeco;
    chbt = data.chbo + data.chbr;
        
    onset = runningBlocks(run).onset;
    
    a = onset-nframes;
    b = onset+nframes;
    
    rho_jRGECO = corr(jrgeco(:,a:b)');
    rho_HbT = corr(detrend(chbt(:,a+delay:b+delay)'));
    
    runningBlocks(run).correlations_neural = rho_jRGECO;
    runningBlocks(run).correlations_HbT = rho_HbT;
    
    
    clear jrgeco chbt
    clear rho_jRGECO rho_HbT
    clear onset
    
end

%
save(strcat(corrDIR,'onset/',mousename),'runningBlocks','-v7.3')

end