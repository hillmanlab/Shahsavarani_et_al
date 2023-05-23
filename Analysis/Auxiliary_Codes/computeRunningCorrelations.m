function runningBlocks = computeRunningCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
% This function computes the correlation map over the running window for each
% running block. The running window is defined as 10 s middle of locomotion.
% For each mouse, it returns a struct of running blocks with two new fields: 
% correlations_neural and correlations_HbT, which contain correlation maps
% over the running winow for neural and hemo data.
% 
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu
     
% remove running bouts with less than 20 s duration
duration = [runningBlocks.duration];
durationIDX = duration < 400;
runningBlocks(durationIDX) = [];

% for FC maps of 10 seconds, set it to 100 for half a window
nframes = 100;

     
for run = 1:length(runningBlocks)
    sessionName = runningBlocks(run).session;
    
    fprintf('run %d ... \n', run)
    % load the data
    load(strcat(dataDIR,sessionName))
    
    jrgeco = data.jrgeco;
    chbt = data.chbo + data.chbr;
        
    duration = runningBlocks(run).duration;
    onset = runningBlocks(run).onset;
    
    a = onset+floor(duration/2)-nframes;
    b = onset+floor(duration/2)+nframes;
    
        
    rho_jRGECO = corr(jrgeco(:,a:b)');
    rho_HbT = corr(detrend(chbt(:,a+delay:b+delay)'));
    
    runningBlocks(run).correlations_neural = rho_jRGECO;
    runningBlocks(run).correlations_HbT = rho_HbT;
    
    
    clear jrgeco chbt
    clear rho_jRGECO rho_HbT
    clear onset
    
    
end

%
save(strcat(corrDIR,'running/',mousename),'runningBlocks','-v7.3')

end