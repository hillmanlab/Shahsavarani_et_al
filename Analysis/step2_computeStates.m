% This script generates behavioral states for each mouse, similar to the 
% ones depicted in Figures 3B, 7A, S4, S6, and S7. The correlation maps are 
% calculated over a 10-second window, corresponding to 200 frames.
% 
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

%% initialize the directories
clear;clc;

%{
rootDIR = '/home/ss6238/hillman_servers/';
dataDIR = strcat(rootDIR,'enterprise/3/Bahar/organizing_data/data/');  % ROI data
boutsDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/runningBlocks/'); % running blocks created in step1
corrDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/states/runningBlocksWCorrMaps/');  % to save running blocks with correlation maps
stateDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/states/'); % to save behavioral states
%}

rootDIR = ' ';
dataDIR = ' ';  % ROI data
boutsDIR = ' '; % running blocks created in step1
corrDIR = ' ';  % to save running blocks with correlation maps
stateDIR = ' '; % to save behavioral states


mice = ["cm124","cm125","cm126","cm127","cm128"];

%% Compute correlation maps related to each state
%
% Correlation maps over each behavioral state are calculated for both neural
% and hemodynamic acitivty.
% The average of these maps is considered as the states.
% Five states are identified based on the mouse behavior: locomotion onset,
% locomotion (runnig), locomotion offset, resting states (post-offset) including
% initial rest (right after mouse stops running) and sustained state (when
% mouse has not been running for a while, i.e., 40 s), locomotion onset
%
% Within corrDIR, you need to have 4 folders called: onset, running,
% offset, and post-offset. The struct running blocks with new field related
% to the correlation maps for each state will be saved here.
%
% To be able to run this section, you need following auxiliary codes:
%       - computeRunningCorrelations()
%       - computeOffsetCorrelations()
%       - computePostOffsetCorrelations()
%       - computeOnsetCorrelations()


delay = 30; % 1.5 s delay for hemodynamic data
for i = 1:length(mice)
    mousename = mice(i);
    
    %
    % load running bouts
    load(strcat(boutsDIR,mousename,'.mat'))
    
    %
    % compute running correlation maps
    computeRunningCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
    
    % compute offset correlation maps
    computeOffsetCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
    %
    % comput post-running correlation maps
    % 10-s window
    computePostOffsetCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
    %
    % compute onset correlation maps
    computeOnsetCorrelations(mousename,runningBlocks,corrDIR,dataDIR,delay)
    %
    
    clear runningBlocks
end
%}


%% Compute each state
% This section computes each behavioral state, averaging the correlation
% maps computed in the previous section. It computes the behavioral states
% for both neural and hemodynamic activity. You need to set the modality by
% variable mod where the value of jRGECO computes the neuronal states and the
% value of Hbt computes the hemodynamic states.
%
% To be able to run this section, you need following auxiliary codes:
%       - computeRunningState()
%       - computeOffsetState()
%       - computePostOffsetStates()
%       - computeOnsetState()

% choose the modality
mod = 'HbT'; %{'jRGECO','HbT'}

f = waitbar(0,'Processing your data');
for i = 1:length(mice)
    
    mousename = mice(i);
    
    % locomotion state
    clear runningBlocks
    load(strcat(corrDIR,'running/',mousename,'.mat'))
    runningState = computeRunningState(runningBlocks,mod);
    
    % offset state
    clear runningBlocks
    load(strcat(corrDIR,'offset/',mousename,'.mat'))
    offsetState = computeOffsetState(runningBlocks,mod);
    
    
    % initial-rest and sustained-rest states
    clear runningBlocks
    windows = [1,2];
    load(strcat(corrDIR,'post-offset/',mousename,'.mat'))
    postRunningStates = computePostOffsetStates(runningBlocks,windows,mod);
    
    % onset state
    clear runningBlocks
    load(strcat(corrDIR,'onset/',mousename,'.mat'))
    onsetState = computeOnsetState(runningBlocks,mod);
   
    states(:,:,1) = runningState;
    states(:,:,2) = offsetState;
    states(:,:,3) = postRunningStates(:,:,1);
    states(:,:,4) = postRunningStates(:,:,2);
    states(:,:,5) = onsetState;
    
    finalStates(i).mouseName = mousename;
    finalStates(i).states = states;
    
    clear mousename onsetState runningState postRunningStates states
    clear offsetstate
    
    waitbar(i/length(mice),f,'Processing your data');
    
end
save(strcat(stateDIR,'states_',mod),'finalStates')
delete(f)

