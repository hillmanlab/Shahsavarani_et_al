% This script generates running blocks for each mouse. To run it 
% successfully, you will need the auxiliary code called "getrunningpulses2()".
% You will generate a struc for each mouse, which consists of 10 fields:
% session, mouse_name, day, run, pulse signal corresponding to the
% locomotion signal, onset of running, offset of running, duration of
% running, time to the preveiouse run, and time to the nex run.
%
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

%% initialize the directories
clear;clc

dataDIR = ' '; % ROI data

mousename = 'cm124'; %{'cm124','cm125','cm126','cm127','cm128'}
runs = dir(strcat(dataDIR,mousename,'*.mat'));


% we use the variable idx to count the number of running blocks
idx = 0;

%% load one experimental session, one recording session at a time
for n = 1:length(runs)
    
    mouse = runs(n).name;
    
    mousename = mouse(1:5);
    day = mouse(7);
    run = mouse(9:12);
    
    
    % load data
    load(mouse)
    
    %rotf = lowpass(abs(m.rotf),0.5,20);
    rotf = m.rotf;
    
    rz = getrunningpulses2(rotf',0.1,20,100,"off");
    %figure;plot(m.rotf);hold on;plot(rz)
    
    
    % detect the running blocks (pulses)
    rz_diff = diff(rz);
    % detect rising edges
    rise_edge_idx = find(rz_diff == 1);
    % detect falling edges
    fall_edge_idx = find(rz_diff == -1);
    
    
    % There are four possible conditions for running pulses:
    % (1) all pulses have the rise edge and fall edge ==>   ____|----|______|----|______|----|_______
    % (2) the first pulse does not have the rise edge ==>   ----|______|----|______|----|_______
    % (3) the last pulse does not have the fall edge  ==> ____|----|______|----|______|---------
    % (4) the first pulse does not have the rise edge and the last pulse
    % does not have the fall edge ==> -------|______|----|______|--------
    
    % We treat the conditions in a way as they are condition (4).
    % That is we ignore the first and the last pulse of each run.
    
    % number of blocks to consider
    if length(rise_edge_idx) == length(fall_edge_idx)
        if rise_edge_idx(1) > fall_edge_idx(1) % condition (4)
            rise_edge_idx = [0;rise_edge_idx];
            fall_edge_idx = [fall_edge_idx;0];
        end
    else
        if rise_edge_idx(1) > fall_edge_idx(1)
            rise_edge_idx = [0;rise_edge_idx];
        else
            fall_edge_idx = [fall_edge_idx;0];
        end
    end
    numBlocks = length(rise_edge_idx);
    
    
    onset = rise_edge_idx + 1;
    offset = fall_edge_idx;
    
    for j = 2:numBlocks-1 % the first and the last running blocks are ignored
        
        startpoint = offset(1) + 1;
        endpoint = onset(end) - 1;
        
        
        idx = idx + 1;
        
        runningBlocks(idx).session = mouse;
        runningBlocks(idx).mouse_name = mousename;
        runningBlocks(idx).day = day;
        runningBlocks(idx).run = run;
        
        runningBlocks(idx).rz = rz;
        
        runningBlocks(idx).onset = onset(j);
        runningBlocks(idx).offset = offset(j);
        runningBlocks(idx).duration = offset(j) - onset(j);
        runningBlocks(idx).time2previousrun = onset(j) - offset(j-1);
        runningBlocks(idx).time2nextrun = onset(j+1) - offset(j);
        
    end
    
    % combine the blocks with less than 5-s rest (100 frame)
    % this part approximately corrects for the multiple pulses idetified
    % within one single running bout
    for i = 1:length(runningBlocks)
        if i > length(runningBlocks)
            break
        end
        time2nextrun = runningBlocks(i).time2nextrun;
        
        while time2nextrun < 100
            % update offset
            try
            runningBlocks(i).offset = runningBlocks(i+1).offset;
            catch
                runningBlocks(i) = [];
                break
            end
            
            %update time2nextrung
            runningBlocks(i).time2nextrun = runningBlocks(i+1).time2nextrun;
            
            runningBlocks(i+1) = [];
            time2nextrun = runningBlocks(i).time2nextrun;
        end
    end
    idx = length(runningBlocks);
            
    
    clearvars -except Hd mousename runs idx n runningBlocks T saveDIR
end

% save runningBlocks
save(strcat(saveDIR,mousename), 'runningBlocks', '-v7.3')

