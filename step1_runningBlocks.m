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

rootDIR = ' ';
dataDIR = ' '; % ROI data
saveDIR = ' ';

mousename = 'cm128'; %{'cm124','cm125','cm126','cm127','cm128'}
runs = dir(strcat(dataDIR,mousename,'*.mat'));


% we use the variable idx to count the number of running blocks
idx = 0;
runningBlocks = [];

%% load one experimental session, one recording session at a time
for n = 1:length(runs)
    
    mouse = runs(n).name;
    
    %
    % remove these runs due to excessive grooming
    if (strcmp(mouse(1:12),'cm126_6_runC') || strcmp(mouse(1:12),'cm126_6_runD') ...
            || strcmp(mouse(1:12),'cm126_6_runI') || strcmp(mouse(1:12),'cm126_6_runJ'))
        continue
        
    elseif strcmp(mouse(1:12),'cm127_1_runH')
        continue
        
    elseif (strcmp(mouse(1:12),'cm128_1_runH') || strcmp(mouse(1:12),'cm128_1_runJ') ...
            || strcmp(mouse(1:12),'cm128_7_runI') || strcmp(mouse(1:12),'cm128_5_runI'))
        continue
        
    end
    %}
    
    mousename = mouse(1:5);
    day = mouse(7);
    run = mouse(9:12);
    
    
    % load data
    load(strcat(dataDIR,mouse))
    
    %rotf = lowpass(abs(m.rotf),0.5,20);
    rotf = info.behavior.wheelVelocity;
    
    rz = getrunningpulses2(rotf,0.1,20,100,"off");
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
    
    prev_length = length(runningBlocks);
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
    for i = prev_length+1:length(runningBlocks)
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
    
    % clean up and remove grooming
    switch mousename
        case 'cm126'
            if strcmp(mouse(1:12),'cm126_2_runJ')
                runningBlocks = runningBlocks(1:prev_length+1); % the rest is grooming
            elseif strcmp(mouse(1:12),'cm126_3_runC')
                runningBlocks = runningBlocks([1:prev_length,prev_length+4:length(runningBlocks)]); % the first three blocks are grooming
            elseif strcmp(mouse(1:12),'cm126_3_runH')
                runningBlocks = runningBlocks([1:prev_length,prev_length+5:length(runningBlocks)-1]); % excluding grooming
            elseif strcmp(mouse(1:12),'cm126_3_runI')
                runningBlocks = runningBlocks([1:prev_length,prev_length+4:length(runningBlocks)]); % excluding grooming
            elseif strcmp(mouse(1:12),'cm126_5_runD')
                runningBlocks = runningBlocks([1:prev_length,prev_length+1:length(runningBlocks)-3]); % excluding grooming at the end
            end
          
            
        case 'cm128'
            if strcmp(mouse(1:12),'cm128_4_runC')
                runningBlocks = runningBlocks(1:prev_length+3); % the rest is grooming
            elseif strcmp(mouse(1:12),'cm128_5_runD')
                runningBlocks = runningBlocks(1:prev_length+4); % the rest is grooming
            elseif strcmp(mouse(1:12),'cm128_5_runH')
                runningBlocks = runningBlocks(1:prev_length+8); % the rest is grooming
            elseif strcmp(mouse(1:12),'cm128_7_runH')
                runningBlocks = runningBlocks(1:prev_length+6); % the rest is grooming
            elseif strcmp(mouse(1:12),'cm128_8_runH')
                runningBlocks = runningBlocks(1:prev_length+4); % the rest is grooming
            end
    end
    idx = length(runningBlocks);
    
    
    clearvars -except Hd mousename runs idx n runningBlocks T saveDIR dataDIR
end

% save runningBlocks
save(strcat(saveDIR,mousename), 'runningBlocks', '-v7.3')

