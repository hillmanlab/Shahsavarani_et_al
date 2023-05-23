% This script generates correlation maps across every recording session using 
% 10-second moving windows. These correlation maps will be utilized in 
% subsequent steps, specifically step4 and step5, for non-negative least 
% squares (NNLS) fits.
%
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu

%% initialize the directories
clear;clc;

%{
rootDIR = '/home/ss6238/hillman_servers/';
dataDIR = strcat(rootDIR,'enterprise/3/Bahar/organizing_data/data/'); % ROI data
corrDIR = strcat(rootDIR,'enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/step3_correlationMaps/'); % to save correlation maps
%}

rootDIR = ' ';
dataDIR = ' '; % ROI data
corrDIR = ' '; % to save correlation maps


mice = ["cm124","cm125","cm126","cm127","cm128"];


% 10-s window (200 frames)
windowlength = 200; % in frames
windowstep = 1;
%%
tic
for num = 5%:length(mice)
    
    % read all recording sessions for one mouse
    runs = dir(strcat(dataDIR,mice(num),'*.mat'));
    f = waitbar(0,'Processing your data...');
    for r = 12:length(runs)
        waitbar(r/length(runs),f, 'Processing your data...')
        sessionName = runs(r).name;
        fprintf('Processing %s ... \n', sessionName(1:12))
        
        % load the data
        load(strcat(dataDIR,sessionName))
        
        jRGECO = data.jrgeco;
        HbT = data.chbo + data.chbr;
        
        ss = size(jRGECO);
        
        m = 0;
        for i = 1:windowstep:ss(2)-windowlength
            m = m + 1;
            [rho_jRGECO(:,:,m),~] = corr(jRGECO(:,i+[0:windowlength-1])');
            [rho_HbT(:,:,m),~] = corr(detrend(HbT(:,i+[0:windowlength-1])'));
        end
        
        correlations.jRGECO = rho_jRGECO;
        correlations.HbT = rho_HbT;
        
        save(strcat(corrDIR,sessionName(1:12)),'-struct','correlations','-v7.3')
        clear sessionName jRGECO HbT ss 
        clear correlations 
    end
    close(f)
    
end
waitbar(1, 'Finished!')
toc