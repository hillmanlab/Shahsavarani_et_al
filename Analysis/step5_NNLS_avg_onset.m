% This script calculates the NNLS coefficients for the running blocks with
% a 60-second rest period before running and a 20-second running duration. 
% Additionally, it plots the average of these coefficients, similar to 
% Figure 3E. To run this script, you need the running blocks created in 
% step1, the states computed in step2, and the correlation maps generated 
% in step3.
%
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu
%% initialize the directories
clearvars -except nnlsqDistance
clc

%{
behaviorDIR = '/home/ss6238/hillman_servers/chronos/2/cmdata_analysis/RS_analysis/';     % behavior data
boutsDIR = '/home/ss6238/hillman_servers/enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/step2_runningBlocks/'; % running blocks created in step1
stateDIR = '/home/ss6238/hillman_servers/enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/step1_states/';         % states computed in step2 
corrMapsDIR  = '/home/ss6238/hillman_servers/enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/step3_correlationMaps/';    % correlation maps generated in step3
%}

%
behaviorDIR = ' ';     % behavior data
boutsDIR = ' ';        % running blocks created in step1
stateDIR = ' ';        % states computed in step2 
corrMapsDIR  = ' ';    % correlation maps generated in step3
%}

mice = ["cm124","cm125","cm126","cm127","cm128"];

% 10-s window (200 frames)
windowlength = 200; % in frames

stateMode = 'onset';      %{'onset','offset','running'};
mod = 'jRGECO';           %{'jRGECO','HbT'}

%% estimate NNLS coefficients
%load states
% (1) running, (2) offset, (3) initial rest,
% (4) sustained rest, (5) onset
load(strcat(stateDIR,'states_',mod))

% load behavior
load(strcat(behaviorDIR,'behavior2.mat'))


for i = 1:length(mice)
    mousename = mice(i);
    
    % load running blocks
    load(strcat(boutsDIR,mousename,'.mat'))
    
    
    % filter running blocks to select runs that meet certain specifications
    % remove running blocks with less than 20-s duration
    duration = [runningBlocks.duration];
    durationIDX = duration < 400;
    runningBlocks(durationIDX) = [];
    
    % remove running blocks with pre-running rest time less than 60 s
    time2previousrun = [runningBlocks.time2previousrun];
    time2previousIDX = time2previousrun < 1200;
    runningBlocks(time2previousIDX) = [];
    
    
    % states
    states = finalStates(i).states;
    
    % flatten each state
    for j = 1:size(states,3)
        mytemp = squeeze(states(:,:,j));
        C(:,j) = mytemp(:);
        clear mytemp
    end
    
    sessionNameInit = ' ';
    count = 0;
    for run = 1:length(runningBlocks)
        fprintf('Processing run %d ...\n',run)
        sessionName = runningBlocks(run).session;
        onset = runningBlocks(run).onset;
        
        if ~strcmp(sessionNameInit,sessionName)
            clear whisk pupil correlations whisking pupilRadius
            
            behavior = getfield(B2,sessionName(1:12));
            
            %pupil, whisking, and locomotion data
            pupil = behavior.pupil;
            whisk = behavior.whisk;
            rotf = behavior.rotf;
            
            % load FC maps
            load(strcat(corrMapsDIR,sessionName(1:12)),mod)
            whisking = whisk;
            pupilRadius = pupil;
        end
        
        sessionNameInit = sessionName;
        if strcmp(mod,'jRGECO')
            FCmaps = jRGECO;
        elseif strcmp(mod,'HbT')
            FCmaps = HbT;
        end
        
        kk = 0;
        % estimate state coefficients using non-negative least squares
        a = -1200; b = 400;
        for k = a:b
            kk = kk + 1;
            
            
            whisking_ave(kk,1) =  mean(whisking(onset+k:onset+k+windowlength-1),'omitnan');
            pupil_ave(kk,1) = mean(pupilRadius(onset+k:onset+k+windowlength-1),'omitnan');
            
            
            % non-negative linear least-squares problem
            d = squeeze(FCmaps(:,:,onset+k));
            d = d(:);
            [x,resnorm,residual] = lsqnonneg(C,d);
            
            GoF(:,kk) = resnorm/norm(d.^2,1);
            X(:,kk) = x; %coefficients
            clear d x resnorm residual
        end
        
        output.coefficients = X;
        output.GoF = GoF;
        output.pupil_ave = pupil_ave;
        output.whisking_ave = whisking_ave;
        
        count = count + 1;
        results{count} = output;
        clear X GoF pupil_ave whisking_ave
        clear onset
        %}
    end
    nnlsqDistance{i} = results;
    %save(strcat(rootDIR,saveDIR,stateMode,'/',mod,'/',mice(i)),'nnlsqDistance','-v7.3')
    clear results
    
end

%% summarize results
for num = 1:length(mice)
    
    disp(num)
    %load(strcat(rootDIR,nnlsqResultsDIR,stateMode,'/',mod,'/',mice(num)))
    results = nnlsqDistance{num};
    for i = 1:length(results)
        
        output = results{i};
        
        gof_temp(:,i) = output.GoF;
        nnlsqCoefficients_temp(:,:,i) = output.coefficients;
        
    end
    %
    
    nnlsqCoefficients{num} = nnlsqCoefficients_temp;
    gof{num} = gof_temp;
    
    clear output gof_temp nnlsqCoefficients_temp

end

%% plot the average coefficients

% colors of std/se
colfill = [[0.5,0.5,0.5];[255/255,0/255,0/255];...
    [236/255,176/255,31/255];[119/255,172/255,48/255];...
    [77/255,190/255,238/255]];

%plot one for all
mycoe = [];
myres = [];
for n = 1:5
    
    % residual
    mytemp1 = gof{n};
    myres = [myres mytemp1];
    
    % nnlsq coefficients
    mytemp2 = nnlsqCoefficients{n};
    mycoe = cat(3,mycoe,mytemp2);
    
    clear mytemp1 mytemp2

end
%}

% make the last row as the first row; we do this to put onset as the first
% state
mycoe = [mycoe(5,:,:);mycoe];
mycoe(6,:,:) = [];

%
myres_mean = mean(myres,2);
myres_std = std(myres,[],2);
myres_se = myres_std / sqrt(size(myres,2));


mycoe_mean = mean(mycoe,3);
mycoe_std = std(mycoe,[],3);
mycoe_se = mycoe_std / sqrt(size(mycoe,3));

t = 1:1:size(mycoe,2);
%
figure
yyaxis right

f_ave = myres_mean';
f_sem = myres_se';

fill([t fliplr(t)], [f_ave+(f_sem) fliplr(f_ave - (f_sem ))],...
    [0.85,0.33,0.10],'FaceAlpha', 0.2,'EdgeColor','none');
hold on
h(6) = plot(t,f_ave,'LineWidth',2,'LineStyle','--');
xlabel('time (s)','fontname','Arial','FontSize',16,'fontweight','normal')
%ylabel('Residual','fontname','Arial','FontSize',16,'fontweight','bold')
box off
ylim([0 1.1])
%xlim([1 1401])
set(gca,'xtick',[1,400,800,1201],'xticklabel',...
     {'-60','-40','-20','0'},'fontname','Arial','FontSize',14) 
set(gca,'ytick',[],'yticklabel',[])
set(gca,'YColor','w')


clear f_ave f_sem

yyaxis left
for i = 1:size(mycoe,1) % for each state
    
    f_ave = mycoe_mean(i,:);
    f_sem = mycoe_se(i,:);
    
    fill([t fliplr(t)], [f_ave+(f_sem) fliplr(f_ave - (f_sem ))],...
        colfill(i,:),'FaceAlpha', 0.2,'EdgeColor','none');
    hold on
    h(i) = plot(t,f_ave,'LineStyle','-','Marker','none','color',colfill(i,:),'LineWidth',3);   
end
set(gca,'YColor','k')
set(gca,'ytick',0:0.2:1,'yticklabel',{'0','0.2','0.4','0.6','0.8','1'})
ylim([-0.8 1.1])
ylabel('nnlsq coefficients','fontname','Arial','FontSize',16,'fontweight','normal')
xline(1201,'k:','linewidth',3)
xlim([1 1401])
box off
grid on

%title(mod{k2},'FontSize',20)
%
legend(h,'c_1: onset','c_2: locomotion','c_3: offset',...
    'c_4: initial rest','c_5: sustained rest', 'residual','box','off',...
    'fontname','Arial','fontsize',14,'orientation','horizontal','fontweight','normal',...
    'Position',[0.140231643279603 0.85779411878419 0.637062473990354 0.188235291843384])

 %}
%

set(gcf,'Unit','Inches','position',[0, 0, 4*3.25,4*2],'color','w')

    

%}

% Create line: rest
annotation(gcf,'line',[0.12993762993763 0.794178794178794],...
    [0.883445945945946 0.883445945945946],...
    'Color',[0.670588235294118 0.807843137254902 0.858823529411765],...
    'LineWidth',8);

% Create textbox: rest
annotation(gcf,'textbox',...
    [0.437629937629939 0.89695945945946 0.047817047817047 0.0405405405405411],...
    'String','Rest',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);

% Create line: running
annotation(gcf,'line',[0.794178794178795 0.904365904365905],...
    [0.883445945945946 0.883445945945946],...
    'Color',[0.819607843137255 0.196078431372549 0.580392156862745],...
    'LineWidth',8);

% Create textbox: C1
annotation(gcf,'textbox',...
    [0.66735966735967 0.814189189189191 0.0291060291060283 0.0489864864864867],...
    'String','C_1',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);

% Create textbox: number of trials
annotation(gcf,'textbox',...
    [0.153846153846156 0.641891891891895 0.175675675675674 0.0489864864864868],...
    'String','n = 54 trials, 5 mice',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);

% Create textbox:  running start
annotation(gcf,'textbox',...
    [0.626819126819127 0.251689189189189 0.111226611226613 0.0489864864864865],...
    'String','running start',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'BackgroundColor',[1 1 1]);

% Create arrow
annotation(gcf,'arrow',[0.744282744282745 0.787941787941788],...
    [0.282094594594595 0.282094594594595],'LineWidth',2);


clear myR2 mycoe mycoe_mean mycoe_std mycoe_se myres myres_mean myres_std myres_se

