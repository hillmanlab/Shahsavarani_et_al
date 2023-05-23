% This script calculates the NNLS coefficients for an entire example 
% recording session using 10-second moving windows, similar to Figure 3C. 
% The correlation map within each 10-second window is reconstructed as a 
% linear combination of the behavioral states. To run this script, you need 
% to set the modality using the variable "mod". Additionally, the states 
% computed in step2 and the correlation maps generated in step3 are 
% prerequisites for this script.
%
% Author: Somayeh "Bahar" Shahsavarani
% email: bahar@huskers.unl.edu
%% initialize the directories
clear;clc;

%{
behaviorDIR = '/home/ss6238/hillman_servers/chronos/2/cmdata_analysis/RS_analysis/';     % behavior data
stateDIR = '/home/ss6238/hillman_servers/enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/step1_states/';         % states computed in step2 
corrMapsDIR  = '/home/ss6238/hillman_servers/enterprise/3/Bahar/states_025/24-codes4share/Analysis/Sanity_Testing/step3_correlationMaps/';    % correlation maps generated in step3
%}

%
behaviorDIR = '  ';     % behavior data
stateDIR = ' ';         % states computed in the previous step2 
corrMapsDIR  = '  ';    % correlation data computed in step3
%}

mice = ["cm124","cm125","cm126","cm127","cm128"];

% 10-s window (200 frames)
windowlength = 200; % in frames

% choose modality
modality ={'jRGECO','HbT'};
mod = modality{1};

%% estimate NNLS coefficients
%load states
% (1) running, (2) offset, (3) initial rest,
% (4) sustained rest, (5) onset
load(strcat(stateDIR,'states_',mod))

% load behavior
load(strcat(behaviorDIR,'behavior2.mat'))

% by num, choose different mice; by r, choose different recording sessions
for num=2%:length(mice)
    states = finalStates(num).states;
    
    corrs = dir(strcat(corrMapsDIR,mice(num),'*.mat'));
    for r = 1%:length(corrs)
        sessionName = corrs(r).name;
        behavior = getfield(B2,sessionName(1:12));
        
        %pupil, whisking, and locomotion data
        pupil = behavior.pupil;
        whisk = behavior.whisk;
        rotf = behavior.rotf;
        
        % load correlations
        load(strcat(corrMapsDIR,sessionName),mod)
        
        if strcmp(mod,'jRGECO')
            FCmaps = jRGECO;
        elseif strcmp(mod,'HbT')
            FCmaps = HbT;
        end
        
        % flatten each state
        for i = 1:size(states,3)
            mytemp = squeeze(states(:,:,i));
            C(:,i) = mytemp(:);
            clear mytemp
        end
        
        % construct every correlation map as a linear combination of the
        % states using an NNLS fit
        for j = 1:size(FCmaps,3)
            
            % non-negative linear least-squares problem
            % flatten the correlation map
            d = squeeze(FCmaps(:,:,j));
            d = d(:);
            [x,resnorm,residual] = lsqnonneg(C,d);
            
            D(:,j) = d; %actual FC
            Y(:,j) = C*x; %estimated FC
            R_norm(j,1) = resnorm;
            R(:,j) = residual; %residual d-C*x
            GoF(:,j) = resnorm/norm(d.^2,1);
            X(:,j) = x; %coefficients
            clear d x resnorm residua
        end
        
        % average behavior data within the windows
        for i = 1:size(FCmaps,3)
            rotf_ave(i) = mean(rotf(i+[0:windowlength-1]),'omitnan');
            whisk_ave(i) = mean(whisk(i+[0:windowlength-1]),'omitnan');
            pupil_ave(i) = mean(pupil(i+[0:windowlength-1]),'omitnan');
        end
    end
end

%
%% plot coefficients aligned with behavior

lineW = 2;
cols = {'#FF0000','#ECB01F','#77AC30','#4DBEEE','#808080',};

figure

% behavior
subplot(8,1,1)
plot(rotf/max(rotf),'linewidth',lineW,'color',[0.3,0.75,0.93]); axis tight
hold on
plot(whisk/max(whisk),'linewidth',lineW,'color',[0.93,0.69,0.13]); axis tight
plot(pupil/max(pupil),'linewidth',lineW,'color','k'); axis tight
leg = legend('running','whisking','pupil');
set(leg,...
    'FontSize',16,...
    'FontName','Arial','Orientation','horizontal','box','off','fontweight','normal',...
    'Position',[0.462996034806811 0.933185185442935 0.461607139796363 0.0253333330755835]);
ax = gca;
ax.XAxis.Visible = 'off';
ax.FontSize = 12;
set(ax,'YTick',[0 1],'fontname','Aria')
ylabel({'normalized','amplitude'},'fontname','Arial','fontweight','normal','fontsize',12)
box off

% average behavior
subplot(8,1,2)
plot(rotf_ave/max(rotf_ave),'linewidth',lineW,'color',[0.3,0.75,0.93]); axis tight
hold on
plot(whisk_ave/max(whisk_ave),'linewidth',lineW,'color',[0.93,0.69,0.13]); axis tight
plot(pupil_ave/max(pupil_ave),'linewidth',lineW,'color','k'); axis tight
leg = legend('ave running','ave whisking','ave pupil');
set(leg,...
    'FontSize',16,...
    'FontName','Arial','Orientation','horizontal','box','off','fontweight','normal',...
    'Position',[0.447639769450543 0.824799984158992 0.461607139796362 0.0356472791870313]);
ax = gca;
ax.XAxis.Visible = 'off';
ax.FontSize = 12;
set(ax,'YTick',[0 1],'fontname','Aria')
ylabel({'normalized','amplitude'},'fontname','Arial','fontweight','normal','fontsize',12)
box off


% locomotion onset
subplot(8,1,3)
color = cols{5};
plot(X(5,:),'linewidth',lineW,'color',sscanf(color(2:end),'%2x%2x%2x',[1 3])/255)
xlim([0 size(FCmaps,3)])
ylim([0,1.3])
legend('Locomotion onset','box','off','fontsize',16,'fontweight','normal','fontname','Arial',...
    'position',[0.427719337051601 0.690231509903939 0.155407709227951 0.0356472791870311])
ylabel('c1','fontname','Arial','fontweight','normal','fontsize',12)
ax = gca;
ax.XAxis.Visible = 'off';
set(ax,'YTick',[0 1],'fontname','Aria','fontsize',12)
box off

% locomotion
subplot(8,1,4)
color = cols{1};
plot(X(1,:),'linewidth',lineW,'color',sscanf(color(2:end),'%2x%2x%2x',[1 3])/255)
xlim([0 size(FCmaps,3)])
ylim([0,1.3])
legend('Locomotion','box','off','fontsize',16,'fontweight','normal','fontname','Arial',...
    'Position',[0.401621740466869 0.59613789095286 0.221805188693086 0.035547146380326])
ylabel('c2','fontname','Arial','fontweight','normal','fontsize',12)
ax = gca;
ax.XAxis.Visible = 'off';
set(ax,'YTick',[0 1],'fontname','Aria','fontsize',12)
box off

% locomotion offset
subplot(8,1,5)
color = cols{2};
plot(X(2,:),'linewidth',lineW,'color',sscanf(color(2:end),'%2x%2x%2x',[1 3])/255)

xlim([0 size(FCmaps,3)])
ylim([0,1.3])
legend('Locomotion offset','box','off','fontsize',16,'fontweight','normal','fontname','Arial',...
    'Position',[0.442597528448644 0.473819826499965 0.221805188693086 0.0355471463803261])
ylabel('c3','fontname','Arial','fontweight','normal','fontsize',12)
ax = gca;
ax.XAxis.Visible = 'off';
set(ax,'YTick',[0 1],'fontname','Aria','fontsize',12)
box off

% initial rest
subplot(8,1,6)
color = cols{3};
plot(X(3,:),'linewidth',lineW,'color',sscanf(color(2:end),'%2x%2x%2x',[1 3])/255)
xlim([0 size(FCmaps,3)])
ylim([0,1.3])
legend('initial rest','box','off','fontsize',16,'fontweight','normal','fontname','Arial',...
    'Position',[0.446692532543649 0.385081087761224 0.221805188693085 0.0355471463803259])
ylabel('c4','fontname','Arial','fontweight','normal','fontsize',12)
ax = gca;
ax.XAxis.Visible = 'off';
set(ax,'YTick',[0 1],'fontname','Aria','fontsize',12)
box off

% sustained rest
subplot(8,1,7)
color = cols{4};
plot(X(4,:),'linewidth',lineW,'color',sscanf(color(2:end),'%2x%2x%2x',[1 3])/255)
xlim([0 size(FCmaps,3)])
ylim([0,1.3])
legend('sustained rest','box','off','fontsize',16,'fontweight','normal','fontname','Arial',...
    'position',[0.534811283631174 0.284195690908382 0.133420295736499 0.035547146380326])
ylabel('c5','fontname','Arial','fontweight','normal','fontsize',12)
ax = gca;
ax.XAxis.Visible = 'off';
set(ax,'YTick',[0 1],'fontname','Aria','fontsize',12)
box off

% residual
subplot(8,1,8)
plot(GoF,'linewidth',lineW,'color','k')
xlim([0 size(FCmaps,3)])
ylim([0,1])
legend('nnlsq residual','box','off','fontsize',16,'fontweight','normal','fontname','Arial',...
    'Position',[0.57429578118481 0.173571514730287 0.131866256102478 0.0356472791870312])
xlabel('Frames','fontname','Arial','fontweight','normal','fontsize',16)
ylabel('res','fontname','Arial','fontweight','normal','fontsize',12)
set(gca,'YTick',[0 1],'fontname','Aria','fontsize',12)
box off

set(gcf,'Unit','Inches','position',[0,0, 3.5*6.8,3.5*2.7],'color','w')

