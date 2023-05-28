% This script outlines the fundamental steps for preprocessing the raw WFOM
% data. It is divided into two main parts. PART 1 is associated with 
% generating Figures 1 and S2, as well as Videos S1 and S2. PART 2 is 
% related to correlation analysis.
%% initial loader (not provided) loads raw data from camera, reading metadata and config files. Parameters for the load are stored in m (provided)
% [m, data3] = justloadWFOM(m);
% From dataset provided:
load('cm128_day8_runB_all_data.mat','info','WFOM_data_raw', 'IDX');

%% load behavior data to match raw data loaded
mouse = 'cm128';
day = '8';
run = 'runB'; 
load('behavior2.mat');
eval(sprintf('B = B2.%s_%s_%s;', mouse, day, run));

%% PART 1 --> FOR DATA USED TO GENERATE MOVIES (generates WFOM_data_processed variable from load('cm128_day8_runB_all_data.mat','info','WFOM_data_processed');


%% processing (used for videos) to remove heart rate contamination
clear corfact
ss = size(data3.green)
for i = [1:3] % step through each LED
    flkcor = squeeze(nanmean(nanmean(WFOM_data_raw.(info.acquisition.LEDs{i}),2),1)); % average each frame over space
    if strcmp(info.acquisition.LEDs{i},'red') || strcmp(info.acquisition.LEDs{i},'green')
        corfact(i,1,:) = nanmean(flkcor)+flkcor-smooth(flkcor,floor(str2double(info.acquisition.frameRate(1:2))/6)*2+1); % high pass filter (reflectance)
    else
        corfact(i,1,:) = nanmean(flkcor)+flkcor-smooth(flkcor,floor(str2double(info.acquisition.frameRate(1:2))/12)*2+1); % high pass filter (fluorescence)
    end
    data2.(info.acquisition.LEDs{i}) = nanmean(corfact(i,:,:))*data3.(info.acquisition.LEDs{i})(:,:,:)./repmat(corfact(i,:,:),[ss(1),ss(2),1]); % divide by contamination
end
disp('done')
clear WFOM_data_raw

%% Make a mask for the data by drawing around perimeter (important for denoising)
figure
imagesc(data2.lime(:,:,100)); axis image; colormap gray
mask = roipoly;

%% PCA denoise (keep a reduced number of PCs)
tlim = ss(3);
pcakeep = [300 200 200];
for i = 1:3
    eval(sprintf('tmp = reshape(data2.%s(:,:,1:tlim),[256*256,tlim]);',info.acquisition.LEDs{i}));;
    [coeff, score, latent, tsquared, explained] = pca(tmp(find(mask>0),:),'Centered','off');
    zz = zeros(size(tmp(:,1:pcakeep(i))));
    zz(find(mask>0),1:pcakeep(i)) = score(:,1:pcakeep(i));
    eval(sprintf('WFOM_data_processed.%s = reshape(zz*coeff(:,1:pcakeep(i))'',[256,256,tlim]);',info.acquisition.LEDs{i}));
    disp(['PCA... ', info.acquisition.LEDs{i}])
end
clear data2;
disp('done with PCA denoising');

%% Convert data to hemodynamic and apply hemo correction
mm.spectrafile = 'Spectra_share.mat'; % download
mm.hbcols = 'gr'; % use red and green for hemo conversion
mm.Dg = 0.4; % factors for jrgeco correction
mm.Dr = 0.8; % factors for jrgeco correction

% --- figure out baseline ----
tmp = find(B.rotf>2);
starts = tmp(find(diff(tmp)>50));
ends = tmp(find(diff(tmp)>50)+1);
for i = 1:length(starts)
    if ends(i)-starts(i)>200;
        mm.baseinterval = [starts(i):ends(i)];
        figure
        plot(info.behavior.wheelVelocity);
        hold on;
        plot(mm.baseinterval, ones(size(mm.baseinterval)),'r','LineWidth',2);
        title('red marks resting period chosen');
        legend('rotary (running)','chosen baseline');
    end
end

% ----- JRgecko correction
WFOM_data_processed.jrgeco = WFOM_data_processed.lime./((WFOM_data_processed.red.^mm.Dr).*(WFOM_data_processed.green.^mm.Dg));
WFOM_data_processed.jrgeco = WFOM_data_processed.jrgeco./repmat(mean(WFOM_data_processed.jrgeco(:,:, mm.baseinterval),3),[1,1,size(data_PCA.lime,3)])-1;

% ---- Do Hemodynamic conversion and correct jrgeco1a signal for hemo contamination
load(mm.spectrafile)
% NOTE: this file includes hemoglobin absorption spectra, as well as LED
% spectral measured on our individual WFOM system. To use this code on your
% own data, you must make your own measurements of your LED spectra.

for i = 1:2
    if strcmp(mm.hbcols(i),'g'); tmp = 'green'; nm = 530;end
    if strcmp(mm.hbcols(i),'r'); tmp = 'red'; nm = 630;end
    eval(sprintf('EHb(i) = sum((1/sum(spectra_%s))*spectra_%s.*splineyHb);', tmp, tmp));
    eval(sprintf('EHbO(i) = sum((1/sum(spectra_%s))*spectra_%s.*splineyHbO);', tmp, tmp));
    eval(sprintf('DPF(i) = sum((1/sum(spectra_%s))*spectra_%s.*splineydpff_%d);', tmp, tmp,nm));
    eval(sprintf('mua(i,:,:,:)=-(1/DPF(i))*log(squeeze(WFOM_data_processed.%s(:,:,:)./repmat(mean(WFOM_data_processed.%s(:,:,round(mm.baseinterval)),3),[1,1,size(WFOM_data_processed.%s,3)])));', tmp,tmp,tmp));
end
WFOM_data_processed.chbo = squeeze((EHb(2)*mua(1,:,:,:)-EHb(1)*mua(2,:,:,:))/(EHbO(1)*EHb(2)-EHbO(2)*EHb(1)));
WFOM_data_processed.chbr = squeeze((EHbO(2)*mua(1,:,:,:)-EHbO(1)*mua(2,:,:,:))/(EHb(1)*EHbO(2)-EHb(2)*EHbO(1)));
clear mua;



%% extract signals from ROIs (generating 92 x time H matrix)

figure; imagesc(IDX)
tmp1 = reshape(flipdim(WFOM_data_processed.jrgeco,1),[256*256,11979]);
tmp2 = reshape(flipdim(WFOM_data_processed.chbo,1),[256*256,11979]);
tmp3 = reshape(flipdim(WFOM_data_processed.chbr,1),[256*256,11979]);
tmp4 = reshape(flipdim(WFOM_data_processed.chbr+WFOM_data_processed.chbo,1),[256*256,11979]);

for i = 1:92
    H.jrgeco(i,:) = nanmean(tmp1(IDX==i,:),1);
    H.chbo(i,:) = nanmean(tmp2(IDX==i,:),1);
    H.chbr(i,:) = nanmean(tmp3(IDX==i,:),1);
    H.chbt(i,:) = nanmean(tmp4(IDX==i,:),1);
    disp(i)
end
clear tmp1 tmp2 tmp3 tmp4




%% PART 2 --> FOR DATA USED IN CORRELATION ANALYSIS

%% For correlation analysis, RAW camera underwent ROI timecourse extraction, then filtering and then conversion / correction
load('cm128_day8_runB_all_data.mat','info','WFOM_data_raw','IDX');

%%
IDX = info.rois; %IDX here has already been registered to each mouse / day
figure; imagesc(IDX)
tmp1 = reshape(flipdim(WFOM_data_rawlime,1),[256*256,11979]);
tmp2 = reshape(flipdim(WFOM_data_rawgreen,1),[256*256,11979]);
tmp3 = reshape(flipdim(WFOM_data_rawred,1),[256*256,11979]);

for i = 1:92
    H.lime(i,:) = nanmean(tmp1(IDX==i,:),1);
    H.green(i,:) = nanmean(tmp2(IDX==i,:),1);
    H.red(i,:) = nanmean(tmp3(IDX==i,:),1);
    disp(i)
end
clear tmp1 tmp2 tmp3 

%% Low pass filter
zlen = 1000; n = size(H.red,1);
H.red = [zeros(n,zlen),H.red,zeros(n,zlen)];
H.green = [zeros(n,zlen),H.green,zeros(n,zlen)];
H.lime = [zeros(n,zlen),H.lime,zeros(n,zlen)];

lpFilt = designfilt('lowpassiir','FilterOrder',20,...
'PassbandFrequency',6.5,'PassbandRipple',0.2,'SampleRate',19.9668);

for n = 1:size(H.red,1)
    H.red(n,:) = filtfilt(lpFilt,H.red(n,:));
    H.green(n,:) = filtfilt(lpFilt,H.green(n,:));
    H.lime(n,:) = filtfilt(lpFilt,H.lime(n,:));
end

H.red = H.red(:,zlen+1:end-zlen);
H.green = H.green(:,zlen+1:end-zlen);
H.lime = H.lime(:,zlen+1:end-zlen);

%% conversion and correction
mm.spectrafile = 'Spectra_share.mat'; % download
mm.hbcols = 'gr'; % use red and green for hemo conversion
mm.Dg = 0.4; % factors for jrgeco correction
mm.Dr = 0.8; % factors for jrgeco correction

% --- figure out baseline ----
tmp = find(B.rotf>2);
starts = tmp(find(diff(tmp)>50));
ends = tmp(find(diff(tmp)>50)+1);
for i = 1:length(starts)
    if ends(i)-starts(i)>200;
        mm.baseinterval = [starts(i):ends(i)];
        figure
        plot(info.behavior.wheelVelocity);
        hold on;
        plot(mm.baseinterval, ones(size(mm.baseinterval)),'r','LineWidth',2);
        title('red marks resting period chosen');
        legend('rotary (running)','chosen baseline');
    end
end

% ----- JRgecko correction
H.jrgeco = H.lime./((H.red.^mm.Dr).*(H.green.^mm.Dg));
H.jrgeco = H.jrgeco./repmat(mean(H.jrgeco(:,:, mm.baseinterval),3),[1,1,size(H.lime,3)])-1;

% ---- Do Hemodynamic conversion and correct jrgeco1a signal for hemo contamination
load(mm.spectrafile)
% NOTE: this file includes hemoglobin absorption spectra, as well as LED
% spectral measured on our individual WFOM system. To use this code on your
% own data, you must make your own measurements of your LED spectra.

for i = 1:2
    if strcmp(mm.hbcols(i),'g'); tmp = 'green'; nm = 530;end
    if strcmp(mm.hbcols(i),'r'); tmp = 'red'; nm = 630;end
    eval(sprintf('EHb(i) = sum((1/sum(spectra_%s))*spectra_%s.*splineyHb);', tmp, tmp));
    eval(sprintf('EHbO(i) = sum((1/sum(spectra_%s))*spectra_%s.*splineyHbO);', tmp, tmp));
    eval(sprintf('DPF(i) = sum((1/sum(spectra_%s))*spectra_%s.*splineydpff_%d);', tmp, tmp,nm));
    eval(sprintf('mua(i,:,:,:)=-(1/DPF(i))*log(squeeze(H.%s(:,:,:)./repmat(mean(H.%s(:,:,round(mm.baseinterval)),3),[1,1,size(H.%s,3)])));', tmp,tmp,tmp));
end
H.chbo = squeeze((EHb(2)*mua(1,:,:,:)-EHb(1)*mua(2,:,:,:))/(EHbO(1)*EHb(2)-EHbO(2)*EHb(1)));
H.chbr = squeeze((EHbO(2)*mua(1,:,:,:)-EHbO(1)*mua(2,:,:,:))/(EHb(1)*EHbO(2)-EHb(2)*EHbO(1)));
clear mua;

% the variable H generated here is equivalenet to "data" in our shared datasets (e.g. cm128_8_runB.mat)
