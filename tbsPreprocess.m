function tbsPreprocess(subjDir, subjFile,subjFilename)

load(subjFile);

% GET DAT VAR
vars = who();
temp = contains(vars, 'Trialed') | contains(vars, 'Voltage');
datString = vars{temp};
datLoad = eval(datString);
clear Trialed* *Voltage*

ParcellationValues(isnan(ParcellationValues)) = 49; % FOR SUBJ1 incomplete labels

% PARCELLATION LABEL
parcLabel = cell(1,length(ChannelList));
for c = 1:length(ChannelList)
    parcLabel{c} = TargetLabelsParc{ParcellationValues(c,8)};
end

coordElec = ParcellationValues(:,5:7);

% LOAD MNI DATA
mniDir = fullfile(subjDir, 'MNIMappedLocations');
mniFiles=dir([mniDir,'\MNICoordTrialed_Sess*']);
mniNames = cell(length(mniFiles),1);
matchIdx = 0;
for n = 1:length(mniNames)
    strtemp = strsplit(subjFilename, 'Trialed');
    if contains(erase(erase(mniFiles(n).name,'0'),'_'),erase(erase(strtemp{2},'0'),'_'))
        matchIdx = n; break;
    end
end
mniElec = load(fullfile(mniDir, mniFiles(matchIdx).name));
parcLabelVol = mniElec.VolumetricMapping.ParcellationVolLabel;

% GET WHITE MATTER INDEX
% 2 is probability elec is in the label
% 3 is probability elec is in WM
% 4 is the log ratio between the two
% mni elec in grey matter volumetric style
wmMetrics = [ParcellationValues(:,2:4) mniElec.VolumetricMapping.Volume_RatioGreyElec];


% GET STIM CHANNEL
stimSession = 1;
splitName = strsplit(subjFilename, '_');
if length(splitName) == 5
    stimChan = {[splitName{2} num2str(str2num(splitName{3})) splitName{4} num2str(str2num(splitName{5}))]};
elseif length(splitName) == 7 % for subj1 with different sessions
    stimChan = {[splitName{4} num2str(str2num(splitName{5})) splitName{6} num2str(str2num(splitName{7}))]};
    stimSession = num2str(splitName{2}(5));
else
    stimChan = {[splitName{2} splitName{3}]};
    if ~contains(stimChan, '10') %S2 contains 0 in names
        stimChan = erase(stimChan,'0');
    end
    if contains(stimChan, '09')
        stimChan = strrep(stimChan,'09','9');
    end
end
stimIdx = match_str(erase(ChannelList', " "), [stimChan]);

% AVG DATA FOR PEAK LOC
fs = int32(1/(TIME(2)-TIME(1)));
if contains(subjFilename, 'Interval8') || contains(subjFile, '9u0z')
    avgDat = squeeze(datLoad(stimIdx,:,1));
    [pks,locs] = findpeaks(abs(diff(avgDat,2)),'MinPeakHeight',500,'MinPeakDistance',5);
    baselineIdx = double(locs(end)+fs*1.5:1:locs(end)+fs*2.5); % cut very short
else
    avgDat = nanmean(squeeze(nanmean(datLoad(:,:,:),3)));
    [pks,locs] = findpeaks(abs(diff(avgDat,2)),'MinPeakHeight',1,'MinPeakDistance',5);
    baselineIdx = double(locs(1)-fs:1:locs(1)-0.1*fs); % -1s to -100ms
end

% BASELINE CORRECT TO PRIOR THE START OF THE BURST TRIAL
datLoad = permute(datLoad,[1,3,2]);
DatBaselineCor = (datLoad - squeeze(nanmean(datLoad(:,:,baselineIdx),3))) ./ std(datLoad(:,:,baselineIdx),[],3);
datLoad = permute(DatBaselineCor,[1,3,2]);

%GET NOISY CHANNEL BASED ON STD >10
datLoadAvgTrial = nanmean(datLoad,3);
noisyChan = find(std(datLoadAvgTrial(:,locs(end):locs(end)+fs),[],2) >=5);

% EPOCH DATA BASED ON BURST
burstLoc = [locs(diff(locs)>50) locs(end)]; % 4 BURST + 1 LAST BURST
if contains(subjFilename, 'Interval8')
    epochDat = zeros(length(burstLoc), size(datLoad,1), int32(length(burstLoc(1)-130:burstLoc(1)+2.5*fs-1)), size(datLoad,3));
    epochTime = -130/(double(fs)):1/double(fs):2.5-0.0005;
else
    epochDat = zeros(length(burstLoc), size(datLoad,1), int32(2*fs), size(datLoad,3));
    epochTime = -1:1/double(fs):1-0.0005;
end

for b = 1:length(burstLoc)
    if contains(subjFilename, 'Interval8')
        tempDat = datLoad(:,burstLoc(b)-130:burstLoc(b)+2.5*fs-1,:);
    else
        tempDat = datLoad(:,burstLoc(b)-fs:burstLoc(b)+fs-1,:);
    end
    epochDat(b,:,:,:) =  tempDat;
end

% GET CURRENT
splitName = strsplit(subjFilename, 'mA');
stimCurrent = str2double(splitName{1}(end));

% Save
outdata = [];
outdata.epochTime = epochTime;
outdata.epochDat = epochDat;
outdata.noisyChan = noisyChan;
outdata.stimIdx = stimIdx;
outdata.parcLabel = parcLabel;
outdata.parcLabelVol = parcLabelVol;
outdata.ChannelList = ChannelList;
outdata.stimChan = parcLabel(stimIdx);
outdata.stimChanVol = parcLabelVol(stimIdx);
outdata.stimChanLabel = ChannelList(stimIdx);
outdata.stimCurrent = stimCurrent;
outdata.coordElec = coordElec;
outdata.wmMetrics = wmMetrics;
outdata.stimSession = stimSession;
outdata.mniElec = mniElec;

save(['processedDat'], '-struct', 'outdata', '-v7.3');



