function quantDatHeatmap(epochDat,ccepDat,win,dtype)

intervalUse = contains(epochDat, 'Interval8');

load(epochDat);
load(ccepDat);

% MATCH CCEPs
chanList = erase(ChannelList, " ");
ccepChan = [];
for s = 1:length(CCEPMetrics)
    ccepChan = [ccepChan; CCEPMetrics{s}.stimCCEPChan];
end
% MATCH THE STIM CHAN
matchStimIdx = match_str(ccepChan, erase(ChannelList{stimIdx}," "));
ccepMatched = CCEPMetrics{matchStimIdx};
% MATCH THE TBS CHANNELS
matchChanIdx = match_str(erase(ChannelList," "), erase(ccepMatched.chanList," "));
% GET MATCHED CCEP DATA
% P-VALUE
crp_p = nan(length(ChannelList),1); 
crp_p(matchChanIdx) = ccepMatched.crp(matchChanIdx,1);
ccepChanNotMatched = find(isnan(crp_p));
% AMPLITUDE
crp_a = nan(length(ChannelList),1); 
crp_a(matchChanIdx) = ccepMatched.crp(matchChanIdx,3);
% DURATION
crp_l = nan(length(ChannelList),1); 
crp_l(matchChanIdx) = ccepMatched.crp(matchChanIdx,2);

% Get quant window idx
if intervalUse
    baselineWin = [2.4 2.5];
else
    baselineWin = [-0.9 -0.8];
end

baselineWinIdx = epochTime > baselineWin(1) & epochTime <= baselineWin(2);
quantWinIdx = epochTime > win(1) & epochTime <= win(2);

numTrial = size(epochDat,4);
numBurst = size(epochDat,1);
numChan = size(epochDat,2);

% QUANT EVOKED RESPONSE OVEDRALL
epochDatRegroup = permute(epochDat,[2,3,1,4]);
epochDatRegroup = reshape(epochDatRegroup,numChan,size(epochDat,3),[]);
evokedResponseP = zeros(numChan,1);
evokedAmplitude = zeros(numChan,1);
for c = 1:numChan
    quantDat = squeeze(peak2peak(epochDatRegroup(c,quantWinIdx,:),2));
    baselineDat = squeeze(peak2peak(epochDatRegroup(c,baselineWinIdx,:),2));
    [~,evokedResponseP(c)] = ttest2(rmoutliers(quantDat), rmoutliers(baselineDat));
    evokedAmplitude(c) = nanmean(rmoutliers(quantDat));
end
[~,~,evokedResponseAdjP] = fdr(evokedResponseP);
evokedSig = double(evokedResponseAdjP<=0.05);
evokedSig(noisyChan) = 0;

% CCEP BASED TBS RESPONSES
% CALCULATE TBS ONLY IN CCEP+ CHANNELS
ccepSigByCRP = double(crp_p < 0.05);
ccepSigByCRP(noisyChan) = 0;
[~,~,evokedResponseAdjP] = fdr(evokedResponseP(ccepSigByCRP==1));
evokedSigByCCEP = zeros(length(evokedSig), 1);
evokedSigByCCEP(ccepSigByCRP==1) = evokedResponseAdjP <= 0.05;

% QUANT TRAIN
heatmapQuant = squeeze(peak2peak(epochDat(:,:,quantWinIdx,:),3));
heatmapQuant = permute(heatmapQuant,[2,3,1]);
heatmapStat = zeros(numChan,2,4);
heatmapStatRepeated = zeros(numChan,2,5);
quantSmooth = movmean(heatmapQuant,2,2);
quantSmooth = movmean(quantSmooth,2,3);

for c = 1:numChan
    heatmapQuantChan = squeeze(heatmapQuant(c,:,:));
    [p,Table,stats] = anova2(heatmapQuantChan,numRep,'off');
    %[results,~,~,gnames] = multcompare(stats,"Estimate",'row');
    heatmapStat(c,1,1) = Table{2,5}; % WITHIN BURST
    heatmapStat(c,2,1) = Table{3,5}; % ACROSS BURST
    heatmapStat(c,1,2) = Table{2,6}; % WITHIN BURST
    heatmapStat(c,2,2) = Table{3,6}; % ACROSS BURST
    heatmapStat(c,1,3) = Table{2,2} / Table{5,2};%WITHIN partial eta2
    heatmapStat(c,2,3) = Table{3,2} / Table{5,2}; %ACROSS BURST eta2

    % WITHIN DESIGN 
    factorData = array2table(heatmapQuantChan);
    factorVar = cell(1,numBurst);
    for n = 1:numBurst
        factorVar{n} = ['t' num2str(n)];
    end
    try
        factorData.Properties.VariableNames = {'t1','t2','t3','t4','t5'};
    catch
        factorData.Properties.VariableNames = {'t1','t2','t3','t4','t5', 't6','t7','t8'};
    end
    factorData.Trials = [1:numTrial]';
    withinDsgn = table((1:numBurst)','VariableNames',{'Time'});
    rm = fitrm(factorData,['t' num2str(1) '-' 't' num2str(numBurst) '~Trials'], 'WithinDesign',withinDsgn);
    [ranovatbl] = ranova(rm);

    heatmapStatRepeated(c,1,1) = ranovatbl.F(1); % WITHIN BURST
    heatmapStatRepeated(c,1,2) = ranovatbl.pValue(1); % WITHIN BURST
    heatmapStatRepeated(c,1,3) = ranovatbl.SumSq(1) / (ranovatbl.SumSq(1) + ranovatbl.SumSq(3));%WITHIN partial eta2

    factorData = array2table(heatmapQuantChan');
    factorVar = cell(1,numTrial);
    for n = 1:numTrial
        factorVar{n} = ['t' num2str(n)];
    end
    try 
        factorData.Properties.VariableNames = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10'};
    catch
        try
            factorData.Properties.VariableNames = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10', 't11','t12'};
        catch
            try
                factorData.Properties.VariableNames = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10', 't11'};
            catch
                factorData.Properties.VariableNames = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10',...
                    't11','t12','t13','t14','t15','t16','t17','t18','t19','t20'};
            end
        end
    end
    factorData.Trials = [1:numBurst]';
    withinDsgn = table((1:numTrial)','VariableNames',{'Time'});
    rm = fitrm(factorData,['t' num2str(1) '-' 't' num2str(numTrial) '~Trials'], 'WithinDesign',withinDsgn);
    [ranovatbl] = ranova(rm);

    heatmapStatRepeated(c,2,1) = ranovatbl.F(1); % ACROSS BURST
    heatmapStatRepeated(c,2,2) = ranovatbl.pValue(1); % ACROSS BURST
    heatmapStatRepeated(c,2,3) = ranovatbl.SumSq(1) / (ranovatbl.SumSq(1) + ranovatbl.SumSq(3));
end

% null distribution
nboot = 1000;
heatmapStatRepeatedNull = zeros(numChan,2,nboot);
for c = 1:numChan
    for b = 1:nboot
        heatmapQuantChan = squeeze(heatmapQuant(c,:,:));
        shuffleTempDat = reshape(heatmapQuantChan, [], 1);
        shuffleTempDat = shuffleTempDat(randperm(length(shuffleTempDat)));
        heatmapQuantShuffled = reshape(shuffleTempDat.', size(heatmapQuantChan,2),size(heatmapQuantChan,1))';
    
        % WITHIN DESIGN 
        factorData = array2table(heatmapQuantShuffled);
        factorVar = cell(1,numBurst);
        for n = 1:numBurst
            factorVar{n} = ['t' num2str(n)];
        end
        factorData.Properties.VariableNames = factorVar;
        factorData.Trials = [1:numTrial]';
        withinDsgn = table((1:numBurst)','VariableNames',{'Time'});
        rm = fitrm(factorData,['t' num2str(1) '-' 't' num2str(numBurst) '~Trials'], 'WithinDesign',withinDsgn);
        [ranovatbl] = ranova(rm);
    
        heatmapStatRepeatedNull(c,1,b) = ranovatbl.F(1); % WITHIN BURST
    
        factorData = array2table(heatmapQuantShuffled');
        factorVar = cell(1,numTrial);
        for n = 1:numTrial
            factorVar{n} = ['t' num2str(n)];
        end
        factorData.Properties.VariableNames = factorVar;
        factorData.Trials = [1:numBurst]';
        withinDsgn = table((1:numTrial)','VariableNames',{'Time'});
        rm = fitrm(factorData,['t' num2str(1) '-' 't' num2str(numTrial) '~Trials'], 'WithinDesign',withinDsgn);
        [ranovatbl] = ranova(rm);
    
        heatmapStatRepeatedNull(c,2,b) = ranovatbl.F(1); % ACROSS BURST
    end
end

% get p-value from null dist
nullPVal = zeros(numChan, 2);
for c = 1:numChan
    nullWithin = squeeze(heatmapStatRepeatedNull(c,1,:)); % within burst
    pd = fitdist(nullWithin, 'Gamma');
    pval_temp = cdf(pd, heatmapStatRepeated(c,1,1));
    if pval_temp > 0.5
        nullPVal(c,1) = 1-pval_temp;
    else
        nullPVal(c,1) = pval_temp;
    end

    nullAcross = squeeze(heatmapStatRepeatedNull(c,2,:)); 
    pd = fitdist(nullAcross, 'Gamma');
    pval_temp = cdf(pd, heatmapStatRepeated(c,2,1));
    if pval_temp > 0.5
        nullPVal(c,2) = 1-pval_temp;
    else
        nullPVal(c,2) = pval_temp;
    end
end

% set noisy channel to null
heatmapStat(unique([stimIdx; noisyChan]),:,:) = NaN;
evokedSig(unique([stimIdx; noisyChan])) = NaN;

% FDR CORRECTION FOR F-Pvals
[~,~,withinAdjP] = fdr(squeeze(heatmapStat(evokedSig==1,1,2)));
[~,~,acrossAdjP] = fdr(squeeze(heatmapStat(evokedSig==1,2,2)));
heatmapStat(:,1,4) = heatmapStat(:,1,2);
heatmapStat(evokedSig==1,1,4) = withinAdjP;
heatmapStat(:,2,4) = heatmapStat(:,2,2);
heatmapStat(evokedSig==1,2,4) = acrossAdjP;

[~,~,withinAdjP] = fdr(squeeze(heatmapStatRepeated(evokedSig==1,1,2)));
[~,~,acrossAdjP] = fdr(squeeze(heatmapStatRepeated(evokedSig==1,2,2)));
heatmapStatRepeated(:,1,4) = heatmapStatRepeated(:,1,2);
heatmapStatRepeated(evokedSig==1,1,4) = withinAdjP;
heatmapStatRepeated(:,2,4) = heatmapStatRepeated(:,2,2);
heatmapStatRepeated(evokedSig==1,2,4) = acrossAdjP;

% 5th dimension is p-val calculated using null distribution approach
[~,~,withinAdjP] = fdr(squeeze(nullPVal(evokedSig==1,1)));
[~,~,acrossAdjP] = fdr(squeeze(nullPVal(evokedSig==1,2)));
heatmapStatRepeated(:,1,5) = nullPVal(:,1);
heatmapStatRepeated(evokedSig==1,1,5) = withinAdjP;
heatmapStatRepeated(:,2,5) = nullPVal(:,2);
heatmapStatRepeated(evokedSig==1,2,5) = acrossAdjP;

%% SPECTRAL QUANTIFICATION
epochBaseline = epochDatRegroup(:,baselineWinIdx,:);
epochQuant = epochDatRegroup(:,quantWinIdx,:);
numTotalTrials = size(epochQuant,3);

fs = 2000;
nw = 4;
freqQuant = [25 50];
specQuant = zeros(numChan, numTotalTrials,size(freqQuant,1));
specBaseline = zeros(numChan, numTotalTrials,size(freqQuant,1));

for c = 1:numChan
    for t = 1:numTotalTrials
        for f = 1:size(freqQuant,1)
            sampleDat = squeeze(epochQuant(c,:,t));
            [pxx, freq] = pmtm(sampleDat,nw,[],fs);
            freqIdx = freq>= freqQuant(f,1) & freq <= freqQuant(f,2);
            specQuant(c,t,f) = nanmean(pxx(freqIdx));

            sampleDat = squeeze(epochBaseline(c,:,t));
            [pxx, freq] = pmtm(sampleDat,nw,[],fs);
            freqIdx = freq>= freqQuant(f,1) & freq <= freqQuant(f,2);
            specBaseline(c,t,f) = nanmean(pxx(freqIdx));
        end
    end
end

evokedSpecP = zeros(numChan,size(freqQuant,1));
evokedSpecAmp = zeros(numChan,size(freqQuant,1));

for c = 1:numChan
    for f = 1:size(freqQuant,1)
        [~,evokedSpecP(c,f)] = ttest2(rmoutliers(specQuant(c,:,f)), rmoutliers(specBaseline(c,:,f)));
        evokedSpecAmp(c,f) = nanmean(specQuant(c,:,f));
    end
end

evokedSigSpec = zeros(numChan,size(freqQuant,1)); 
for f = 1:size(freqQuant,1)
    [~,~,evokedResponseAdjP] = fdr(evokedSpecP(:,f));
    evokedSigSpec(:,f) = double(evokedResponseAdjP<=0.05);
    evokedSigSpec(noisyChan,f) = 0;
end

%% DYNAMIC SPECTRAL TRAIN QUANT NEED TO WORK ON
specQuantHeatmap = zeros(numChan, numTrial,numBurst);
for c= 1:numChan
    for t=1:numTrial
        for b = 1:numBurst
            quantTemp = squeeze(epochDat(b,c,quantWinIdx,t));
            [pxx, freq] = pmtm(quantTemp,nw,[],fs);
            freqIdx = freq>= freqQuant(f,1) & freq <= freqQuant(f,2);
            specQuantHeatmap(c,t,b) = nanmean(pxx(freqIdx));
        end
    end
end

heatmapQuant_spec = specQuantHeatmap;
heatmapStatRepeated_spec = zeros(numChan,2,4);

for c = 1:numChan
    heatmapQuantChan = squeeze(heatmapQuant_spec(c,:,:));

    % WITHIN DESIGN 
    factorData = array2table(heatmapQuantChan);
    factorVar = cell(1,numBurst);
    for n = 1:numBurst
        factorVar{n} = ['t' num2str(n)];
    end
    factorData.Properties.VariableNames = factorVar;
    factorData.Trials = [1:numTrial]';
    withinDsgn = table((1:numBurst)','VariableNames',{'Time'});
    rm = fitrm(factorData,['t' num2str(1) '-' 't' num2str(numBurst) '~Trials'], 'WithinDesign',withinDsgn);
    [ranovatbl] = ranova(rm);

    heatmapStatRepeated_spec(c,1,1) = ranovatbl.F(1); % WITHIN BURST
    heatmapStatRepeated_spec(c,1,2) = ranovatbl.pValue(1); % WITHIN BURST
    heatmapStatRepeated_spec(c,1,3) = ranovatbl.SumSq(1) / (ranovatbl.SumSq(1) + ranovatbl.SumSq(3));%WITHIN partial eta2

    factorData = array2table(heatmapQuantChan');
    factorVar = cell(1,numTrial);
    for n = 1:numTrial
        factorVar{n} = ['t' num2str(n)];
    end
    factorData.Properties.VariableNames = factorVar;
    factorData.Trials = [1:numBurst]';
    withinDsgn = table((1:numTrial)','VariableNames',{'Time'});
    rm = fitrm(factorData,['t' num2str(1) '-' 't' num2str(numTrial) '~Trials'], 'WithinDesign',withinDsgn);
    [ranovatbl] = ranova(rm);

    heatmapStatRepeated_spec(c,2,1) = ranovatbl.F(1); % ACROSS BURST
    heatmapStatRepeated_spec(c,2,2) = ranovatbl.pValue(1); % ACROSS BURST
    heatmapStatRepeated_spec(c,2,3) = ranovatbl.SumSq(1) / (ranovatbl.SumSq(1) + ranovatbl.SumSq(3));
end

% set noisy channel to null
heatmapStatRepeated_spec(unique([stimIdx; noisyChan]),:,:) = NaN;

% FDR CORRECTION FOR F-Pvals
[~,~,withinAdjP] = fdr(squeeze(heatmapStatRepeated_spec(evokedSigSpec==1,1,2)));
[~,~,acrossAdjP] = fdr(squeeze(heatmapStatRepeated_spec(evokedSigSpec==1,2,2)));
heatmapStatRepeated_spec(:,1,4) = heatmapStatRepeated_spec(:,1,2);
heatmapStatRepeated_spec(evokedSigSpec==1,1,4) = withinAdjP;
heatmapStatRepeated_spec(:,2,4) = heatmapStatRepeated_spec(:,2,2);
heatmapStatRepeated_spec(evokedSigSpec==1,2,4) = acrossAdjP;

%% Save
outdata.evokedSig = evokedSig;
outdata.evokedSigSpec = evokedSigSpec;
outdata.specQuantHeatmap = specQuantHeatmap;
outdata.heatmapStat = heatmapStat;
outdata.heatmapStatRepeated = heatmapStatRepeated;
outdata.heatmapStatRepeated_spec = heatmapStatRepeated_spec;
outdata.heatmapQuant = heatmapQuant;
outdata.quantSmooth = quantSmooth;
outdata.evokedAmp = evokedAmplitude;
outdata.evokedSigByCCEP = evokedSigByCCEP;
outdata.ccepSigByCRP = ccepSigByCRP;
outdata.ccepChanNotMatched = ccepChanNotMatched;
outdata.crpa = crp_a;
outdata.crpl = crp_l;
outdata.crpp = crp_p;
save(['quant_' dtype], '-struct', 'outdata', '-v7.3');
