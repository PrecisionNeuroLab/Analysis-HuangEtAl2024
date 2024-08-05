function tbsGroupQuant(groupDat)

load(groupDat);

subjLabel = []; subjLabelMatched = [];
numTBS_low = [];
numWithin_low = [];
numAcross_low = [];
numBoth_low = [];
numNonMod_low = [];
numNonTBS_low = [];
numTBS_high = [];
numWithin_high = [];
numAcross_high = [];
numBoth_high = [];
numNonMod_high = [];
numNonTBS_high = [];
etaWithin_low = [];
etaWithin_high = [];
etaAcross_low = [];
etaAcross_high = [];
lowCurrentWithin = [];
lowCurrentAcross = [];
lowCurrentBoth = [];
highCurrentWithin = [];
highCurrentAcross = [];
highCurrentBoth = [];
lowCurrentChanMod = [];
highCurrentChanMod = [];
lowCurrentPerChanMod = [];
highCurrentPerChanMod = [];
lowCurrentChanDyna = [];
highCurrentChanDyna = [];
lowCurrentPercChanDyan = [];
highCurrentPercChanDyan = [];
lowCurrentStimIdx = [];
lowCurrentStimLabel = [];
highCurrentStimIdx = [];
highCurrentStimLabel = [];
lowCurrentEtaWithinGroup = [];
lowCurrentEtaAcrossGroup = [];
highCurrentEtaWithinGroup = [];
highCurrentEtaAcrossGroup = [];
lowCurrentNetwork = [];
highCurrentNetwork = [];
groupEC = []; groupDIST = []; groupAMP = []; groupPLV = [];
groupBW = []; groupPART = []; groupWDZ = []; groupPEA = [];
groupPlasticityPerc = []; groupEWITHIN = []; groupLAT = [];
groupWM = []; groupQexp = [];
groupStimWM = []; groupEvokedAmp = [];
groupEvokedAmpLow = []; groupEvokedAmpHigh = [];

etaWithinDLPFC_low = []; etaWithinDLPFC_high = []; tbsDLPFC_high = []; 
etaAcrossDLPFC_low = []; etaAcrossDLPFC_high = []; tbsDLPFC_low = [];
DLPFC_label_low = []; DLPFC_label_high = []; 
etaWithinVLPFC_low = []; etaWithinVLPFC_high = []; tbsVLPFC_high = [];
etaAcrossVLPFC_low = []; etaAcrossVLPFC_high = []; tbsVLPFC_low = [];
VLPFC_label_low = []; VLPFC_label_high = [];
etaWithinOFC_low = []; etaWithinOFC_high = []; tbsOFC_high = [];
etaAcrossOFC_low = []; etaAcrossOFC_high = []; tbsOFC_low = [];
OFC_label_low = []; OFC_label_high = [];
etaWithinTEMPORAL_low = []; etaWithinTEMPORAL_high = []; tbsTEMPORAL_high = [];
etaAcrossTEMPORAL_low = []; etaAcrossTEMPORAL_high = [];  tbsTEMPORAL_low = [];
TEMPORAL_label_low = []; TEMPORAL_label_high = [];
etaWithinINS_low = []; etaWithinINS_high = []; tbsINS_high = [];
etaAcrossINS_low = []; etaAcrossINS_high = [];  tbsINS_low = [];
INS_label_low = []; INS_label_high = [];
etaWithinCING_low = []; etaWithinCING_high = []; tbsCING_high = [];
etaAcrossCING_low = []; etaAcrossCING_high = [];  tbsCING_low = [];
CING_label_low = []; CING_label_high = [];
etaWithinCENTRAL_low = []; etaWithinCENTRAL_high = []; tbsCENTRAL_high = [];
etaAcrossCENTRAL_low = []; etaAcrossCENTRAL_high = [];  tbsCENTRAL_low = [];
CENTRAL_label_low = []; CENTRAL_label_high = [];
cnt_dlpfc = [0 0]; cnt_opercular = [0 0]; cnt_temporal = [0 0]; cnt_ofc = [0 0];
cnt_ins = [0 0]; cnt_cing = [0 0]; cnt_central = [0 0];

DLPFC_unilabel_low = []; DLPFC_uniMod_low = [];
VLPFC_unilabel_low = []; VLPFC_uniMod_low = [];
INS_unilabel_low = []; INS_uniMod_low = [];
CING_unilabel_low = []; CING_uniMod_low = [];
TEMPORAL_unilabel_low = []; TEMPORAL_uniMod_low = [];
CENTRAL_unilabel_low = []; CENTRAL_uniMod_low = [];
OFC_unilabel_low = []; OFC_uniMod_low = [];

DLPFC_unilabel_high = []; DLPFC_uniMod_high = [];
VLPFC_unilabel_high = []; VLPFC_uniMod_high = [];
INS_unilabel_high = []; INS_uniMod_high = [];
CING_unilabel_high = []; CING_uniMod_high = [];
TEMPORAL_unilabel_high = []; TEMPORAL_uniMod_high = [];
CENTRAL_unilabel_high = []; CENTRAL_uniMod_high = [];
OFC_unilabel_high = []; OFC_uniMod_high = [];

DLPFC_coord_low = []; DLPFC_coord_high = []; 
VLPFC_coord_low = []; VLPFC_coord_high = []; 
INS_coord_low = []; INS_coord_high = []; 
CING_coord_low = []; CING_coord_high = []; 
TEMPORAL_coord_low = []; TEMPORAL_coord_high = []; 
CENTRAL_coord_low = []; CENTRAL_coord_high = []; 
OFC_coord_low = []; OFC_coord_high = []; 

DLPFC_elecLabel_low = []; DLPFC_elecLabel_high = []; 
VLPFC_elecLabel_low = []; VLPFC_elecLabel_high = []; 
INS_elecLabel_low = []; INS_elecLabel_high = []; 
CING_elecLabel_low = []; CING_elecLabel_high = []; 
TEMPORAL_elecLabel_low = []; TEMPORAL_elecLabel_high = []; 
CENTRAL_elecLabel_low = []; CENTRAL_elecLabel_high = []; 
OFC_elecLabel_low = []; OFC_elecLabel_high = []; 

DLPFC_modStatusWithin_low = []; DLPFC_modStatusWithin_high = []; 
VLPFC_modStatusWithin_low = []; VLPFC_modStatusWithin_high = []; 
INS_modStatusWithin_low = []; INS_modStatusWithin_high = []; 
CING_modStatusWithin_low = []; CING_modStatusWithin_high = []; 
TEMPORAL_modStatusWithin_low = []; TEMPORAL_modStatusWithin_high = []; 
CENTRAL_modStatusWithin_low = []; CENTRAL_modStatusWithin_high = []; 
OFC_modStatusWithin_low = []; OFC_modStatusWithin_high = []; 

DLPFC_modStatusAcross_low = []; DLPFC_modStatusAcross_high = []; 
VLPFC_modStatusAcross_low = []; VLPFC_modStatusAcross_high = []; 
INS_modStatusAcross_low = []; INS_modStatusAcross_high = []; 
CING_modStatusAcross_low = []; CING_modStatusAcross_high = []; 
TEMPORAL_modStatusAcross_low = []; TEMPORAL_modStatusAcross_high = []; 
CENTRAL_modStatusAcross_low = []; CENTRAL_modStatusAcross_high = []; 
OFC_modStatusAcross_low = []; OFC_modStatusAcross_high = []; 

DLPFC_modStatusBoth_low = []; DLPFC_modStatusBoth_high = []; 
VLPFC_modStatusBoth_low = []; VLPFC_modStatusBoth_high = []; 
INS_modStatusBoth_low = []; INS_modStatusBoth_high = []; 
CING_modStatusBoth_low = []; CING_modStatusBoth_high = []; 
TEMPORAL_modStatusBoth_low = []; TEMPORAL_modStatusBoth_high = []; 
CENTRAL_modStatusBoth_low = []; CENTRAL_modStatusBoth_high = []; 
OFC_modStatusBoth_low = []; OFC_modStatusBoth_high = []; 
temp = 1; storeS = 0; storeI = 0;

groupTrainHeatmap = [];
groupBurstHeatmap = [];

% DECODING VARIABLES
decodeFeaturesGroup = [];
decodeDistGroup = [];
decodeOutcomesGroup = [];

% SPK VARIABLES
spkTGroup = [];
spkLogGroup = [];
spkRestingConGroup = [];
spkCCEPGroup = [];
spkTBSAmpGroup = [];
spkTBSDynamicGroup = [];
spkTBSSig= [];
spkTLow = [];
spkTHigh = [];

% demographics
subjCombined = [];

for subj = 1:length(subjGroup)
    if subj == 10 % SUBJ9 has both 8 pulses and 5 pulses data
        runSubjFiles = 2:2:length(subjGroup{subj});
    elseif subj == 1
        runSubjFiles = 1:length(subjGroup{subj}); %1:3 is session 1
    else
        runSubjFiles = 1:length(subjGroup{subj});
    end
    for i = runSubjFiles

        temp = temp+1;
        if temp == 35
            storeS = subj; storeI = i;
        end
        subjLabel = [subjLabel; subj];
        sessionNum = subjGroup{subj}{i}.demo.stimSession;
        stimCurrent = subjGroup{subj}{i}.demo.stimCurrent;
        stimIdx = subjGroup{subj}{i}.demo.stimIdx;
        chanName = subjGroup{subj}{i}.demo.ChannelList;
        pWithin = squeeze(subjGroup{subj}{i}.heatmapStatRepeated(:,1,4)); % use 2 for non-FDR correction and 4 for FDR correction
        pAcross = squeeze(subjGroup{subj}{i}.heatmapStatRepeated(:,2,4));
        crpp = subjGroup{subj}{i}.crpp;
        elecCoords = subjGroup{subj}{i}.demo.mniElec.RASRecMNIBipolar;
        elecLabel = subjGroup{subj}{i}.demo.ChannelList;

        % DEFINE GROUPS BASED ON CCEP RESPONSES
        useCCEPDef = 0;
        if useCCEPDef == 0
            TBSidx = subjGroup{subj}{i}.evokedSig==1; % DEFINE GROUPS BASED ON TBS
            noTBSidx = ~TBSidx; noTBSidx(subjGroup{subj}{i}.demo.noisyChan) = 0;
            totalChanNum = length(chanName);
        else
            TBSidx = subjGroup{subj}{i}.evokedSigByCCEP==1; % DEFINE GROUPS BASED ON CCEP
            noTBSidx = crpp<0.05 & ~TBSidx; noTBSidx(subjGroup{subj}{i}.demo.noisyChan) = 0;
            totalChanNum = length(find(crpp<0.05));
        end
        plasticityIdx = (pWithin <= 0.05 | pAcross <= 0.05) & TBSidx;
        nonplasticityIdx = (pWithin > 0.05 & pAcross > 0.05) & TBSidx;
        plasticityWithin = pWithin <= 0.05 & TBSidx;
        plasticityAcross = pAcross <= 0.05 & TBSidx;
        plastWithinExclu = pWithin <= 0.05 & pAcross > 0.05 & TBSidx;
        plastAcrossExclu = pWithin > 0.05 & pAcross <= 0.05 & TBSidx;
        plastBoth = pWithin <= 0.05 & pAcross <= 0.05 & TBSidx;
        

        % GET GROUP VALUES
        etaWithin = squeeze(subjGroup{subj}{i}.heatmapStatRepeated(:,1,3));
        etaAcross = squeeze(subjGroup{subj}{i}.heatmapStatRepeated(:,2,3));
        etaWithin_TBS = squeeze(subjGroup{subj}{i}.heatmapStatRepeated(TBSidx,1,3));
        etaAcross_TBS = squeeze(subjGroup{subj}{i}.heatmapStatRepeated(TBSidx,2,3));
        numChanModDyna = length(find(plasticityIdx));
        numChanMod = length(find(TBSidx));
        perChanModDyna = numChanModDyna / totalChanNum;
        perChanMod = numChanMod / totalChanNum;
        nonModNum = length(find(nonplasticityIdx));
        withinNum = length(find(plastWithinExclu));
        acrossNum = length(find(plastAcrossExclu));
        bothNum = length(find(plastBoth));

        withinNumPct = length(find(plastWithinExclu)) / totalChanNum;
        acrossNumPct = length(find(plastAcrossExclu)) / totalChanNum;
        bothNumPct = length(find(plastBoth)) / totalChanNum;

        % STIM WM
        %wmMetric = subjGroup{subj}{i}.demo.mniElec.VolumetricMapping.ProbabilityELAGrey;
        %wmMetric = subjGroup{subj}{i}.demo.mniElec.VolumetricMapping.VolumesWithinGreyMatter;
        %wmMetric = subjGroup{subj}{i}.demo.mniElec.VolumetricMapping.DistanceGreyWhite;
        wmMetric = log(subjGroup{subj}{i}.demo.mniElec.VolumetricMapping.ProbabilityELAGrey ./ subjGroup{subj}{i}.demo.mniElec.VolumetricMapping.ProbabilityELAWhite);
        wmStim = wmMetric(stimIdx);
        %wmStim = subjGroup{subj}{i}.demo.wmMetrics(stimIdx,4); % 3 is log ratio, 2 is WM%, 1 is label%
        groupStimWM = [groupStimWM wmStim];
        
        % PLASTICITY
        evokedAmp = subjGroup{subj}{i}.evokedAmp(((TBSidx | noTBSidx)));
        tbsStatus = TBSidx((TBSidx | noTBSidx));
        modStatusBoth = plastBoth((TBSidx | noTBSidx));
        modStatusAcross = plastAcrossExclu((TBSidx | noTBSidx));
        modStatusWithin = plastWithinExclu((TBSidx | noTBSidx));
        groupEvokedAmp = [groupEvokedAmp nanmean(evokedAmp)];

        % ETA WITHIN
        EWITHIN = zeros(1,2);
        EWITHIN(1) = nanmean(etaWithin(plasticityIdx));
        EWITHIN(2) = nanmean(etaWithin(nonplasticityIdx));
        groupEWITHIN = [groupEWITHIN; EWITHIN];

        % HEATMAP DYNAMICS QUANT
        modHeatmap = subjGroup{subj}{i}.heatmapQuant(plastWithinExclu,:,:);
        groupBurstHeatmap = cat(1, groupBurstHeatmap, modHeatmap(:,1:10,1:5));
        modHeatmap = subjGroup{subj}{i}.heatmapQuant(plastAcrossExclu,:,:);
        groupTrainHeatmap = cat(1, groupTrainHeatmap, modHeatmap(:,1:10,1:5));
        
        subjRep = ones(length(chanName),1)*stimCurrent;
        subjCombined = cat(1,subjCombined,subjRep(plastWithinExclu));

        % DISTANCE
        distanceFromSite = subjGroupRest{subj}.restPre.dist(stimIdx,:);
        DIST = zeros(1,7);
        DIST(1) = nanmean(distanceFromSite(plasticityIdx));
        DIST(2) = nanmean(distanceFromSite(nonplasticityIdx));
        DIST(3) = nanmean(distanceFromSite(plastWithinExclu));
        DIST(4) = nanmean(distanceFromSite(plastAcrossExclu));
        DIST(5) = nanmean(distanceFromSite(TBSidx));
        DIST(6) = nanmean(distanceFromSite(noTBSidx));
        DIST(7) = nanmean(distanceFromSite(plastBoth));
        groupDIST = [groupDIST; DIST];

        % CCEP DURATION
        TBSLat = subjGroup{subj}{i}.crpl;
        LAT = zeros(1,7);
        LAT(1) = nanmean(TBSLat(plasticityIdx));
        LAT(2) = nanmean(TBSLat(nonplasticityIdx));
        LAT(3) = nanmean(TBSLat(plastWithinExclu));
        LAT(4) = nanmean(TBSLat(plastAcrossExclu));
        LAT(5) = nanmean(TBSLat(TBSidx));
        LAT(6) = nanmean(TBSLat(noTBSidx));
        LAT(7) = nanmean(TBSLat(plastBoth));
        groupLAT = [groupLAT; LAT];


        % CCEP AMPLITUDE
        %TBSAmp = subjGroup{subj}{i}.evokedAmp;
        TBSAmp = hampel(subjGroup{subj}{i}.crpa,length(subjGroup{subj}{i}.crpa),80);
        AMP = zeros(1,7);
        AMP(1) = nanmean(TBSAmp(plasticityIdx));
        AMP(2) = nanmean(TBSAmp(nonplasticityIdx));
        AMP(3) = nanmean(TBSAmp(plastWithinExclu));
        AMP(4) = nanmean(TBSAmp(plastAcrossExclu));
        AMP(5) = nanmean(TBSAmp(TBSidx));
        AMP(6) = nanmean(TBSAmp(noTBSidx));
        AMP(7) = nanmean(TBSAmp(plastBoth));
        groupAMP = [groupAMP; AMP];

        % PLV
        PLV = zeros(1,7);
        freqUse = 2;
        plvUse = squeeze(subjGroupRest{subj}.restPre.epochRestingPre_plv(freqUse,stimIdx,:));
        PLV(1) = nanmean(plvUse(plasticityIdx));
        PLV(2) = nanmean(plvUse(nonplasticityIdx));
        PLV(3) = nanmean(plvUse(plastWithinExclu));
        PLV(4) = nanmean(plvUse(plastAcrossExclu));
        PLV(5) = nanmean(plvUse(TBSidx));
        PLV(6) = nanmean(plvUse(noTBSidx));
        PLV(7) = nanmean(plvUse(plastBoth));
        groupPLV = [groupPLV; PLV];

        % voltage correlations
        cohUse = squeeze(subjGroupRest{subj}.restPre.volPea(stimIdx,:))';
        PEA = zeros(1,7);
        PEA(1) = nanmean(cohUse(plasticityIdx));
        PEA(2) = nanmean(cohUse(nonplasticityIdx));
        PEA(3) = nanmean(cohUse(plastWithinExclu));
        PEA(4) = nanmean(cohUse(plastAcrossExclu));
        PEA(5) = nanmean(cohUse(TBSidx));
        PEA(6) = nanmean(cohUse(noTBSidx));
        PEA(7) = nanmean(cohUse(plastBoth));
        groupPEA = [groupPEA; PEA];

        % WM

        WM = zeros(1,7);
        wmUse = wmMetric;
        wmUse(wmUse==49) = NaN; % 49 was set prev to hide NaN for S1
        WM(1) = nanmean(wmUse(plasticityIdx));
        WM(2) = nanmean(wmUse(nonplasticityIdx));
        WM(3) = nanmean(wmUse(plastWithinExclu));
        WM(4) = nanmean(wmUse(plastAcrossExclu));
        WM(5) = nanmean(wmUse(TBSidx));
        WM(6) = nanmean(wmUse(noTBSidx));
        WM(7) = nanmean(wmUse(plastBoth));
        groupWM = [groupWM; WM];

        % group Qexp
        QEXP = zeros(1,7);
        freqUse = 2;
        QEXPUse = nanmean(subjGroupRest{subj}.restPre.qexp(freqUse,:,:),3); 
        QEXP(1) = nanmean(QEXPUse(plasticityIdx));
        QEXP(2) = nanmean(QEXPUse(nonplasticityIdx));
        QEXP(3) = nanmean(QEXPUse(plastWithinExclu));
        QEXP(4) = nanmean(QEXPUse(plastAcrossExclu));
        QEXP(5) = nanmean(QEXPUse(TBSidx));
        QEXP(6) = nanmean(QEXPUse(noTBSidx));
        QEXP(7) = nanmean(QEXPUse(plastBoth));
        groupQexp = [groupQexp; QEXP];

       

        % DECODING VARIABLES
        decodeFeaturesSubj = [];
        decodeOutcomesSubj = [];
        decodeDistSubj = [];
        decodeDistSubj = cat(2,decodeDistSubj,distanceFromSite(TBSidx | noTBSidx)');
        decodeFeaturesSubj = cat(2,decodeFeaturesSubj, (distanceFromSite(TBSidx | noTBSidx)'));
        decodeFeaturesSubj = cat(2,decodeFeaturesSubj, (TBSLat(TBSidx | noTBSidx)));
        decodeFeaturesSubj = cat(2,decodeFeaturesSubj, (TBSAmp(TBSidx | noTBSidx)));
        decodeFeaturesSubj = cat(2,decodeFeaturesSubj, (plvUse(TBSidx | noTBSidx)));
        decodeFeaturesSubj = cat(2,decodeFeaturesSubj, (cohUse(TBSidx | noTBSidx)));
        decodeFeaturesSubj = cat(2,decodeFeaturesSubj, (wmUse(TBSidx | noTBSidx)));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, TBSidx(TBSidx | noTBSidx));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, plastWithinExclu(TBSidx | noTBSidx));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, plastAcrossExclu(TBSidx | noTBSidx));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, plastBoth(TBSidx | noTBSidx));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, plasticityIdx(TBSidx | noTBSidx));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, plasticityWithin(TBSidx | noTBSidx));
        decodeOutcomesSubj = cat(2, decodeOutcomesSubj, plasticityAcross(TBSidx | noTBSidx));

        decodeDistGroup = cat(1,decodeDistGroup,decodeDistSubj);
        decodeFeaturesGroup = cat(1,decodeFeaturesGroup,decodeFeaturesSubj);
        decodeOutcomesGroup = cat(1,decodeOutcomesGroup,decodeOutcomesSubj);
        % NETWORK 
        EC = zeros(1,3);
        EC(1) = nanmean(subjGroupRest{subj}.restPre.ec(plasticityIdx,2));
        EC(2) = nanmean(subjGroupRest{subj}.restPre.ec(nonplasticityIdx,2));
        EC(3) = subjGroupRest{subj}.restPre.ec(stimIdx,2);
        groupEC = [groupEC; EC];
        BW = zeros(1,2);
        BW(1) = nanmean(subjGroupRest{subj}.restPre.bw(plasticityIdx,2));
        BW(2) = nanmean(subjGroupRest{subj}.restPre.bw(nonplasticityIdx,2));
        groupBW = [groupBW; BW];
        PART = zeros(1,2);
        PART(1) = nanmean(subjGroupRest{subj}.restPre.part(plasticityIdx,2));
        PART(2) = nanmean(subjGroupRest{subj}.restPre.part(nonplasticityIdx,2));
        groupPART = [groupPART; PART];
        WDZ = zeros(1,2);
        WDZ(1) = nanmean(subjGroupRest{subj}.restPre.wdz(plasticityIdx,2));
        WDZ(2) = nanmean(subjGroupRest{subj}.restPre.wdz(nonplasticityIdx,2));
        groupWDZ = [groupWDZ;WDZ];

        groupPlasticityPerc = [groupPlasticityPerc; perChanMod];

        % DLPFC
        DLPFC = {'rostralmiddlefrontal','caudalmiddlefrontal','superiorfrontal'};
        % OPERCULAR
        VLPFC = {'parsopercularis','parstriangularis', 'parsorbitalis'};
        % CINGULATE
        CING = {'posteriorcingulate','caudalanteriorcingulate', 'medialorbitofrontal'};
        % INSULA
        INS = {'insula', 'putamen', 'precentral'};
        % OFC
        OFC = {'lateralorbitofrontal'};
        % TEMPORAL
        TEMPORAL = {'middletemporal','superiortemporal'};
        % CENTRAL
        CENTRAL = {'postcentral'};

        subjAnatName = subjGroup{subj}{i}.demo.parcLabel((TBSidx | noTBSidx))';
        subjElecCoord = elecCoords((TBSidx | noTBSidx),:);
        subjElecLabel = elecLabel((TBSidx | noTBSidx));

        if subjGroup{subj}{i}.demo.stimCurrent == 1
            etaWithin_low = [etaWithin_low; etaWithin_TBS];
            etaAcross_low = [etaAcross_low; etaAcross_TBS];
            lowCurrentWithin = [lowCurrentWithin; withinNum / totalChanNum];
            lowCurrentAcross = [lowCurrentAcross; acrossNum / totalChanNum];
            lowCurrentBoth = [lowCurrentBoth; bothNum / totalChanNum];
            lowCurrentChanDyna = [lowCurrentChanDyna; numChanModDyna];
            lowCurrentPercChanDyan = [lowCurrentPercChanDyan; perChanModDyna];
            lowCurrentChanMod = [lowCurrentChanMod; numChanMod];
            lowCurrentPerChanMod = [lowCurrentPerChanMod; perChanMod];
            lowCurrentStimIdx = [lowCurrentStimIdx; {[num2str(subj) '-' 'sesh' num2str(sessionNum) '-' chanName{stimIdx}]}];
            lowCurrentStimLabel = [lowCurrentStimLabel; subjGroup{subj}{i}.demo.parcLabel(stimIdx)];
            lowCurrentEtaWithinGroup = [lowCurrentEtaWithinGroup; etaWithin];
            lowCurrentEtaAcrossGroup = [lowCurrentEtaAcrossGroup; etaAcross];
            lowCurrentNetwork = [lowCurrentNetwork; subjGroupRest{subj}.restPre.cw(stimIdx,2)];
            groupEvokedAmpLow = [groupEvokedAmpLow; nanmean(evokedAmp)];

            % RAW NUMBERS
            numNonMod_low = [numNonMod_low; nonModNum];
            numTBS_low = [numTBS_low; numChanMod];
            numWithin_low = [numWithin_low; withinNum];
            numAcross_low = [numAcross_low; acrossNum];
            numNonTBS_low = [numNonTBS_low; length(find(noTBSidx))];
            numBoth_low = [numBoth_low; bothNum];

            % GROUP BY STIM
            if any(match_str(DLPFC,subjGroup{subj}{i}.demo.stimChan)) 
                etaWithinDLPFC_low = [etaWithinDLPFC_low; etaWithin];
                etaAcrossDLPFC_low = [etaAcrossDLPFC_low; etaAcross];
                tbsDLPFC_low = [tbsDLPFC_low; evokedAmp];
                DLPFC_label_low = [DLPFC_label_low; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                DLPFC_unilabel_low = [DLPFC_unilabel_low; subjAnatName];
                DLPFC_uniMod_low = [DLPFC_uniMod_low; tbsStatus];
                DLPFC_coord_low = [DLPFC_coord_low; subjElecCoord];
                DLPFC_elecLabel_low = [DLPFC_elecLabel_low subjElecLabel];
                DLPFC_modStatusWithin_low = [DLPFC_modStatusWithin_low; modStatusWithin];
                DLPFC_modStatusAcross_low = [DLPFC_modStatusAcross_low; modStatusAcross];
                DLPFC_modStatusBoth_low = [DLPFC_modStatusBoth_low; modStatusBoth];

                cnt_dlpfc(1) = cnt_dlpfc(1) + 1;
            elseif any(match_str(VLPFC,subjGroup{subj}{i}.demo.stimChan))
                etaWithinVLPFC_low = [etaWithinVLPFC_low; etaWithin];
                etaAcrossVLPFC_low = [etaAcrossVLPFC_low; etaAcross];
                tbsVLPFC_low = [tbsVLPFC_low; evokedAmp];
                VLPFC_label_low = [VLPFC_label_low; subjAnatName];

                 % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                VLPFC_unilabel_low = [VLPFC_unilabel_low; subjAnatName];
                VLPFC_uniMod_low = [VLPFC_uniMod_low; tbsStatus];
                VLPFC_coord_low = [VLPFC_coord_low; subjElecCoord];
                VLPFC_elecLabel_low = [VLPFC_elecLabel_low subjElecLabel];
                VLPFC_modStatusWithin_low = [VLPFC_modStatusWithin_low; modStatusWithin];
                VLPFC_modStatusAcross_low = [VLPFC_modStatusAcross_low; modStatusAcross];
                VLPFC_modStatusBoth_low = [VLPFC_modStatusBoth_low; modStatusBoth];


                cnt_opercular(1) = cnt_opercular(1) + 1;
            elseif any(match_str(OFC,subjGroup{subj}{i}.demo.stimChan))
                etaWithinOFC_low = [etaWithinOFC_low; etaWithin];
                etaAcrossOFC_low = [etaAcrossOFC_low; etaAcross];
                tbsOFC_low = [tbsOFC_low; evokedAmp];
                OFC_label_low = [OFC_label_low; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                OFC_unilabel_low = [OFC_unilabel_low; subjAnatName];
                OFC_uniMod_low = [OFC_uniMod_low; tbsStatus];
                OFC_coord_low = [OFC_coord_low; subjElecCoord];
                OFC_elecLabel_low = [OFC_elecLabel_low subjElecLabel];
                OFC_modStatusWithin_low = [OFC_modStatusWithin_low; modStatusWithin];
                OFC_modStatusAcross_low = [OFC_modStatusAcross_low; modStatusAcross];
                OFC_modStatusBoth_low = [OFC_modStatusBoth_low; modStatusBoth];


                cnt_ofc(1) = cnt_ofc(1) + 1;
            elseif any(match_str(TEMPORAL,subjGroup{subj}{i}.demo.stimChan))
                etaWithinTEMPORAL_low = [etaWithinTEMPORAL_low; etaWithin];
                etaAcrossTEMPORAL_low = [etaAcrossTEMPORAL_low; etaAcross];
                tbsTEMPORAL_low = [tbsTEMPORAL_low; evokedAmp];
                TEMPORAL_label_low = [TEMPORAL_label_low; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                TEMPORAL_unilabel_low = [TEMPORAL_unilabel_low; subjAnatName];
                TEMPORAL_uniMod_low = [TEMPORAL_uniMod_low; tbsStatus];
                TEMPORAL_coord_low = [TEMPORAL_coord_low; subjElecCoord];
                TEMPORAL_elecLabel_low = [TEMPORAL_elecLabel_low subjElecLabel];
                TEMPORAL_modStatusWithin_low = [TEMPORAL_modStatusWithin_low; modStatusWithin];
                TEMPORAL_modStatusAcross_low = [TEMPORAL_modStatusAcross_low; modStatusAcross];
                TEMPORAL_modStatusBoth_low = [TEMPORAL_modStatusBoth_low; modStatusBoth];


                cnt_temporal(1) = cnt_temporal(1) + 1;
            elseif any(match_str(INS,subjGroup{subj}{i}.demo.stimChan))
                etaWithinINS_low = [etaWithinINS_low; etaWithin];
                etaAcrossINS_low = [etaAcrossINS_low; etaAcross];
                tbsINS_low = [tbsINS_low; evokedAmp];
                INS_label_low = [INS_label_low; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                INS_unilabel_low = [INS_unilabel_low; subjAnatName];
                INS_uniMod_low = [INS_uniMod_low; tbsStatus];
                INS_coord_low = [INS_coord_low; subjElecCoord];
                INS_elecLabel_low = [INS_elecLabel_low subjElecLabel];
                INS_modStatusWithin_low = [INS_modStatusWithin_low; modStatusWithin];
                INS_modStatusAcross_low = [INS_modStatusAcross_low; modStatusAcross];
                INS_modStatusBoth_low = [INS_modStatusBoth_low; modStatusBoth];


                cnt_ins(1) = cnt_ins(1) + 1;
            elseif any(match_str(CING,subjGroup{subj}{i}.demo.stimChan))
                etaWithinCING_low = [etaWithinCING_low; etaWithin];
                etaAcrossCING_low = [etaAcrossCING_low; etaAcross];
                tbsCING_low = [tbsCING_low; evokedAmp];
                CING_label_low = [CING_label_low; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                CING_unilabel_low = [CING_unilabel_low; subjAnatName];
                CING_uniMod_low = [CING_uniMod_low; tbsStatus];
                CING_coord_low = [CING_coord_low; subjElecCoord];
                CING_elecLabel_low = [CING_elecLabel_low subjElecLabel];
                CING_modStatusWithin_low = [CING_modStatusWithin_low; modStatusWithin];
                CING_modStatusAcross_low = [CING_modStatusAcross_low; modStatusAcross];
                CING_modStatusBoth_low = [CING_modStatusBoth_low; modStatusBoth];

                cnt_cing(1) = cnt_cing(1) + 1;
            elseif any(match_str(CENTRAL,subjGroup{subj}{i}.demo.stimChan))
                etaWithinCENTRAL_low = [etaWithinCENTRAL_low; etaWithin];
                etaAcrossCENTRAL_low = [etaAcrossCENTRAL_low; etaAcross];
                tbsCENTRAL_low = [tbsCENTRAL_low; evokedAmp];
                CENTRAL_label_low = [CENTRAL_label_low; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                CENTRAL_unilabel_low = [CENTRAL_unilabel_low; subjAnatName];
                CENTRAL_uniMod_low = [CENTRAL_uniMod_low; tbsStatus];
                CENTRAL_coord_low = [CENTRAL_coord_low; subjElecCoord];
                CENTRAL_elecLabel_low = [CENTRAL_elecLabel_low subjElecLabel];
                CENTRAL_modStatusWithin_low = [CENTRAL_modStatusWithin_low; modStatusWithin];
                CENTRAL_modStatusAcross_low = [CENTRAL_modStatusAcross_low; modStatusAcross];
                CENTRAL_modStatusBoth_low = [CENTRAL_modStatusBoth_low; modStatusBoth];

                cnt_central(1) = cnt_central(1) + 1;
            end
        else
            etaWithin_high = [etaWithin_high; etaWithin_TBS];
            etaAcross_high = [etaAcross_high; etaAcross_TBS];
            highCurrentWithin = [highCurrentWithin; withinNum / totalChanNum];
            highCurrentAcross = [highCurrentAcross; acrossNum / totalChanNum];
            highCurrentBoth =[highCurrentBoth; bothNum / totalChanNum];
            highCurrentChanDyna = [highCurrentChanDyna; numChanModDyna];
            highCurrentPercChanDyan = [highCurrentPercChanDyan; perChanModDyna];
            highCurrentChanMod = [highCurrentChanMod; numChanMod];
            highCurrentPerChanMod = [highCurrentPerChanMod; perChanMod];
            highCurrentStimIdx = [highCurrentStimIdx; {[num2str(subj) '-' 'sesh' num2str(sessionNum) '-' chanName{stimIdx}]}];
            highCurrentStimLabel = [highCurrentStimLabel; subjGroup{subj}{i}.demo.parcLabel(stimIdx)];
            highCurrentEtaWithinGroup = [highCurrentEtaWithinGroup; etaWithin];
            highCurrentEtaAcrossGroup = [highCurrentEtaAcrossGroup; etaAcross];
            highCurrentNetwork = [highCurrentNetwork; subjGroupRest{subj}.restPre.cw(stimIdx,2)];
            subjLabelMatched = [subjLabelMatched;subj];
            groupEvokedAmpHigh = [groupEvokedAmpHigh; nanmean(evokedAmp)];

            % RAW NUMBERS
            numNonMod_high = [numNonMod_high; nonModNum];
            numTBS_high = [numTBS_high; numChanMod];
            numWithin_high = [numWithin_high; withinNum];
            numAcross_high = [numAcross_high; acrossNum];
            numNonTBS_high = [numNonTBS_high; length(find(noTBSidx))];
            numBoth_high = [numBoth_high; bothNum];

            % GROUP BY STIM
            if any(match_str(DLPFC,subjGroup{subj}{i}.demo.stimChan)) 
                etaWithinDLPFC_high = [etaWithinDLPFC_high; etaWithin];
                etaAcrossDLPFC_high = [etaAcrossDLPFC_high; etaAcross];
                tbsDLPFC_high = [tbsDLPFC_high; evokedAmp];
                DLPFC_label_high = [DLPFC_label_high; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                DLPFC_unilabel_high = [DLPFC_unilabel_high; subjAnatName];
                DLPFC_uniMod_high = [DLPFC_uniMod_high; tbsStatus];
                DLPFC_coord_high = [DLPFC_coord_high; subjElecCoord];
                DLPFC_elecLabel_high = [DLPFC_elecLabel_high subjElecLabel];
                DLPFC_modStatusWithin_high = [DLPFC_modStatusWithin_high; modStatusWithin];
                DLPFC_modStatusAcross_high = [DLPFC_modStatusAcross_high; modStatusAcross];
                DLPFC_modStatusBoth_high = [DLPFC_modStatusBoth_high; modStatusBoth];

                cnt_dlpfc(2) = cnt_dlpfc(2) + 1;
            elseif any(match_str(VLPFC,subjGroup{subj}{i}.demo.stimChan))
                etaWithinVLPFC_high = [etaWithinVLPFC_high; etaWithin];
                etaAcrossVLPFC_high = [etaAcrossVLPFC_high; etaAcross];
                tbsVLPFC_high = [tbsVLPFC_high; evokedAmp];
                VLPFC_label_high = [VLPFC_label_high; subjAnatName];
                

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                VLPFC_unilabel_high = [VLPFC_unilabel_high; subjAnatName];
                VLPFC_uniMod_high = [VLPFC_uniMod_high; tbsStatus];
                VLPFC_coord_high = [VLPFC_coord_high; subjElecCoord];
                VLPFC_elecLabel_high = [VLPFC_elecLabel_high subjElecLabel];
                VLPFC_modStatusWithin_high = [VLPFC_modStatusWithin_high; modStatusWithin];
                VLPFC_modStatusAcross_high = [VLPFC_modStatusAcross_high; modStatusAcross];
                VLPFC_modStatusBoth_high = [VLPFC_modStatusBoth_high; modStatusBoth];

                cnt_opercular(2) = cnt_opercular(2) + 1;
            elseif any(match_str(OFC,subjGroup{subj}{i}.demo.stimChan))
                etaWithinOFC_high = [etaWithinOFC_high; etaWithin];
                etaAcrossOFC_high = [etaAcrossOFC_high; etaAcross];
                tbsOFC_high = [tbsOFC_high; evokedAmp];
                OFC_label_high = [OFC_label_high; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                OFC_unilabel_high = [OFC_unilabel_high; subjAnatName];
                OFC_uniMod_high = [OFC_uniMod_high; tbsStatus];
                OFC_coord_high = [OFC_coord_high; subjElecCoord];
                OFC_elecLabel_high = [OFC_elecLabel_high subjElecLabel];
                OFC_modStatusWithin_high = [OFC_modStatusWithin_high; modStatusWithin];
                OFC_modStatusAcross_high = [OFC_modStatusAcross_high; modStatusAcross];
                OFC_modStatusBoth_high = [OFC_modStatusBoth_high; modStatusBoth];

                cnt_ofc(2) = cnt_ofc(2) + 1;
            elseif any(match_str(TEMPORAL,subjGroup{subj}{i}.demo.stimChan))
                etaWithinTEMPORAL_high = [etaWithinTEMPORAL_high; etaWithin];
                etaAcrossTEMPORAL_high = [etaAcrossTEMPORAL_high; etaAcross];
                tbsTEMPORAL_high = [tbsTEMPORAL_high; evokedAmp];
                TEMPORAL_label_high = [TEMPORAL_label_high; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                TEMPORAL_unilabel_high = [TEMPORAL_unilabel_high; subjAnatName];
                TEMPORAL_uniMod_high = [TEMPORAL_uniMod_high; tbsStatus];
                TEMPORAL_coord_high = [TEMPORAL_coord_high; subjElecCoord];
                TEMPORAL_elecLabel_high = [TEMPORAL_elecLabel_high subjElecLabel];
                TEMPORAL_modStatusWithin_high = [TEMPORAL_modStatusWithin_high; modStatusWithin];
                TEMPORAL_modStatusAcross_high = [TEMPORAL_modStatusAcross_high; modStatusAcross];
                TEMPORAL_modStatusBoth_high = [TEMPORAL_modStatusBoth_high; modStatusBoth];

                cnt_temporal(2) = cnt_temporal(2) + 1;
            elseif any(match_str(INS,subjGroup{subj}{i}.demo.stimChan))
                etaWithinINS_high = [etaWithinINS_high; etaWithin];
                etaAcrossINS_high = [etaAcrossINS_high; etaAcross];
                tbsINS_high= [tbsINS_high; evokedAmp];
                INS_label_high = [INS_label_high; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                INS_unilabel_high = [INS_unilabel_high; subjAnatName];
                INS_uniMod_high = [INS_uniMod_high; tbsStatus];
                INS_coord_high = [INS_coord_high; subjElecCoord];
                INS_elecLabel_high = [INS_elecLabel_high subjElecLabel];
                INS_modStatusWithin_high = [INS_modStatusWithin_high; modStatusWithin];
                INS_modStatusAcross_high = [INS_modStatusAcross_high; modStatusAcross];
                INS_modStatusBoth_high = [INS_modStatusBoth_high; modStatusBoth];

                cnt_ins(2) = cnt_ins(2) + 1;
            elseif any(match_str(CING,subjGroup{subj}{i}.demo.stimChan))
                etaWithinCING_high = [etaWithinCING_high; etaWithin];
                etaAcrossCING_high = [etaAcrossCING_high; etaAcross];
                tbsCING_high = [tbsCING_high; evokedAmp];
                CING_label_high = [CING_label_high; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                CING_unilabel_high = [CING_unilabel_high; subjAnatName];
                CING_uniMod_high = [CING_uniMod_high; tbsStatus];
                CING_coord_high = [CING_coord_high; subjElecCoord];
                CING_elecLabel_high = [CING_elecLabel_high subjElecLabel];
                CING_modStatusWithin_high = [CING_modStatusWithin_high; modStatusWithin];
                CING_modStatusAcross_high = [CING_modStatusAcross_high; modStatusAcross];
                CING_modStatusBoth_high = [CING_modStatusBoth_high; modStatusBoth];

                cnt_cing(2) = cnt_cing(2) + 1;
            elseif any(match_str(CENTRAL,subjGroup{subj}{i}.demo.stimChan))
                etaWithinCENTRAL_high = [etaWithinCENTRAL_high; etaWithin];
                etaAcrossCENTRAL_high = [etaAcrossCENTRAL_high; etaAcross];
                tbsCENTRAL_high = [tbsCENTRAL_high; evokedAmp];
                CENTRAL_label_high = [CENTRAL_label_high; subjAnatName];

                % calc electrode modulated
                regionModulation = zeros(length(unique(subjAnatName)), 1);
                subjUniLabel = unique(subjAnatName);
                for a = 1:length(subjUniLabel)
                    labelIdx = match_str(subjAnatName, subjUniLabel{a});
                    regionModulation(a) = nanmean(tbsStatus(labelIdx));
                end

                CENTRAL_unilabel_high = [CENTRAL_unilabel_high; subjAnatName];
                CENTRAL_uniMod_high = [CENTRAL_uniMod_high; tbsStatus];
                CENTRAL_coord_high = [CENTRAL_coord_high; subjElecCoord];
                CENTRAL_elecLabel_high = [CENTRAL_elecLabel_high subjElecLabel];
                CENTRAL_modStatusWithin_high = [CENTRAL_modStatusWithin_high; modStatusWithin];
                CENTRAL_modStatusAcross_high = [CENTRAL_modStatusAcross_high; modStatusAcross];
                CENTRAL_modStatusBoth_high = [CENTRAL_modStatusBoth_high; modStatusBoth];
                cnt_central(2) = cnt_central(2) + 1;
            end
        end
    end
end

matchedStimIdx = match_str(erase(lowCurrentStimIdx,' '),erase(highCurrentStimIdx,' '));
matchedChanName  = lowCurrentStimIdx(matchedStimIdx);
matchedPerChanMod = zeros(length(matchedStimIdx),2);
matchedPerChanMod(:,1) = lowCurrentPerChanMod(matchedStimIdx);
matchedPerChanMod(:,2) = highCurrentPerChanMod;
matchedPerChanModDyna = zeros(length(matchedStimIdx),2);
matchedPerChanModDyna(:,1) = lowCurrentPercChanDyan(matchedStimIdx);
matchedPerChanModDyna(:,2) = highCurrentPercChanDyan;
matchedNumChanDyna = zeros(length(matchedStimIdx),2);
matchedNumChanDyna(:,1) = lowCurrentChanDyna(matchedStimIdx);
matchedNumChanDyna(:,2) = highCurrentChanDyna;
matchedNumChanMod = zeros(length(matchedStimIdx),2);
matchedNumChanMod(:,1) = lowCurrentChanMod(matchedStimIdx);
matchedNumChanMod(:,2) = highCurrentChanMod;
matchedNumWithin = zeros(length(matchedStimIdx),2);
matchedNumWithin(:,1) = lowCurrentWithin(matchedStimIdx);
matchedNumWithin(:,2) = highCurrentWithin;
matchedNumAcross = zeros(length(matchedStimIdx),2);
matchedNumAcross(:,1) = lowCurrentAcross(matchedStimIdx);
matchedNumAcross(:,2) = highCurrentAcross;
matchedNumBoth = zeros(length(matchedStimIdx),2);
matchedNumBoth(:,1) = lowCurrentBoth(matchedStimIdx);
matchedNumBoth(:,2) = highCurrentBoth;
matchedEvokedAmp = zeros(length(matchedStimIdx),2);
matchedEvokedAmp(:,1) = groupEvokedAmpLow(matchedStimIdx);
matchedEvokedAmp(:,2) = groupEvokedAmpHigh;

% COLOR
clr = [];
clr{1} = {'Color', [52, 152, 219]/255}; 
clr{2} = {'Color', [39, 174, 96]/255}; 
clr{3} = {'Color', [231, 76, 60]/255}; 
clr{4} = {'Color', [108, 92, 231]/255};
clr{5} = {'Color', [253, 203, 110]/255};
clr{6} = {'Color', [127, 140, 141]/255};
clr{7} = {'Color', [211, 84, 0]/255};
clr{8} = {'Color', [46, 204, 113]/255};
clr{9} = {'Color', [39, 174, 96]/255};
clr{10} = {'Color', [75, 75, 75]/255};

%% XBURST DYNAMICS

xBurstPattern = normalize(groupBurstHeatmap,3);
xBurstPattern = squeeze(nanmean(xBurstPattern,2));
kCluster = 6;

calinskiEvaluation = evalclusters(xBurstPattern,"kmeans","CalinskiHarabasz", ...
    "KList",1:kCluster);
daviesEvaluation = evalclusters(xBurstPattern,"kmeans","DaviesBouldin", ...
    "KList",1:kCluster);
gapEvaluation = evalclusters(xBurstPattern,"kmeans","gap","KList",1:kCluster);
silhouetteEvaluation = evalclusters(xBurstPattern,"kmeans","silhouette", ...
    "KList",1:kCluster);

figure;
title("Optimal Number of Clusters for Different Criteria")
colors = lines(4);

% Calinski-Harabasz Criterion Plot
subplot(2,2,1)
h1 = plot(calinskiEvaluation);
h1.Color = colors(1,:);
hold on
xline(calinskiEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

% Davies-Bouldin Criterion Plot
subplot(2,2,2)
h2 = plot(daviesEvaluation);
h2.Color = colors(2,:);
hold on
xline(daviesEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

% Gap Criterion Plot
subplot(2,2,3)
h3 = plot(gapEvaluation);
h3.Color = colors(3,:);
hold on
xline(gapEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

% Silhouette Criterion Plot
subplot(2,2,4)
h4 = plot(silhouetteEvaluation);
h4.Color = colors(4,:);
hold on
xline(silhouetteEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

ckSTIM_saveFig(['burstClusterEval'],20,20,300,'',1,[10,10],[]);
close all;

% PLOT DIFFERENT CLUSTERED PATTERNS
nClust = 3;
idx = kmeans(xBurstPattern,nClust);

figure;
for n = 1:nClust
    errorbar(nanmean(xBurstPattern(idx==n,:)),std(xBurstPattern(idx==n,:),[],1)/sqrt(length(find(idx==n)))); hold on;
end
box off; xlim([0.5 5.5]); xticks([1 2 3 4 5]); xlabel('Burst Number'); ylabel('Post-Burst Response (Z)');
legend({'Cluster 1', 'Cluster 2', 'Cluster 3'});

ckSTIM_saveFig(['burstClusters'],20,20,300,'',1,[10,10],[]);
close all;

% PIE DISTRIBUTRION
posTrend = length(find(idx==1));
negTrend = length(find(idx==2));
inbetweenTrend = length(find(idx==3));
totalChan = length(idx);

figure;
X = [posTrend negTrend inbetweenTrend];
pieObj = pie(X, [0 0 0], {[num2str(posTrend) ' (' num2str(round(posTrend/totalChan,2)*100) '%)' ' clust1'],[num2str(negTrend),...
    ' (' num2str(round(negTrend/totalChan,2)*100) '%)' ' clust2'],[num2str(inbetweenTrend) ' (' num2str(round(inbetweenTrend/totalChan,2)*100) '%)' ' clust3']});
title('Direction of Burst modulation');

ckSTIM_saveFig(['burstDynamic-pie'],10,10,0,'',1,[8, 8],[]);
close all;

%% XTRAIN DYNAMICS

% GET OVERALL DIRECTION OF CHANGE
xTrainPattern = normalize(groupTrainHeatmap,2);
heatmapReshape = zeros(size(groupTrainHeatmap,1), size(groupTrainHeatmap,2)*size(groupTrainHeatmap,3));
cnt = 1;
for t = 1:5:size(heatmapReshape,2)
    heatmapReshape(:,t:t+size(groupTrainHeatmap,3)-1) = xTrainPattern(:,cnt,:);
    cnt = cnt + 1;
end

intervalOne = 1:10; intervalThree = 40:50;
interval_Z = zeros(length(heatmapReshape),1);
interval_p = zeros(length(heatmapReshape),1);

for c = 1:length(heatmapReshape)
    [p, ~, stats] = ranksum(heatmapReshape(c,intervalOne), heatmapReshape(c,intervalThree));
    interval_Z(c) = stats.zval;
    interval_p(c) = p;
end

% PLOT DIFFERENT CLUSTERED PATTERNS
nClust = 3;
xTrainPatternMean = nanmean(xTrainPattern,3);

figure;
for n = 1:nClust
    if n == 1
        idx = interval_p <= 0.05 & interval_Z>0;
    elseif n == 2
        idx = interval_p <= 0.05 & interval_Z<0;
    else
        idx = interval_p > 0.05;
    end
    errorbar(nanmean(xTrainPatternMean(idx,:)),std(xTrainPatternMean(idx,:),[],1)/sqrt(length(find(idx)))); hold on;
end
box off; xlim([0.5 10.5]); xticks([1:10]); xlabel('Train Number'); ylabel('Post-Burst Response (Z)');
legend({'Neg trend', 'Pos trend', 'Inbetween'});

ckSTIM_saveFig(['trainClustersManual'],20,20,300,'',1,[10,10],[]);
close all;

% PIE DISTRIBUTRION
posTrend = length(find(interval_p <= 0.05 & interval_Z<0));
negTrend = length(find(interval_p <= 0.05 & interval_Z>0));
inbetweenTrend = length(find(interval_p > 0.05));
totalChan = length(interval_p);

figure;
X = [posTrend negTrend inbetweenTrend];
pieObj = pie(X, [0 0 0], {[num2str(posTrend) ' (' num2str(round(posTrend/totalChan,2)*100) '%)' ' pos trend'],[num2str(negTrend),...
    ' (' num2str(round(negTrend/totalChan,2)*100) '%)' ' neg trend'],[num2str(inbetweenTrend) ' (' num2str(round(inbetweenTrend/totalChan,2)*100) '%)' ' Inbetween']});
title('Direction of Train modulation');

ckSTIM_saveFig(['trainDynamic-pie'],10,10,0,'',1,[8, 8],[]);
close all;


% CLUSTERING EVAL

xTrainPattern = normalize(groupTrainHeatmap,2);
xTrainPattern = squeeze(nanmean(xTrainPattern,3));
kCluster = 20;

calinskiEvaluation = evalclusters(xTrainPattern,"kmeans","CalinskiHarabasz", ...
    "KList",1:kCluster);
daviesEvaluation = evalclusters(xTrainPattern,"kmeans","DaviesBouldin", ...
    "KList",1:kCluster);
gapEvaluation = evalclusters(xTrainPattern,"kmeans","gap","KList",1:kCluster);
silhouetteEvaluation = evalclusters(xTrainPattern,"kmeans","silhouette", ...
    "KList",1:kCluster);

figure;
title(t,"Optimal Number of Clusters for Different Criteria")
colors = lines(4);

% Calinski-Harabasz Criterion Plot
subplot(2,2,1)
h1 = plot(calinskiEvaluation);
h1.Color = colors(1,:);
hold on
xline(calinskiEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

% Davies-Bouldin Criterion Plot
subplot(2,2,2)
h2 = plot(daviesEvaluation);
h2.Color = colors(2,:);
hold on
xline(daviesEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

% Gap Criterion Plot
subplot(2,2,3)
h3 = plot(gapEvaluation);
h3.Color = colors(3,:);
hold on
xline(gapEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

% Silhouette Criterion Plot
subplot(2,2,4)
h4 = plot(silhouetteEvaluation);
h4.Color = colors(4,:);
hold on
xline(silhouetteEvaluation.OptimalK,"--","Optimal K", ...
    "LabelVerticalAlignment","middle")
hold off

ckSTIM_saveFig(['trainClustersEval'],20,20,300,'',1,[10,10],[]);
close all;




%% DECODING

% TBS+ vs TBS - CHANNELS
tStatGroup = [];
rocX = cell(1,3); rocY = cell(1,3); rocAUC = cell(1,3);
for d = 1:2
    if d == 1
        distLimitIdx = decodeDistGroup < 30;
        pred = decodeFeaturesGroup(distLimitIdx,:);
        resp = logical(decodeOutcomesGroup(distLimitIdx,1));
    elseif d == 2
        distLimitIdx = decodeDistGroup > 30;
        pred = decodeFeaturesGroup(distLimitIdx,:);
        resp = logical(decodeOutcomesGroup(distLimitIdx,1));
    else
        distLimitIdx = decodeDistGroup < 30;
        %pred = 1./exp(decodeDistGroup/100); 
        pred = decodeFeaturesGroup(distLimitIdx,1);
        resp = logical(decodeOutcomesGroup(distLimitIdx,1));
    end
    
    % CV
    num_folds = 10;
    indices = crossvalind('Kfold', resp, num_folds);
    scores = [];
    newResp = [];
    CVFold = 10;
    outcomeCV = cell(CVFold,1);
    scoreCV = cell(CVFold,1);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
        
        [post_score] = predict(mdl, pred(test,:));
        
        scores = cat(1, scores, post_score);
        newResp = cat(1, newResp, resp(test));
        outcomeCV{i} = logical(resp(test));
        scoreCV{i} = post_score;

        tStatGroup = cat(2, tStatGroup, mdl.Coefficients.tStat(2:7));
    end
    [rocX{d},rocY{d},~,rocAUC{d}] = perfcurve(outcomeCV, scoreCV,'true');
end

figure;
for d = 1:2
    plot(rocX{d}(:,1),rocY{d}(:,1), 'LineWidth', 2,'Color', clr{d}{2}); hold on;
    %ck_shadedErrorBar(rocX{d}(:,1), rocY{d}(:,1), rocY{d}(:,3)-rocY{d}(:,2),clr{d},1); hold on;
end
plot(0:1, 0:1, 'k:');
box off;
for d = 1:2
    if d == 1
        rocText{d} = ['Dist <30mm' num2str(round(rocAUC{d}(1),2)) ' (' num2str(round(rocAUC{d}(2),2)) '-' num2str(round(rocAUC{d}(3),2)) ')' ' AUC'];
    elseif d == 2
        rocText{d} = ['Dist >30mm' num2str(round(rocAUC{d}(1),2)) ' (' num2str(round(rocAUC{d}(2),2)) '-' num2str(round(rocAUC{d}(3),2)) ')' ' AUC'];
    else
        rocText{d} = ['Dist only' num2str(round(rocAUC{d}(1),2)) ' (' num2str(round(rocAUC{d}(2),2)) '-' num2str(round(rocAUC{d}(3),2)) ')' ' AUC'];
    end
    text(0.5, 0.25-d*0.05, rocText{d}, 'Units', 'Normalized' , 'Color', clr{d}{2}); 
end
title('Prediction of TBS+channels');
box off; xlabel('False positive rate'); ylabel('True positive rate');

ckSTIM_saveFig(['TBSROC'],20,20,300,'',1,[10,10],[]);
close all;

% T-STATS
featuresLabel = {'distance', 'CCEP duration', 'CCEP amplitude', 'PLV', 'PEA', 'WM'};

figure;
h = barwitherr(std(tStatGroup,[],2)/sqrt(size(tStatGroup,1)), 1:size(tStatGroup,1), nanmean(tStatGroup,2));
box off;
title('Model T-statistic');
box off; xlabel('Features'); ylabel('Feature Tstat');
xticklabels(featuresLabel); ylim([-4 16]);
ckSTIM_saveFig(['TBS_Feature_tStat'],20,20,300,'',1,[10,10],[]);
close all;

% LEAVE ONE OUT DECODING
aucDelta = zeros(3,size(decodeFeaturesGroup,2)+1);
for d = 1:size(decodeFeaturesGroup,2)+1
    featSelect = 1:size(decodeFeaturesGroup,2);
    featSelect(featSelect==d) = [];
    pred = decodeFeaturesGroup(:,featSelect);
    resp = logical(decodeOutcomesGroup(:,1));
    % CV
    num_folds = 10;
    indices = crossvalind('Kfold', resp, num_folds);
    scores = [];
    newResp = [];
    CVFold = 10;
    outcomeCV = cell(CVFold,1);
    scoreCV = cell(CVFold,1);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
        
        [post_score] = predict(mdl, pred(test,:));
        
        scores = cat(1, scores, post_score);
        newResp = cat(1, newResp, resp(test));
        outcomeCV{i} = logical(resp(test));
        scoreCV{i} = post_score;
    end
    [~,~,~,aucDelta(:,d)] = perfcurve(outcomeCV, scoreCV,'true');
end
aucDelta = aucDelta(:,end) - aucDelta;
aucErr = aucDelta(3,:) - aucDelta(2,:);
aucChange = aucDelta(1,:);

figure;
h = barwitherr(aucErr(1:end-1), 1:size(aucErr,2)-1, aucChange(1:end-1));
box off;
title('Features LOOCV');
box off; xlabel('Features'); ylabel('Change in AUC');
xticklabels(featuresLabel)
ckSTIM_saveFig(['TBS_Feature_LOOCV'],20,20,300,'',1,[10,10],[]);
close all;


% SINGLE DECODING
for d = 1:size(decodeFeaturesGroup,2)
    rocX = cell(1,2); rocY = cell(1,2); rocAUC = cell(1,2);
    for n = 1:2
        if n == 1
            distLimitIdx = decodeDistGroup < 30;
        else
            distLimitIdx = decodeDistGroup >= 30;
        end
        pred = decodeFeaturesGroup(distLimitIdx,d);
        resp = logical(decodeOutcomesGroup(distLimitIdx,1));
    
        % CV
        num_folds = 10;
        indices = crossvalind('Kfold', resp, num_folds);
        scores = [];
        newResp = [];
        CVFold = 10;
        outcomeCV = cell(CVFold,1);
        scoreCV = cell(CVFold,1);
        for i = 1:num_folds
            test = (indices == i); train = ~test;
            mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
            
            [post_score] = predict(mdl, pred(test,:));
            
            scores = cat(1, scores, post_score);
            newResp = cat(1, newResp, resp(test));
            outcomeCV{i} = logical(resp(test));
            scoreCV{i} = post_score;
        end
        [rocX{n},rocY{n},~,rocAUC{n}] = perfcurve(outcomeCV, scoreCV,'true');
    end

    figure;
    for n = 1:2
        plot(rocX{n}(:,1),rocY{n}(:,1), 'LineWidth', 2,'Color', clr{n}{2}); hold on;
    end
    plot(0:1, 0:1, 'k:');
    box off;
    for n = 1:2
        if n == 1
            rocText{n} = ['Dist <30mm' num2str(round(rocAUC{n}(1),2)) ' (' num2str(round(rocAUC{n}(2),2)) '-' num2str(round(rocAUC{n}(3),2)) ')' ' AUC'];
        else
            rocText{n} = ['Dist >30mm' num2str(round(rocAUC{n}(1),2)) ' (' num2str(round(rocAUC{n}(2),2)) '-' num2str(round(rocAUC{n}(3),2)) ')' ' AUC'];
        end
        text(0.5, 0.25-n*0.05, rocText{n}, 'Units', 'Normalized' , 'Color', clr{n}{2}); 
    end
    title(featuresLabel{d});
    box off; xlabel('False positive rate'); ylabel('True positive rate');
    
    ckSTIM_saveFig(['TBSROC-' featuresLabel{d}],20,20,300,'',1,[10,10],[]);
    close all;
end



%% MOD+ VS MOD- CHANNELS DECODING
tStatGroup = [];
rocX = cell(1,2); rocY = cell(1,2); rocAUC = cell(1,2);
for d = 1:2
    if d == 1
        distLimitIdx = decodeDistGroup < 30;
    else
        distLimitIdx = decodeDistGroup > 30;
    end

    pred = decodeFeaturesGroup(logical(decodeOutcomesGroup(:,1)) & distLimitIdx,:);
    resp = logical(decodeOutcomesGroup(logical(decodeOutcomesGroup(:,1)) & distLimitIdx,5));
    % CV
    num_folds = 10;
    indices = crossvalind('Kfold', resp, num_folds);
    scores = [];
    newResp = [];
    CVFold = 10;
    outcomeCV = cell(CVFold,1);
    scoreCV = cell(CVFold,1);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
        
        [post_score] = predict(mdl, pred(test,:));
        
        scores = cat(1, scores, post_score);
        newResp = cat(1, newResp, resp(test));
        outcomeCV{i} = logical(resp(test));
        scoreCV{i} = post_score;

        tStatGroup = cat(2, tStatGroup, mdl.Coefficients.tStat(2:7));
    end
    [rocX{d},rocY{d},~,rocAUC{d}] = perfcurve(outcomeCV, scoreCV,'true');
end

figure;
for d = 1:2
    plot(rocX{d}(:,1),rocY{d}(:,1), 'LineWidth', 2,'Color', clr{d}{2}); hold on;
    %ck_shadedErrorBar(rocX{d}(:,1), rocY{d}(:,1), rocY{d}(:,3)-rocY{d}(:,2),clr{d},1); hold on;
end
plot(0:1, 0:1, 'k:');
box off;
rocText{1} = ['Dist <30mm' num2str(round(rocAUC{1}(1),2)) ' (' num2str(round(rocAUC{1}(2),2)) '-' num2str(round(rocAUC{1}(3),2)) ')' ' AUC'];
rocText{2} = ['Dist >30mm' num2str(round(rocAUC{2}(1),2)) ' (' num2str(round(rocAUC{2}(2),2)) '-' num2str(round(rocAUC{2}(3),2)) ')' ' AUC'];
text(0.5, 0.25, rocText{1}, 'Units', 'Normalized' , 'Color', clr{1}{2}); title('Prediction of MOD+channels');
text(0.5, 0.15, rocText{2}, 'Units', 'Normalized', 'Color', clr{2}{2}); 
box off; xlabel('False positive rate'); ylabel('True positive rate');

ckSTIM_saveFig(['MODROC'],20,20,300,'',1,[10,10],[]);
close all;

% T-STATS
featuresLabel = {'distance', 'CCEP duration', 'CCEP amplitude', 'PLV', 'PEA', 'WM'};

figure;
h = barwitherr(std(tStatGroup,[],2)/sqrt(size(tStatGroup,1)), 1:size(tStatGroup,1), nanmean(tStatGroup,2));
box off;
title('MODvsMOD- T-statistic');
box off; xlabel('Features'); ylabel('Feature Tstat');
xticklabels(featuresLabel);
ckSTIM_saveFig(['MOD_Feature_tStat'],20,20,300,'',1,[10,10],[]);
close all;

% LEAVE ONE OUT DECODING
featuresLabel = {'distance', 'CCEP duration', 'CCEP amplitude', 'PLV', 'PEA', 'WM'};
aucDelta = zeros(3,size(decodeFeaturesGroup,2)+1);
for d = 1:size(decodeFeaturesGroup,2)+1
    featSelect = 1:size(decodeFeaturesGroup,2);
    featSelect(featSelect==d) = [];
    pred = decodeFeaturesGroup(logical(decodeOutcomesGroup(:,1)),featSelect);
    resp = logical(decodeOutcomesGroup(logical(decodeOutcomesGroup(:,1)),5));

    % CV
    num_folds = 10;
    indices = crossvalind('Kfold', resp, num_folds);
    scores = [];
    newResp = [];
    CVFold = 10;
    outcomeCV = cell(CVFold,1);
    scoreCV = cell(CVFold,1);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
        
        [post_score] = predict(mdl, pred(test,:));
        
        scores = cat(1, scores, post_score);
        newResp = cat(1, newResp, resp(test));
        outcomeCV{i} = logical(resp(test));
        scoreCV{i} = post_score;
    end
    [~,~,~,aucDelta(:,d)] = perfcurve(outcomeCV, scoreCV,'true');
end
aucDelta = aucDelta(:,end) - aucDelta;
aucErr = aucDelta(3,:) - aucDelta(2,:);
aucChange = aucDelta(1,:);

figure;
h = barwitherr(aucErr(1:end-1), 1:size(aucErr,2)-1, aucChange(1:end-1));
box off;
title('Features LOOCV');
box off; xlabel('Features'); ylabel('Change in AUC');
xticklabels(featuresLabel);
ckSTIM_saveFig(['MOD_Feature_LOOCV'],20,20,300,'',1,[10,10],[]);
close all;


% SINGLE DECODING
for d = 1:size(decodeFeaturesGroup,2)
    rocX = cell(1,2); rocY = cell(1,2); rocAUC = cell(1,2);
    for n = 1:2
        if n == 1
            distLimitIdx = decodeDistGroup < 30;
        else
            distLimitIdx = decodeDistGroup >= 30;
        end
        pred = decodeFeaturesGroup(logical(decodeOutcomesGroup(:,1)) & distLimitIdx,d);
        resp = logical(decodeOutcomesGroup(logical(decodeOutcomesGroup(:,1)) & distLimitIdx,5));
    
        % CV
        num_folds = 10;
        indices = crossvalind('Kfold', resp, num_folds);
        scores = [];
        newResp = [];
        CVFold = 10;
        outcomeCV = cell(CVFold,1);
        scoreCV = cell(CVFold,1);
        for i = 1:num_folds
            test = (indices == i); train = ~test;
            mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
            
            [post_score] = predict(mdl, pred(test,:));
            
            scores = cat(1, scores, post_score);
            newResp = cat(1, newResp, resp(test));
            outcomeCV{i} = logical(resp(test));
            scoreCV{i} = post_score;
        end
        [rocX{n},rocY{n},~,rocAUC{n}] = perfcurve(outcomeCV, scoreCV,'true');
    end

    figure;
    for n = 1:2
        plot(rocX{n}(:,1),rocY{n}(:,1), 'LineWidth', 2,'Color', clr{n}{2}); hold on;
    end
    plot(0:1, 0:1, 'k:');
    box off;
    for n = 1:2
        if n == 1
            rocText{n} = ['Dist <30mm' num2str(round(rocAUC{n}(1),2)) ' (' num2str(round(rocAUC{n}(2),2)) '-' num2str(round(rocAUC{n}(3),2)) ')' ' AUC'];
        else
            rocText{n} = ['Dist >30mm' num2str(round(rocAUC{n}(1),2)) ' (' num2str(round(rocAUC{n}(2),2)) '-' num2str(round(rocAUC{n}(3),2)) ')' ' AUC'];
        end
        text(0.5, 0.25-n*0.05, rocText{n}, 'Units', 'Normalized' , 'Color', clr{n}{2}); 
    end
    title(featuresLabel{d});
    box off; xlabel('False positive rate'); ylabel('True positive rate');
    
    ckSTIM_saveFig(['MODROC-' featuresLabel{d}],20,20,300,'',1,[10,10],[]);
    close all;
end

%% WITHIN VS ACROSS CHANNELS
tStatGroup = [];
rocX = cell(1,2); rocY = cell(1,2); rocAUC = cell(1,2);
for d = 1:2
    if d == 1
        distLimitIdx = decodeDistGroup < 30;
    else
        distLimitIdx = decodeDistGroup > 30;
    end

    pred = decodeFeaturesGroup(logical(decodeOutcomesGroup(:,5)) & distLimitIdx,:);
    resp = ~logical(decodeOutcomesGroup(logical(decodeOutcomesGroup(:,5)) & distLimitIdx,3));
    % CV
    num_folds = 10;
    indices = crossvalind('Kfold', resp, num_folds);
    scores = [];
    newResp = [];
    CVFold = 10;
    outcomeCV = cell(CVFold,1);
    scoreCV = cell(CVFold,1);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
        
        [post_score] = predict(mdl, pred(test,:));
        
        scores = cat(1, scores, post_score);
        newResp = cat(1, newResp, resp(test));
        outcomeCV{i} = logical(resp(test));
        scoreCV{i} = post_score;

        tStatGroup = cat(2, tStatGroup, mdl.Coefficients.tStat(2:7));
    end
    [rocX{d},rocY{d},~,rocAUC{d}] = perfcurve(outcomeCV, scoreCV,'true');
end

figure;
for d = 1:2
    plot(rocX{d}(:,1),rocY{d}(:,1), 'LineWidth', 2,'Color', clr{d}{2}); hold on;
end
plot(0:1, 0:1, 'k:');
box off;
rocText{1} = ['Dist <30mm' num2str(round(rocAUC{1}(1),2)) ' (' num2str(round(rocAUC{1}(2),2)) '-' num2str(round(rocAUC{1}(3),2)) ')' ' AUC'];
rocText{2} = ['Dist >30mm' num2str(round(rocAUC{2}(1),2)) ' (' num2str(round(rocAUC{2}(2),2)) '-' num2str(round(rocAUC{2}(3),2)) ')' ' AUC'];
text(0.5, 0.25, rocText{1}, 'Units', 'Normalized' , 'Color', clr{1}{2}); title('Prediction of WITHIN VS ACROSS channels');
text(0.5, 0.15, rocText{2}, 'Units', 'Normalized', 'Color', clr{2}{2}); 
box off; xlabel('False positive rate'); ylabel('True positive rate');

ckSTIM_saveFig(['WITHINvACROSSROC'],20,20,300,'',1,[10,10],[]);
close all;

% T-STATS
featuresLabel = {'distance', 'CCEP duration', 'CCEP amplitude', 'PLV', 'PEA', 'WM'};

figure;
h = barwitherr(std(tStatGroup,[],2)/sqrt(size(tStatGroup,1)), 1:size(tStatGroup,1), nanmean(tStatGroup,2));
box off;
title('Types of Mod T-statistic');
box off; xlabel('Features'); ylabel('Feature Tstat');
xticklabels(featuresLabel);
ckSTIM_saveFig(['TYPESMOD_Feature_tStat'],20,20,300,'',1,[10,10],[]);
close all;


% LEAVE ONE OUT DECODING
featuresLabel = {'distance', 'CCEP duration', 'CCEP amplitude', 'PLV', 'PEA', 'WM'};
aucDelta = zeros(3,size(decodeFeaturesGroup,2)+1);
for d = 1:size(decodeFeaturesGroup,2)+1
    featSelect = 1:size(decodeFeaturesGroup,2);
    featSelect(featSelect==d) = [];
    pred = decodeFeaturesGroup(logical(decodeOutcomesGroup(:,5)) & distLimitIdx,featSelect);
    resp = logical(decodeOutcomesGroup(logical(decodeOutcomesGroup(:,5)) & distLimitIdx,3));

    % CV
    num_folds = 10;
    indices = crossvalind('Kfold', resp, num_folds);
    scores = [];
    newResp = [];
    CVFold = 10;
    outcomeCV = cell(CVFold,1);
    scoreCV = cell(CVFold,1);
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
        
        [post_score] = predict(mdl, pred(test,:));
        
        scores = cat(1, scores, post_score);
        newResp = cat(1, newResp, resp(test));
        outcomeCV{i} = logical(resp(test));
        scoreCV{i} = post_score;
    end
    [~,~,~,aucDelta(:,d)] = perfcurve(outcomeCV, scoreCV,'true');
end
aucDelta = aucDelta(:,end) - aucDelta;
aucErr = aucDelta(3,:) - aucDelta(2,:);
aucChange = aucDelta(1,:);

figure;
h = barwitherr(aucErr(1:end-1), 1:size(aucErr,2)-1, aucChange(1:end-1));
box off;
title('Features LOOCV');
box off; xlabel('Features'); ylabel('Change in AUC');
xticklabels(featuresLabel);
ckSTIM_saveFig(['MODType_Feature_LOOCV'],20,20,300,'',1,[10,10],[]);
close all;

% SINGLE DECODING
for d = 1:size(decodeFeaturesGroup,2)
    rocX = cell(1,2); rocY = cell(1,2); rocAUC = cell(1,2);
    for n = 1:2
        if n == 1
            distLimitIdx = decodeDistGroup < 30;
        else
            distLimitIdx = decodeDistGroup >= 30;
        end
        pred = decodeFeaturesGroup(logical(decodeOutcomesGroup(:,5)) & distLimitIdx,d);
        resp = logical(decodeOutcomesGroup(logical(decodeOutcomesGroup(:,5)) & distLimitIdx,3));

        % CV
        num_folds = 10;
        indices = crossvalind('Kfold', resp, num_folds);
        scores = [];
        newResp = [];
        CVFold = 10;
        outcomeCV = cell(CVFold,1);
        scoreCV = cell(CVFold,1);
        for i = 1:num_folds
            test = (indices == i); train = ~test;
            mdl = fitglm(pred(train,:), resp(train), 'Distribution', 'binomial', 'Link', 'logit');
            
            [post_score] = predict(mdl, pred(test,:));
            
            scores = cat(1, scores, post_score);
            newResp = cat(1, newResp, resp(test));
            outcomeCV{i} = logical(resp(test));
            scoreCV{i} = post_score;
        end
        [rocX{n},rocY{n},~,rocAUC{n}] = perfcurve(outcomeCV, scoreCV,'true');
    end

    figure;
    for n = 1:2
        plot(rocX{n}(:,1),rocY{n}(:,1), 'LineWidth', 2,'Color', clr{n}{2}); hold on;
    end
    plot(0:1, 0:1, 'k:');
    box off;
    for n = 1:2
        if n == 1
            rocText{n} = ['Dist <30mm' num2str(round(rocAUC{n}(1),2)) ' (' num2str(round(rocAUC{n}(2),2)) '-' num2str(round(rocAUC{n}(3),2)) ')' ' AUC'];
        else
            rocText{n} = ['Dist >30mm' num2str(round(rocAUC{n}(1),2)) ' (' num2str(round(rocAUC{n}(2),2)) '-' num2str(round(rocAUC{n}(3),2)) ')' ' AUC'];
        end
        text(0.5, 0.25-n*0.05, rocText{n}, 'Units', 'Normalized' , 'Color', clr{n}{2}); 
    end
    title(featuresLabel{d});
    box off; xlabel('False positive rate'); ylabel('True positive rate');
    
    ckSTIM_saveFig(['MODTYPESROC-' featuresLabel{d}],20,20,300,'',1,[10,10],[]);
    close all;
end



%% PIE CHART SUMMARY
% 1mA+2mA combined
totalChanLowHigh = sum([numTBS_low; numTBS_high]);
pctNonMod = num2str(sum([numNonMod_low; numNonMod_high]) / totalChanLowHigh * 100);
pctWithin = num2str(sum([numWithin_low; numWithin_high]) / totalChanLowHigh * 100);
pctAcross = num2str(sum([numAcross_low; numAcross_high]) / totalChanLowHigh * 100);
pctBoth = num2str(sum([numBoth_low; numBoth_high]) / totalChanLowHigh * 100);
figure;
X = [sum([numWithin_low; numWithin_high]) sum([numAcross_low; numAcross_high]) sum([numBoth_low; numBoth_high]) sum([numNonMod_low; numNonMod_high])];
pieObj = pie(X, [1 1 1 0], {[num2str(sum([numWithin_low; numWithin_high])) ' (' pctWithin '%)' ' Within Plasticity'],[num2str(sum([numAcross_low; numAcross_high])),...
    ' (' pctAcross '%)' ' Across Plasticity'],[num2str(sum([numBoth_low; numBoth_high])) ' (' pctBoth '%)' ' Both'],[num2str(sum([numNonMod_low; numNonMod_high])),...
    ' (' pctNonMod '%)' ' No Plasticity']});
pieObj(1).FaceColor = [26 188 156]./255; pieObj(7).FaceColor = [243 241 239]./255; title('Plasticity %');

ckSTIM_saveFig(['plasticity-breakdown-single'],10,10,0,'',1,[8, 8],[]);
close all;

figure;
% 1mA
subplot(1,2,1)
totalTBSChan = sum(numTBS_low);
pctNonMod = num2str(sum(numNonMod_low) / totalTBSChan * 100);
pctWithin = num2str(sum(numWithin_low) / totalTBSChan * 100);
pctAcross = num2str(sum(numAcross_low) / totalTBSChan * 100);
pctBoth = num2str(sum(numBoth_low) / totalTBSChan * 100);
X = [sum(numWithin_low) sum(numAcross_low) sum(numBoth_low) sum(numNonMod_low)];
pieObj = pie(X, [1 1 1 0], {[num2str(sum(numWithin_low)) ' (' pctWithin '%)' ' Within Plasticity'],[num2str(sum(numAcross_low)),...
    ' (' pctAcross '%)' ' Across Plasticity'],[num2str(sum(numBoth_low)) ' (' pctBoth '%)' ' Both'],[num2str(sum(numNonMod_low)),...
    ' (' pctNonMod '%)' ' No Plasticity']});
pieObj(1).FaceColor = [26 188 156]./255; pieObj(7).FaceColor = [243 241 239]./255; title('1mA');

% 2mA
subplot(1,2,2)
totalTBSChan = sum(numTBS_high);
pctNonMod = num2str(sum(numNonMod_high) / totalTBSChan * 100);
pctWithin = num2str(sum(numWithin_high) / totalTBSChan * 100);
pctAcross = num2str(sum(numAcross_high) / totalTBSChan * 100);
pctBoth = num2str(sum(numBoth_high) / totalTBSChan * 100);
X = [sum(numWithin_high) sum(numAcross_high) sum(numBoth_high) sum(numNonMod_high)];
pieObj = pie(X, [1 1 1 0], {[num2str(sum(numWithin_high)) ' (' pctWithin '%)' ' Within Plasticity'],[num2str(sum(numAcross_high)),...
    ' (' pctAcross '%)' ' Across Plasticity'],[num2str(sum(numBoth_high)) ' (' pctBoth '%)' ' Both'],[num2str(sum(numNonMod_high)),...
    ' (' pctNonMod '%)' ' No Plasticity']});
pieObj(1).FaceColor = [26 188 156]./255; pieObj(7).FaceColor = [243 241 239]./255; title('2mA');

ckSTIM_saveFig(['plasticity-breakdown'],10,10,0,'',1,[10, 8],[]);
close all;

% LOOK AT TBS vs NON-TBS
% 1ma + 2ma combined
totalTBSChan = sum([numTBS_low; numTBS_high]);
totalNonTBSChan = sum([numNonTBS_low; numNonTBS_high]);
totalChan = totalTBSChan + totalNonTBSChan;
pctTBS = num2str(totalTBSChan / totalChan * 100);
pctNonTBS = num2str(totalNonTBSChan / totalChan * 100);
figure;
X = [totalTBSChan totalNonTBSChan];
pieObj = pie(X, [1 0], {[num2str(totalTBSChan) ' (' pctTBS '%)' ' TBS+'],[num2str(totalNonTBSChan),...
    ' (' pctNonTBS '%)' ' No TBS response']});
pieObj(1).FaceColor = [26 188 156]./255; title('1mA+2mA combined');

ckSTIM_saveFig(['TBS-breakdown-single'],10,10,0,'',1,[8, 8],[]);
close all;

figure;
% 1mA
subplot(1,2,1)
totalTBSChan = sum(numTBS_low);
totalNonTBSChan = sum(numNonTBS_low);
totalChan = totalTBSChan + totalNonTBSChan;
pctTBS = num2str(totalTBSChan / totalChan * 100);
pctNonTBS = num2str(totalNonTBSChan / totalChan * 100);
X = [totalTBSChan totalNonTBSChan];
pieObj = pie(X, [1 0], {[num2str(totalTBSChan) ' (' pctTBS '%)' ' TBS+'],[num2str(totalNonTBSChan),...
    ' (' pctNonTBS '%)' ' No TBS response']});
pieObj(1).FaceColor = [26 188 156]./255; title('1mA');

% 2mA
subplot(1,2,2)
totalTBSChan = sum(numTBS_high);
totalNonTBSChan = sum(numNonTBS_high);
totalChan = totalTBSChan + totalNonTBSChan;
pctTBS = num2str(totalTBSChan / totalChan * 100);
pctNonTBS = num2str(totalNonTBSChan / totalChan * 100);
X = [totalTBSChan totalNonTBSChan];
pieObj = pie(X, [1 0], {[num2str(totalTBSChan) ' (' pctTBS '%)' ' TBS+'],[num2str(totalNonTBSChan),...
    ' (' pctNonTBS '%)' ' No TBS response']});
pieObj(1).FaceColor = [26 188 156]./255; title('2mA');

ckSTIM_saveFig(['TBS-breakdown'],10,10,0,'',1,[10, 8],[]);
close all;


%% STIM WM VS TBS
figure;
plot(groupStimWM, groupEvokedAmp, 'o'); hold on;
for s = 1:length(subjGroup)
    plot(groupStimWM(subjLabel==s),groupEvokedAmp(subjLabel==s),'o', clr{s}{:});
end
[R P] = corrcoef(groupStimWM, groupEvokedAmp, 'rows','complete');
c = polyfit(groupStimWM, groupEvokedAmp, 1); y_est = polyval(c,groupStimWM);
plot(groupStimWM, y_est, 'k-', 'LineWidth', 1); box off;
text(0.2,0.8,['R:' num2str(R(2))],'Units','normalized');
text(0.2,0.7,['pval:' num2str(P(2))],'Units','normalized');
xlabel('WM percentage'); ylabel('TBS Amplitude')

ckSTIM_saveFig(fullfile(['StimWM-TBS-Scatter']),10,10,0,'',1,[6,5],[]);
close all;


%% SCATTERPLOT OF ETA

figure;
%1mA
subplot(1,2,1)
within = etaWithin_low;
across = etaAcross_low;
notNanIdx = ~isnan(within);
[R P] = corrcoef(within, across, 'rows','complete');
c = polyfit(within(notNanIdx), across(notNanIdx), 1); y_est = polyval(c,within(notNanIdx));
plot(within, across, 'o', 'MarkerSize', 5, 'Color', [189, 195, 199]/255);  hold on;
plot(within(notNanIdx), y_est, 'k-', 'LineWidth', 1); box off;
text(0.2,0.8,['R:' num2str(R(2))],'Units','normalized');
text(0.2,0.7,['pval:' num2str(P(2))],'Units','normalized');
title('1mA - only tbs channels'); xlabel('WithinBurst eta'); ylabel('AcrossBurst eta');

%2mA
subplot(1,2,2)
within = etaWithin_high;
across = etaAcross_high;
notNanIdx = ~isnan(within);
[R P] = corrcoef(within, across, 'rows','complete');
c = polyfit(within(notNanIdx), across(notNanIdx), 1); y_est = polyval(c,within(notNanIdx));
plot(within, across, 'o', 'MarkerSize', 5, 'Color', [189, 195, 199]/255);  hold on;
plot(within(notNanIdx), y_est, 'k-', 'LineWidth', 1); box off;
text(0.2,0.8,['R:' num2str(R(2))],'Units','normalized');
text(0.2,0.7,['pval:' num2str(P(2))],'Units','normalized');
title('2mA - only tbs channels'); xlabel('WithinBurst eta'); ylabel('AcrossBurst eta');

ckSTIM_saveFig(fullfile(['WithinAcross-Scatter']),10,10,0,'',1,[12,5],[]);
close all;


%% TBS RESPONSE PLOT
figure;
% EVOKED RESPONSE BY DISTANCE
subplot(2,3,1);
X1=[1.2];
Y1=groupDIST(:,6);
X2=[1.3];
Y2=groupDIST(:,5);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('Distance to Stim Site');
for c = 1:length(groupDIST)
    plot([X1,X2],[groupDIST(c,6) groupDIST(c,5)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.025*X1 X2+0.02*X1]);
[pval,~, stat] = signrank(groupDIST(:,6),groupDIST(:,5));
text(0.5,0.8,['Z' '(' num2str(length(groupDIST)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% EVOKED RESPONSE BY RESTING DATA
subplot(2,3,2);
X1=[-5];
Y1=groupPLV(:,6);
X2=[8];
Y2=groupPLV(:,5);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('PLV');
for c = 1:length(groupPLV)
    plot([X1,X2],[groupPLV(c,6) groupPLV(c,5)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+1.2*X1 X2-0.6*X1]);
[pval,~, stat] = signrank(groupPLV(:,6),groupPLV(:,5));
text(0.5,0.8,['Z' '(' num2str(length(groupPLV)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');


% EVOKED RESPONSE BY RESTING DATA CORRELATION
subplot(2,3,3);
X1=[-5];
Y1=groupPEA(:,6);
X2=[35];
Y2=groupPEA(:,5);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('Voltage Correlation');
for c = 1:length(groupPEA)
    plot([X1,X2],[groupPEA(c,6) groupPEA(c,5)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+3*X1 X2-1.5*X1]);
[pval,~, stat] = signrank(groupPEA(:,6),groupPEA(:,5));
text(0.5,0.8,['Z' '(' num2str(length(groupPEA)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% EVOKED RESPONSE BY AMPLITUDE
subplot(2,3,4);
X1=[-2];
Y1=groupAMP(:,6);
X2=[1];
Y2=groupAMP(:,5);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('CCEP Amplitude');
for c = 1:length(groupAMP)
    plot([X1,X2],[groupAMP(c,6) groupAMP(c,5)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+0.7*X1 X2-0.2*X1]);
[pval,~, stat] = signrank(groupAMP(:,6),groupAMP(:,5));
text(0.5,0.8,['Z' '(' num2str(length(groupAMP)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

subplot(2,3,5);
X1=[0];
Y1=groupLAT(:,6);
X2=[10];
Y2=groupLAT(:,5);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('CCEP Duration');
for c = 1:length(groupLAT)
    plot([X1,X2],[groupLAT(c,6) groupLAT(c,5)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.3*X2 X2+0.2*X2]);
[pval,~, stat] = signrank(groupLAT(:,6),groupLAT(:,5));
text(0.5,0.8,['Z' '(' num2str(length(groupLAT)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

subplot(2,3,6);
X1=[0];
Y1=groupWM(:,6);
X2=[1.5];
Y2=groupWM(:,5);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('WM Ratio');
for c = 1:length(groupWM)
    plot([X1,X2],[groupWM(c,6) groupWM(c,5)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.35*X2 X2+0.2*X2]);
[pval,~, stat] = signrank(groupWM(:,6),groupWM(:,5));
text(0.5,0.8,['Z' '(' num2str(length(groupWM)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');


ckSTIM_saveFig(['features-by-TBSresponse'],10,10,0,'',1,[10, 8],[]);
close all;


%% BASELINE FEATURES BY GROUPED PLASTICITY
figure;
% EVOKED RESPONSE BY DISTANCE
subplot(2,3,1);
X1=[1.2];
Y1=groupDIST(~isnan(groupDIST(:,2)),2);
X2=[1.25];
Y2=groupDIST(~isnan(groupDIST(:,1)),1);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'MOD-', 'MOD+'}); ylabel('Distance to Stim Site');
for c = 1:length(groupDIST)
    plot([X1,X2],[groupDIST(c,2) groupDIST(c,1)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.013*X1 X2+0.013*X1]);
[pval,~, stat] = signrank(groupDIST(:,2),groupDIST(:,1));
text(0.5,0.8,['Z' '(' num2str(length(groupDIST)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% EVOKED RESPONSE BY RESTING DATA
subplot(2,3,2);
X1=[-3];
Y1=groupPLV(~isnan(groupPLV(:,2)),2);
X2=[5];
Y2=groupPLV(~isnan(groupPLV(:,1)),1);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'MOD-', 'MOD+'}); ylabel('PLV');
for c = 1:length(groupPLV)
    plot([X1,X2],[groupPLV(c,2) groupPLV(c,1)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+1*X1 X2-0.6*X1]);
[pval,~, stat] = signrank(groupPLV(:,2),groupPLV(:,1));
text(0.5,0.8,['Z' '(' num2str(length(groupPLV)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% EVOKED RESPONSE BY RESTING DATA CORRELATION
subplot(2,3,3);
X1=[-5];
Y1=groupPEA(~isnan(groupPEA(:,2)),2);
X2=[8];
Y2=groupPEA(~isnan(groupPEA(:,1)),1);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'MOD-', 'MOD+'}); ylabel('Voltage correlation');
for c = 1:length(groupPEA)
    plot([X1,X2],[groupPEA(c,2) groupPEA(c,1)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+1*X1 X2-0.8*X1]);
[pval,~, stat] = signrank(groupPEA(:,2),groupPEA(:,1));
text(0.5,0.8,['Z' '(' num2str(length(groupPEA)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% EVOKED RESPONSE BY AMPLITUDE
subplot(2,3,4);
X1=[0];
Y1=groupAMP(~isnan(groupAMP(:,2)),2);
X2=[0.6];
Y2=groupAMP(~isnan(groupAMP(:,1)),1);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'MOD-', 'MOD+'}); ylabel('CCEP Amplitude');
for c = 1:length(groupAMP)
    plot([X1,X2],[groupAMP(c,2) groupAMP(c,1)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.4*X2 X2+0.2*X2]);
[pval,~, stat] = signrank(groupAMP(:,2),groupAMP(:,1));
text(0.5,0.8,['Z' '(' num2str(length(groupAMP)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

subplot(2,3,5);
X1=[0];
Y1=groupLAT(~isnan(groupLAT(:,2)),2);
X2=[5];
Y2=groupLAT(~isnan(groupLAT(:,1)),1);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'MOD-', 'MOD+'}); ylabel('CCEP Duration');
for c = 1:length(groupLAT)
    plot([X1,X2],[groupLAT(c,2) groupLAT(c,1)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.3*X2 X2+0.2*X2]);
[pval,~, stat] = signrank(groupLAT(:,2),groupLAT(:,1));
text(0.5,0.8,['Z' '(' num2str(length(groupLAT)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

subplot(2,3,6);
X1=[0];
Y1=groupWM(~isnan(groupWM(:,2)),2);
X2=[1.5];
Y2=groupWM(~isnan(groupWM(:,1)),1);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'TBS-', 'TBS+'}); ylabel('WM Ratio');
for c = 1:length(groupWM)
    plot([X1,X2],[groupWM(c,2) groupWM(c,1)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.4*X2 X2+0.3*X2]);
[pval,~, stat] = signrank(groupWM(:,2),groupWM(:,1));
text(0.5,0.8,['Z' '(' num2str(length(groupWM)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

ckSTIM_saveFig(['features-by-MODresponse'],10,10,0,'',1,[10, 8],[]);
close all;

%% BASELINE FEATURES BY WITHIN VS ACROSS PLASTICITY
figure;
% EVOKED RESPONSE BY DISTANCE
subplot(2,3,1);
X1=[1.2];
Y1=groupDIST(~isnan(groupDIST(:,3)),3);
X2=[1.25];
Y2=groupDIST(~isnan(groupDIST(:,4)),4);
X3=[1.3];
Y3=groupDIST(~isnan(groupDIST(:,7)),7);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
Hdl3=violinChart(gca,X3,Y3,[1 1 1]);
xticks([X1, X2, X3]); xticklabels({'WITHIN', 'ACROSS', 'BOTH'}); ylabel('Distance to Stim Site');
for c = 1:length(groupDIST)
    plot([X1, X2, X3],[groupDIST(c,3) groupDIST(c,4) groupDIST(c,7)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.02*X1 X3+0.02*X1]);
[p,anovatab,stats] = kruskalwallis([groupDIST(:,3),groupDIST(:,4),groupDIST(:,7)],[],'off');
text(0.5,0.8,['chisqr:' num2str(anovatab{2,5})],'Units','normalized');
text(0.5,0.7,['pval:' num2str(anovatab{2,6})],'Units','normalized');
 
% EVOKED RESPONSE BY RESTING DATA
subplot(2,3,2);
X1=[-1];
Y1=groupPLV(~isnan(groupPLV(:,3)),3);
X2=[3];
Y2=groupPLV(~isnan(groupPLV(:,4)),4);
X3=[7];
Y3=groupPLV(~isnan(groupPLV(:,7)),7);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
Hdl3=violinChart(gca,X3,Y3,[1 1 1]);
xticks([X1,X2,X3]); xticklabels({'WITHIN', 'ACROSS', 'BOTH'}); ylabel('PLV');
for c = 1:length(groupPLV)
    plot([X1,X2,X3],[groupPLV(c,3) groupPLV(c,4) groupPLV(c,7)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+1.5*X1 X3-0.9*X1]);
[p,anovatab,stats] = kruskalwallis([groupPLV(:,3),groupPLV(:,4),groupPLV(:,7)],[],'off');
text(0.2,0.8,['chisqr:' num2str(anovatab{2,5})],'Units','normalized');
text(0.2,0.7,['pval:' num2str(anovatab{2,6})],'Units','normalized');

% EVOKED RESPONSE BY RESTING DATA CORRELATION
subplot(2,3,3);
X1=[-1];
Y1=groupPEA(~isnan(groupPEA(:,3)),3);
X2=[5];
Y2=groupPEA(~isnan(groupPEA(:,4)),4);
X3=[11];
Y3=groupPLV(~isnan(groupPEA(:,7)),7);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
Hdl3=violinChart(gca,X3,Y3,[1 1 1]);
xticks([X1, X2, X3]); xticklabels({'WITHIN', 'ACROSS', 'BOTH'}); ylabel('Coherence');
for c = 1:length(groupPEA)
    plot([X1,X2,X3],[groupPEA(c,3) groupPEA(c,4),groupPEA(c,7)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1+3*X1 X3-1.5*X1]);
[p,anovatab,stats] = kruskalwallis([groupPEA(:,3),groupPEA(:,4),groupPEA(:,7)],[],'off');
text(0.2,0.8,['chisqr:' num2str(anovatab{2,5})],'Units','normalized');
text(0.2,0.7,['pval:' num2str(anovatab{2,6})],'Units','normalized');

% EVOKED RESPONSE BY AMPLITUDE
subplot(2,3,4);
X1=[0];
Y1=groupAMP(~isnan(groupAMP(:,3)),3);
X2=[0.2];
Y2=groupAMP(~isnan(groupAMP(:,4)),4);
X3 = 0.4;
Y3 =groupAMP(~isnan(groupAMP(:,7)),7);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
Hdl3=violinChart(gca,X3,Y3,[1 1 1]);
xticks([X1, X2, X3]); xticklabels({'WITHIN', 'ACROSS', 'BOTH'}); ylabel('CCEP Amplitude');
for c = 1:length(groupAMP)
    plot([X1,X2,X3],[groupAMP(c,3) groupAMP(c,4) groupAMP(c,7)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.4*X2 X3+0.2*X2]);
[p,anovatab,stats] = kruskalwallis([groupAMP(:,3), groupAMP(:,4), groupAMP(:,7)],[],'off');
text(0.2,0.8,['chisqr:' num2str(anovatab{2,5})],'Units','normalized');
text(0.2,0.7,['pval:' num2str(anovatab{2,6})],'Units','normalized');


subplot(2,3,5);
X1=[0];
Y1=groupLAT(~isnan(groupLAT(:,3)),3);
X2=[4];
Y2=groupLAT(~isnan(groupLAT(:,4)),4);
X3=[8];
Y3=groupLAT(~isnan(groupLAT(:,7)),7);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
Hdl2=violinChart(gca,X3,Y3,[1 1 1]);
xticks([X1, X2,X3]); xticklabels({'WITHIN', 'ACROSS', 'BOTH'}); ylabel('CCEP Duration');
for c = 1:length(groupLAT)
    plot([X1,X2,X3],[groupLAT(c,3) groupLAT(c,4) groupLAT(c,7)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.3*X2 X3+0.2*X2]);
[p,anovatab,stats] = kruskalwallis([groupLAT(:,3),groupLAT(:,4),groupLAT(:,7)],[],'off');
text(0.2,0.8,['chisqr:' num2str(anovatab{2,5})],'Units','normalized');
text(0.2,0.7,['pval:' num2str(anovatab{2,6})],'Units','normalized');

subplot(2,3,6);
X1=[0];
Y1=groupWM(~isnan(groupWM(:,3)),3);
X2=[0.5];
Y2=groupWM(~isnan(groupWM(:,4)),4);
X3=[1]
Y3=groupWM(~isnan(groupWM(:,7)),7);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
Hdl3=violinChart(gca,X3,Y3,[1 1 1]);
xticks([X1, X2, X3]); xticklabels({'WITHIN', 'ACROSS', 'BOTH'}); ylabel('WM Ratio');
for c = 1:length(groupWM)
    plot([X1,X2,X3],[groupWM(c,3) groupWM(c,4) groupWM(c,7)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.5*X2 X3+0.2*X3]);
[p,anovatab,stats] = kruskalwallis([groupWM(:,3),groupWM(:,4),groupWM(:,7)],[],'off');
text(0.2,0.8,['chisqr:' num2str(anovatab{2,5})],'Units','normalized');
text(0.2,0.7,['pval:' num2str(anovatab{2,6})],'Units','normalized');

ckSTIM_saveFig(['features-by-WITHINACROSSresponse'],10,10,0,'',1,[10, 8],[]);
close all;


%% TBS VS CURRENT 
figure;
subplot(1,6,1);

X1=[1];
Y1=matchedPerChanMod(:,1);
X2=[10];
Y2=matchedPerChanMod(:,2);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'1mA', '2mA+'}); ylabel('% TBS+ channels');
for c = 1:length(matchedPerChanMod)
    plot([X1,X2],[matchedPerChanMod(c,1) matchedPerChanMod(c,2)], clr{subjLabelMatched(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.4*X2 X2+0.3*X2]);
[pval,~, stat] = signrank(matchedPerChanMod(:,1),matchedPerChanMod(:,2));
text(0.5,0.8,['Z' '(' num2str(length(matchedPerChanMod)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

subplot(1,6,2);

X1=[1];
Y1=matchedEvokedAmp(:,1);
X2=[2.5];
Y2=matchedEvokedAmp(:,2);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'1mA', '2mA+'}); ylabel('TBS Evoked Amplitude');
for c = 1:length(matchedEvokedAmp)
    plot([X1,X2],[matchedEvokedAmp(c,1) matchedEvokedAmp(c,2)], clr{subjLabelMatched(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.25*X2 X2+0.25*X2]);
[pval,~, stat] = signrank(matchedEvokedAmp(:,1),matchedEvokedAmp(:,2));
text(0.5,0.8,['Z' '(' num2str(length(matchedEvokedAmp)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');


% both modulation
subplot(1,6,3);

X1=[1];
Y1=matchedPerChanModDyna(:,1);
X2=[25];
Y2=matchedPerChanModDyna(:,2);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'1mA', '2mA+'}); ylabel('% modulated channels out of total channels');
for c = 1:length(matchedPerChanModDyna)
    plot([X1,X2],[matchedPerChanModDyna(c,1) matchedPerChanModDyna(c,2)], clr{subjLabelMatched(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.45*X2 X2+0.4*X2]);
[pval,~, stat] = signrank(matchedPerChanModDyna(:,1),matchedPerChanModDyna(:,2));
text(0.5,0.8,['Z' '(' num2str(length(matchedPerChanModDyna)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% WITHIN PLASTICITY
subplot(1,6,4);

X1=[1];
Y1=matchedNumWithin(:,1);
X2=[25];
Y2=matchedNumWithin(:,2);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'1mA', '2mA+'}); ylabel('% Across Burst Modulation');
for c = 1:length(matchedNumWithin)
    plot([X1,X2],[matchedNumWithin(c,1) matchedNumWithin(c,2)], clr{subjLabelMatched(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.45*X2 X2+0.3*X2]);
[pval,~, stat] = signrank(matchedNumWithin(:,1),matchedNumWithin(:,2)); if ~isfield(stat, 'zval') stat.zval = 0; end
text(0.5,0.8,['Z' '(' num2str(length(matchedNumWithin)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');


% ACROSS PLASTICITY
subplot(1,6,5);

X1=[1];
Y1=matchedNumAcross(:,1);
X2=[40];
Y2=matchedNumAcross(:,2);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'1mA', '2mA+'}); ylabel('% Across Train Modulation');
for c = 1:length(matchedNumAcross)
    plot([X1,X2],[matchedNumAcross(c,1) matchedNumAcross(c,2)], clr{subjLabelMatched(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.5*X2 X2+0.3*X2]);
[pval,~, stat] = signrank(matchedNumAcross(:,1),matchedNumAcross(:,2)); if ~isfield(stat, 'zval') stat.zval = 0; end
text(0.5,0.8,['Z' '(' num2str(length(matchedNumAcross)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');

% BOTH
subplot(1,6,6);

X1=[1];
Y1=matchedNumBoth(:,1);
X2=[40];
Y2=matchedNumBoth(:,2);
Hdl1=violinChart(gca,X1,Y1,[1 1 1]);
Hdl2=violinChart(gca,X2,Y2,[1 1 1]);
xticks([X1, X2]); xticklabels({'1mA', '2mA+'}); ylabel('% Both Modulation');
for c = 1:length(matchedNumBoth)
    plot([X1,X2],[matchedNumBoth(c,1) matchedNumBoth(c,2)], clr{subjLabelMatched(c)}{:},'Marker','o','MarkerSize',3, 'LineWidth', 0.1); 
end
xlim([X1-0.4*X2 X2+0.3*X2]);
[pval,~, stat] = signrank(matchedNumBoth(:,1),matchedNumBoth(:,2)); if ~isfield(stat, 'zval') stat.zval = 0; end
text(0.5,0.8,['Z' '(' num2str(length(matchedNumBoth)) '): ' num2str(stat.zval)],'Units','normalized');
text(0.5,0.7,['pval:' num2str(pval)],'Units','normalized');


ckSTIM_saveFig(['modulation-1v2ma'],10,10,300,'',1,[15, 5],[]);
close all;

%% PLASTICITY CHANNEL FEATURES
figure;
subplot(2,3,1);
% DISTANCE STRATIFIED BY PLASTIC VS NONPLASTIC CHANNELS
ck_make_bar_graph_quad(groupDIST(:,2),groupDIST(:,1),groupDIST(:,3),groupDIST(:,4),1,1,{'w';'w';'w';'w'},0)
xticklabels({'non-plasticity', '+plasticity', '+withinBurst', '+acrossBurst'}); ylabel('Distance to Stim Site');
for c = 1:length(groupDIST)
    plot(1:4,[groupDIST(c,2) groupDIST(c,1) groupDIST(c,3) groupDIST(c,4)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineStyle', 'none'); 
end

subplot(2,3,2);
% EVOKD AMP STRATIFIED BY PLASTIC VS NONPLASTIC CHANNELS
ck_make_bar_graph_quad(groupAMP(:,2),groupAMP(:,1),groupAMP(:,3),groupAMP(:,4),1,1,{'w';'w';'w';'w'},0)
xticklabels({'non-plasticity', '+plasticity', '+withinBurst', '+acrossBurst'}); ylabel('TBS Amplitude');
for c = 1:length(groupDIST)
    plot(1:4,[groupAMP(c,2) groupAMP(c,1) groupAMP(c,3) groupAMP(c,4)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineStyle', 'none'); 
end

% EVOKD AMP STRATIFIED BY PLASTIC VS NONPLASTIC CHANNELS
subplot(2,3,3);
ck_make_bar_graph_quad(groupLAT(:,2),groupLAT(:,1),groupLAT(:,3),groupLAT(:,4),1,1,{'w';'w';'w';'w'},0)
xticklabels({'non-plasticity', '+plasticity', '+withinBurst', '+acrossBurst'}); ylabel('CCEP Duration');
for c = 1:length(groupDIST)
    plot(1:4,[groupLAT(c,2) groupLAT(c,1) groupLAT(c,3) groupLAT(c,4)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineStyle', 'none'); 
end

subplot(2,3,4);
% RESTING PLV STRATIFIED BY PLASTIC VS NONPLASTIC CHANNELS
ck_make_bar_graph_quad(groupPLV(:,2),groupPLV(:,1),groupPLV(:,3),groupPLV(:,4),1,1,{'w';'w';'w';'w'},0)
xticklabels({'non-plasticity', '+plasticity', '+withinBurst', '+acrossBurst'}); ylabel('Resting PLV');
for c = 1:length(groupPLV)
    plot(1:4,[groupPLV(c,2) groupPLV(c,1) groupPLV(c,3) groupPLV(c,4)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineStyle', 'none'); 
end

subplot(2,3,5);
% RESTING CORRELATION STRATIFIED BY PLASTIC VS NONPLASTIC CHANNELS
ck_make_bar_graph_quad(groupPEA(:,2),groupPEA(:,1),groupPEA(:,3),groupPEA(:,4),1,1,{'w';'w';'w';'w'},0)
xticklabels({'non-plasticity', '+plasticity', '+withinBurst', '+acrossBurst'}); ylabel('Resting coh');
for c = 1:length(groupPEA)
    plot(1:4,[groupPEA(c,2) groupPEA(c,1) groupPEA(c,3) groupPEA(c,4)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineStyle', 'none'); 
end

subplot(2,3,6);
% WM
ck_make_bar_graph_quad(groupWM(:,2),groupWM(:,1),groupWM(:,3),groupWM(:,4),1,1,{'w';'w';'w';'w'},0)
xticklabels({'non-plasticity', '+plasticity', '+withinBurst', '+acrossBurst'}); ylabel('Elec WM');
for c = 1:length(groupWM)
    plot(1:4,[groupWM(c,2) groupWM(c,1) groupWM(c,3) groupWM(c,4)], clr{subjLabel(c)}{:},'Marker','o','MarkerSize',3, 'LineStyle', 'none'); 
end

ckSTIM_saveFig(['features-by-plasticity'],10,10,0,'',1,[10, 8],[]);
close all;

%% ANATOMICAL STIM EFFECTS ON BRAIN
DirVal='C:\Users\danny\Box\HalpernLabSecure\emuECoG\elec\fsaverage\surf\';
pialSuffix='';
if isfile(fullfile(DirVal,'rh.pial.T1'))
    pialSuffix = '.T1';
end
[verticesrh, facesrh] = freesurfer_read_surf( fullfile(DirVal, ['rh.pial' pialSuffix]) );
rpial=[]; rpial.vertices = verticesrh; rpial.faces = facesrh;
[verticeslh, faceslh] = freesurfer_read_surf( fullfile(DirVal, ['lh.pial' pialSuffix]) );
lpial=[]; lpial.vertices = verticeslh; lpial.faces = faceslh;

for anat = 1:7
    quantUse = []; rawLabels = []; quantElec = []; quantElecLabel = [];
    if anat == 1
        % DLPFC
        quantUse = [DLPFC_uniMod_low;DLPFC_uniMod_high];
        rawLabels = [DLPFC_label_low; DLPFC_label_high];
        quantElec = [DLPFC_coord_low; DLPFC_coord_high];
        quantElecLabel = [DLPFC_elecLabel_low DLPFC_elecLabel_high];
        note = 'DLPFC';
    elseif anat == 2
        % VLPFC
        quantUse = [VLPFC_uniMod_low;VLPFC_uniMod_high];
        rawLabels = [VLPFC_unilabel_low; VLPFC_unilabel_high];
        quantElec = [VLPFC_coord_low; VLPFC_coord_high];
        quantElecLabel = [VLPFC_elecLabel_low VLPFC_elecLabel_high];
        note = 'opercular';
    elseif anat == 3
        % TEMPORAL
        quantUse = [TEMPORAL_uniMod_low;TEMPORAL_uniMod_high];
        rawLabels = [TEMPORAL_unilabel_low; TEMPORAL_unilabel_high];
        note = 'temporal';
        quantElec = [TEMPORAL_coord_low; TEMPORAL_coord_high];
        quantElecLabel = [TEMPORAL_elecLabel_low TEMPORAL_elecLabel_high];
    elseif anat == 4
        % OFC
        quantUse = [OFC_uniMod_low;OFC_uniMod_high];
        rawLabels = [OFC_unilabel_low; OFC_unilabel_high];
        note = 'ofc';
        quantElec = [OFC_coord_low; OFC_coord_high];
        quantElecLabel = [OFC_elecLabel_low OFC_elecLabel_high];
    elseif anat == 5
        % INS
        quantUse = [INS_uniMod_low;INS_uniMod_high];
        rawLabels = [INS_unilabel_low; INS_unilabel_high];
        note = 'insula';
        quantElec = [INS_coord_low; INS_coord_high];
        quantElecLabel = [INS_elecLabel_low INS_elecLabel_high];
    elseif anat == 6
        % CING
        quantUse = [CING_uniMod_low;CING_uniMod_high];
        rawLabels = [CING_unilabel_low; CING_unilabel_high];
        note = 'cingulate';
        quantElec = [CING_coord_low; CING_coord_high];
        quantElecLabel = [CING_elecLabel_low CING_elecLabel_high];
    elseif anat == 7
        % CENTRAL
        quantUse = [CENTRAL_uniMod_low;CENTRAL_uniMod_high];
        rawLabels = [CENTRAL_unilabel_low; CENTRAL_unilabel_high];
        note = 'central';
        quantElec = [CENTRAL_coord_low; CENTRAL_coord_high];
        quantElecLabel = [CENTRAL_elecLabel_low CENTRAL_elecLabel_high];
    end

    % SWITCH RIGHT TO LEFT
    isRight = zeros(length(quantElecLabel),1);
    for chan = 1:length(quantElecLabel)
        if quantElecLabel{chan}(1) == 'R'
            isRight(chan) = 1;
        end
    end
    quantElec(logical(isRight),1) = quantElec(logical(isRight),1)*-1;
    nonanIdx = ~isnan(quantElec(:,1));
    quantUse(quantUse==0) = 0.01;
    quantUse(quantUse==1) = 10;

    figure;
    plot_mesh_brain(lpial,[-90 0],[],0.8); 
    plot_data_on_mesh(quantElec(nonanIdx,:), quantUse(nonanIdx),'surface','patchsize',30,'overlapmethod','mean');
    fix_lighting;
    caxis([0 10]);
    colorbar;
    ckSTIM_saveFig([note 'stimEngagementBrain-lateral'],15,15,300,'',4,[20,20],[]);
    close all;
    
    figure;
    plot_mesh_brain(lpial,[90 0],[],0.8); 
    plot_data_on_mesh(quantElec(nonanIdx,:), quantUse(nonanIdx),'surface','patchsize',20,'overlapmethod','mean');
    fix_lighting;
    caxis([0 10]);
    colorbar;
    ckSTIM_saveFig([note 'stimEngagementBrain-medial'],15,15,300,'',4,[20,20],[]);
    close all;
end

%% ANATOMICAL STIM

% GROUP THEM
labelReplace = cell(2,14);
labelReplace{1,1} = {'rostralmiddlefrontal', 'superiorfrontal','caudalmiddlefrontal'};
labelReplace{2,1} = 'DLPFC'; 
labelReplace{1,2} = {'parsopercularis', 'parsorbitalis','parstriangularis'};
labelReplace{2,2} = 'VLPFC';
labelReplace{1,3} = {'medialorbitofrontal'};
labelReplace{2,3} = 'VMPFC';
labelReplace{1,4} = {'lateralorbitofrontal'};
labelReplace{2,4} = 'OFC';
labelReplace{1,5} = {'caudalanteriorcingulate','isthmuscingulate','rostralanteriorcingulate'};
labelReplace{2,5} = 'ACC';
labelReplace{1,6} = {'posteriorcingulate'};
labelReplace{2,6} = 'PCC';
labelReplace{1,7} = {'amygdala', 'hippocamp', 'parahippocampal'};
labelReplace{2,7} = 'MesialTemporal';
labelReplace{1,8} = {'inferiortemporal','middletemporal','superiortemporal','transversetemporal','fusiform'};
labelReplace{2,8} = 'LateralTemporal';
labelReplace{1,9} =  {'supramarginal','superiorparietal','inferiorparietal', 'precuneus'};
labelReplace{2,9} = 'Parietal';
labelReplace{1,10} = {'insula'};
labelReplace{2,10} = 'Insula';
labelReplace{1,11} = {'lingual','cuneus'};
labelReplace{2,11} = 'Occipital';
labelReplace{1,12} = {'paracentral','precentral','postcentral'};
labelReplace{2,12} = 'Central';
labelReplace{1,13} = {'putamen','thalamus','pallidum','caudate','unk'};
labelReplace{2,13} = 'Subcortical';
labelReplace{1,14} = {'cerebellum'};
labelReplace{2,14} = 'Cerebellum';

orderLabel = [];
for l = 1:length(labelReplace)
    orderLabel = cat(1,orderLabel, labelReplace(2,l));
end

% LOOKING AT PROPORTION OF TBS+ CHANNELS & MOD STATUS
for anat = 1:7
    clear quantUse quantUse2 quantUse3 quantUse4 rawLabels
    if anat == 1
        % DLPFC
        quantUse = [DLPFC_uniMod_low;DLPFC_uniMod_high];
        quantUse2 = [DLPFC_modStatusWithin_low;DLPFC_modStatusWithin_high];
        quantUse3 = [DLPFC_modStatusAcross_low;DLPFC_modStatusAcross_high];
        quantUse4 = [DLPFC_modStatusBoth_low;DLPFC_modStatusBoth_high];
        rawLabels = [DLPFC_unilabel_low; DLPFC_unilabel_high];
        note = 'DLPFC';
        stimN = sum(cnt_dlpfc);
    elseif anat == 2
        % OPERCULAR
        quantUse = [VLPFC_uniMod_low;VLPFC_uniMod_high];
        quantUse2 = [VLPFC_modStatusWithin_low;VLPFC_modStatusWithin_high];
        quantUse3 = [VLPFC_modStatusAcross_low;VLPFC_modStatusAcross_high];
        quantUse4 = [VLPFC_modStatusBoth_low;VLPFC_modStatusBoth_high];
        rawLabels = [VLPFC_unilabel_low; VLPFC_unilabel_high];
        note = 'opercular';
        stimN = sum(cnt_opercular);
    elseif anat == 3
        % TEMPORAL
        quantUse = [TEMPORAL_uniMod_low;TEMPORAL_uniMod_high];
        quantUse2 = [TEMPORAL_modStatusWithin_low;TEMPORAL_modStatusWithin_high];
        quantUse3 = [TEMPORAL_modStatusAcross_low;TEMPORAL_modStatusAcross_high];
        quantUse4 = [TEMPORAL_modStatusBoth_low;TEMPORAL_modStatusBoth_high];
        rawLabels = [TEMPORAL_unilabel_low; TEMPORAL_unilabel_high];
        note = 'temporal';
        stimN = sum(cnt_temporal);
    elseif anat == 4
        % OFC
        quantUse = [OFC_uniMod_low;OFC_uniMod_high];
        quantUse2 = [OFC_modStatusWithin_low;OFC_modStatusWithin_high];
        quantUse3 = [OFC_modStatusAcross_low;OFC_modStatusAcross_high];
        quantUse4 = [OFC_modStatusBoth_low;OFC_modStatusBoth_high];
        rawLabels = [OFC_unilabel_low; OFC_unilabel_high];
        note = 'ofc';
        stimN = sum(cnt_ofc);
    elseif anat == 5
        % INS
        quantUse = [INS_uniMod_low;INS_uniMod_high];
        quantUse2 = [INS_modStatusWithin_low;INS_modStatusWithin_high];
        quantUse3 = [INS_modStatusAcross_low;INS_modStatusAcross_high];
        quantUse4 = [INS_modStatusBoth_low;INS_modStatusBoth_high];
        rawLabels = [INS_unilabel_low; INS_unilabel_high];
        note = 'insula';
        stimN = sum(cnt_ins);
    elseif anat == 6
        % CING
        quantUse = [CING_uniMod_low;CING_uniMod_high];
        quantUse2 = [CING_modStatusWithin_low;CING_modStatusWithin_high];
        quantUse3 = [CING_modStatusAcross_low;CING_modStatusAcross_high];
        quantUse4 = [CING_modStatusBoth_low;CING_modStatusBoth_high];
        rawLabels = [CING_unilabel_low; CING_unilabel_high];
        note = 'cingulate';
        stimN = sum(cnt_cing);
    elseif anat == 7
        % CENTRAL
        quantUse = [CENTRAL_uniMod_low;CENTRAL_uniMod_high];
        quantUse2 = [CENTRAL_modStatusWithin_low;CENTRAL_modStatusWithin_high];
        quantUse3 = [CENTRAL_modStatusAcross_low;CENTRAL_modStatusAcross_high];
        quantUse4 = [CENTRAL_modStatusBoth_low;CENTRAL_modStatusBoth_high];
        rawLabels = [CENTRAL_unilabel_low; CENTRAL_unilabel_high];
        note = 'central';
        stimN = sum(cnt_central);
    end

    % LABEL CONVERSION
    newLabels = rawLabels;
    for t = 1:length(labelReplace)
        replaceLabIdx = match_str(rawLabels,labelReplace{1,t});
        newLabels(replaceLabIdx) = {labelReplace{2,t}};
    end
   
    nonanIdx = ~isnan(quantUse);
    rawLabelsNoNaN = newLabels(nonanIdx);
    quantUseNoNaN = quantUse(nonanIdx);
    quantUseNoNaN2 = quantUse2(nonanIdx);
    quantUseNoNaN3 = quantUse3(nonanIdx);
    quantUseNoNaN4 = quantUse4(nonanIdx);
    uniqueLabels = unique(newLabels);
    quantList = []; quantList2 = []; quantList3 = []; quantList4 = [];
    quantListErr = [];
    newLabel = []; cnt = 1;
    for l = 1:length(uniqueLabels)
        labelidx = match_str(rawLabelsNoNaN,uniqueLabels{l});
        if length(labelidx) < 3
            continue;
        end
        quantList(cnt) = sum(quantUseNoNaN(labelidx)) / length(labelidx);
        quantListErr(cnt) = 0; %std(quantUseNoNaN(labelidx))/sqrt(length(labelidx));

        quantList2(cnt) = sum(quantUseNoNaN2(labelidx)) / length(labelidx);
        quantListErr2(cnt) = 0; %std(quantUseNoNaN(labelidx))/sqrt(length(labelidx));

        quantList3(cnt) = sum(quantUseNoNaN3(labelidx)) / length(labelidx);
        quantListErr3(cnt) = 0; %std(quantUseNoNaN(labelidx))/sqrt(length(labelidx));

        quantList4(cnt) = sum(quantUseNoNaN4(labelidx)) / length(labelidx);
        quantListErr4(cnt) = 0;

        newLabel{cnt} = uniqueLabels{l};
        cnt = cnt + 1;
    end
    ordIdx = [];
    for t = 1:length(orderLabel)
        index = find(strcmp(newLabel, orderLabel{t}));
        ordIdx = [ordIdx; index];
    end
    sortidx = ordIdx;
    sortedLabel = newLabel(sortidx);
    
    figure;
    stackedBar = [(quantList-quantList2-quantList3-quantList4);quantList3; quantList2;quantList4]';
    h = bar(stackedBar(sortidx,:), 'stacked');
    set(gca, 'XTick', 1:length(sortedLabel),'XTickLabel', sortedLabel); box off;
    xtickangle(45); title([note ' numStimSession:' num2str(stimN)]);
    ylim([0 1]);
    set(h, 'FaceColor', 'Flat');
    h(1).CData = [236, 240, 241]/255;
    h(2).CData = [39, 174, 96]/255;
    h(3).CData = [41, 128, 185]/255;
    h(4).CData = [231, 76, 60]/255;
   

ckSTIM_saveFig(['chanModByRegion-breakdown-' note],15,15,300,'',1,[10, 13],[]);
close all;
end

figure;
for anat = 1:7
    if anat == 1
        % DLPFC
        quantUse = [tbsDLPFC_low;tbsDLPFC_high];
        rawLabels = [DLPFC_label_low; DLPFC_label_high];
        note = 'DLPFC';
        stimN = sum(cnt_dlpfc);
    elseif anat == 2
        % OPERCULAR
        quantUse = [tbsVLPFC_low;tbsVLPFC_high];
        rawLabels = [VLPFC_label_low; VLPFC_label_high];
        note = 'opercular';
        stimN = sum(cnt_opercular);
    elseif anat == 3
        % TEMPORAL
        quantUse = [tbsTEMPORAL_low;tbsTEMPORAL_high];
        rawLabels = [TEMPORAL_label_low; TEMPORAL_label_high];
        note = 'temporal';
        stimN = sum(cnt_temporal);
    elseif anat == 4
        % OFC
        quantUse = [tbsOFC_low;tbsOFC_high];
        rawLabels = [OFC_label_low; OFC_label_high];
        note = 'ofc';
        stimN = sum(cnt_ofc);
    elseif anat == 5
        % INS
        quantUse = [tbsINS_low;tbsINS_high];
        rawLabels = [INS_label_low; INS_label_high];
        note = 'insula';
        stimN = sum(cnt_ins);
    elseif anat == 6
        % CING
        quantUse = [tbsCING_low;tbsCING_high];
        rawLabels = [CING_label_low; CING_label_high];
        note = 'cingulate';
        stimN = sum(cnt_cing);
    elseif anat == 7
        % CENTRAL
        quantUse = [tbsCENTRAL_low;tbsCENTRAL_high];
        rawLabels = [CENTRAL_label_low; CENTRAL_label_high];
        note = 'central';
        stimN = sum(cnt_central);
    end
   
    
    nonanIdx = ~isnan(quantUse);
    rawLabelsNoNaN = rawLabels(nonanIdx);
    quantUseNoNaN = quantUse(nonanIdx);
    uniqueLabels = unique(rawLabels);
    quantList = [];
    quantListErr = [];
    newLabel = []; cnt = 1;
    for l = 1:length(uniqueLabels)
        labelidx = match_str(rawLabelsNoNaN,uniqueLabels{l});
        if length(labelidx) <= 2
            continue;
        end
        quantList(cnt) = nanmean(quantUseNoNaN(labelidx));
        quantListErr(cnt) = std(quantUseNoNaN(labelidx))/sqrt(length(labelidx));
        newLabel{cnt} = uniqueLabels{l};
        cnt = cnt + 1;
    end
    useIdx = 1:5; % use top 5
    [~,sortidx] = sort(quantList,'descend');
    sortidx = sortidx(useIdx);
    sortedLabel = newLabel(sortidx);
    subplot(1,7,anat)
    h = barwitherr([quantListErr(sortidx)'], [1:length(sortedLabel)], [quantList(sortidx)']);
    set(gca, 'XTick', 1:length(sortedLabel),'XTickLabel', sortedLabel); box off;
    xtickangle(45); title([note ' stimulation' ' N:' num2str(stimN)]);
    set(h(1),'FaceColor','w'); ylim([0 15]);
end
ckSTIM_saveFig(['tbsamp-allstim'],10,10,0,'',1,[20, 5],[]);
close all;



