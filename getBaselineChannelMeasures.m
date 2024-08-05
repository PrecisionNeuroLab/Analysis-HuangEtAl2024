function getBaselineChannelMeasures(datName,dtype)

load(datName);

global ft_default
ft_default.showcallinfo = 'no';

%% GET ELEC DISTANCE
coordElec = ParcellationValues(:,5:7);
d = zeros(size(coordElec,1),size(coordElec,1));
for i=1:size(d,1)
    for j=1:size(d,2)
        d(i,j) = sqrt((coordElec(i,1)-coordElec(j,1)).^2+(coordElec(i,2)-coordElec(j,2)).^2+(coordElec(i,3)-coordElec(j,3)).^2);
    end
end


%% CALCULATE USING MANUAL METHOD
numChan = length(coordElec);

N   = 8;  % Order
Fc1 = 5;  % First Cutoff Frequency
Fc2 = 13;  % Second Cutoff Frequency
% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, FS);
Hd = design(h, 'butter');
[b,a] = sos2tf(Hd.sosMatrix, Hd.ScaleValues);
% Apply filter
filteredRawData = filtfilt(b, a, restBefore');
powerDat = filteredRawData';

% PLOT EXAMPLE
plotExample = 0;
if plotExample
    figure;
    subplot(3,1,1);
    plot(powerDat(1,:), 'k');
    xlim([87424 91709]);
    subplot(3,1,2);
    plot(powerDat(3,:), 'k');
    xlim([87424 91709]);
    subplot(3,1,3);
    plot(powerDat(7,:), 'k');
    xlim([87424 91709]);
    ckSTIM_saveFig(['subj4u3s-resting-filter-example'],5,15,300,'',1,[10,11],[]);
    close all;
elseif plotExample ==2
    % plot various power band for demonstration
    figure;
    freq = [1,4; 4 8; 8 12; 12 25; 25 40; 40 170];
    subplot(7,1,1)
    plot(restBefore(1,87424:91709), 'k');
    xlim([87424 91709]); box off;
    
    for f = 1:6
        subplot(7,1,f+1)
        if f == 1
            N   = 4;  % Order
        else
            N = 8;
        end
        
        Fc1 = freq(f,1);  % First Cutoff Frequency
        Fc2 = freq(f,2);  % Second Cutoff Frequency
        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, FS);
        Hd = design(h, 'butter');
        [b,a] = sos2tf(Hd.sosMatrix, Hd.ScaleValues);
        % Apply filter
        filteredRawData = filtfilt(b, a, restBefore');
        powerDat = filteredRawData';

        plot(powerDat(1,:), 'k');
        xlim([87424 91709]);box off;

        ckSTIM_saveFig(['subj4u3s-resting-filter-example-allpower'],5,15,300,'',1,[10,11],[]);
        close all;
    end


end

% PLV
phis = angle(hilbert(powerDat(:,:)'));
pairs = nchoosek(1:numChan,2);           % all combinations of channels in groups of two
plv = zeros(numChan,numChan);
for pr = 1 : length(pairs)
    plv(pairs(pr,1), pairs(pr,2)) = abs(sum(exp(1i* ( phis(:,pairs(pr,1)) - phis(:,pairs(pr,2)) ))))/length(phis);
end
plv = plv + plv';                   % we take its transpose to create the connectivity matrix

% PEARSON
pea = corrcoef(powerDat(:,:)');
pea = abs(pea - diag(diag(pea)));

plvPre = plv; peaPre = pea;

% GET DISTANCE RESIDUALS
Boolind = tril(repmat(true,numChan,numChan),-1);
Boolind = Boolind(:);
conVec = plvPre(Boolind);
dVec = d(Boolind);

regY = log(conVec./(1-conVec));
regX = 1./exp(dVec/10);
mdl = fitlm(regX, regY);
resVec = mdl.Residuals.Raw;
conSqr = zeros(numChan);
conSqr(Boolind) = resVec;
conSqr = conSqr + conSqr';

plvPreRes = conSqr; 

conVec = peaPre(Boolind);
dVec = d(Boolind);

regY = log(conVec./(1-conVec));
regX = 1./exp(dVec/10);
mdl = fitlm(regX, regY);
resVec = mdl.Residuals.Raw;
conSqr = zeros(numChan);
conSqr(Boolind) = resVec;
conSqr = conSqr + conSqr';

peaPreRes = conSqr; 

% POST
% filter for powerband
filteredRawData = filtfilt(b, a, restAfter');
powerDat = filteredRawData';

% PLV
phis = angle(hilbert(powerDat(:,:)'));
pairs = nchoosek(1:numChan,2);           % all combinations of channels in groups of two
plv = zeros(numChan,numChan);
for pr = 1 : length(pairs)
    plv(pairs(pr,1), pairs(pr,2)) = abs(sum(exp(1i* ( phis(:,pairs(pr,1)) - phis(:,pairs(pr,2)) ))))/length(phis);
end
plv = plv + plv';                   % we take its transpose to create the connectivity matrix

% PEARSON
pea = corrcoef(powerDat(:,:)');
pea = abs(pea - diag(diag(pea)));

plvPost = plv; peaPost = pea;

% GET DISTANCE RESIDUALS
Boolind = tril(repmat(true,numChan,numChan),-1);
Boolind = Boolind(:);
conVec = plvPost(Boolind);
dVec = d(Boolind);

regY = log(conVec./(1-conVec));
regX = 1./exp(dVec/10);
mdl = fitlm(regX, regY);
resVec = mdl.Residuals.Raw;
conSqr = zeros(numChan);
conSqr(Boolind) = resVec;
conSqr = conSqr + conSqr';

plvPostRes = conSqr; 

conVec = peaPost(Boolind);
dVec = d(Boolind);

regY = log(conVec./(1-conVec));
regX = 1./exp(dVec/10);
mdl = fitlm(regX, regY);
resVec = mdl.Residuals.Raw;
conSqr = zeros(numChan);
conSqr(Boolind) = resVec;
conSqr = conSqr + conSqr';

peaPostRes = conSqr; 

%% PROCESS PRE
datTimePre = 0:(1/FS):length(restBefore)/FS; datTimePre = datTimePre(1:end-1);
restingFT_pre = [];
restingFT_pre.trial{1} = restBefore;
restingFT_pre.time{1} = datTimePre;
restingFT_pre.label = ChannelList;
cfgir = [];
cfgir.trl = zeros(size(restingFT_pre.trial,2), 3);
cfgir.trl(1,1) = 1;
cfgir.trl(1,2) = length(restingFT_pre.time{1});
for i = 2:size(restingFT_pre.trial,2)
    cfgir.trl(i,1) = 1 + length(restingFT_pre.time{1})*(i-1);
    cfgir.trl(i,2) = length(restingFT_pre.time{1}) + length(restingFT_pre.time{1})*(i-1);
end
epochRestingPre = ft_redefinetrial(cfgir,restingFT_pre);

% chunk into chunks
cfgi = [];
cfgi.length = 2; 
epochRestingPre_chunked = ft_redefinetrial(cfgi, epochRestingPre);
for t = 1:length(epochRestingPre_chunked.trial)
    epochRestingPre_chunked.time{t} = epochRestingPre_chunked.time{1};
end

% GET VOLTAGE CORRELATIONS
cfgi = [];
cfgi.keeptrials = 'yes';
volDat = ft_timelockanalysis(cfgi,epochRestingPre_chunked);

volPea = zeros(size(volDat.trial,1), length(ChannelList),length(ChannelList));
for t = 1:size(volDat.trial,1)
    pea = corrcoef(squeeze(volDat.trial(t,:,:))');
    pea = abs(pea - diag(diag(pea)));
    volPea(t,:,:) = pea;
end
volPeaPre = squeeze(nanmean(volPea));

% calculate measures
cfgi            = [];
cfgi.output     = 'fourier';
cfgi.method     = 'mtmfft';
cfgi.taper = 'dpss';
%cfgi.foi     = [1:30, 35:5:120];
cfgi.foilim = [0 120];
cfgi.tapsmofrq  = 4;
cfgi.keeptrials = 'yes';
cfgi.channel    = 'all';
freqfourier    = ft_freqanalysis(cfgi, epochRestingPre_chunked);

% frequency bands of interest
freqBands        = [1 4; 5 13;12 25;25 50;70 120];

% calculate coh
epochRestingPre_coh = zeros(length(freqBands), length(ChannelList),length(ChannelList));
epochRestingPre_cohRes = zeros(length(freqBands), length(ChannelList),length(ChannelList));
cfgi             = [];
cfgi.method      = 'coh';
coh_ft           = ft_connectivityanalysis(cfgi, freqfourier);
for k = 1:size(freqBands,1)
    epochRestingPre_coh(k,:,:) = nanmean(coh_ft.cohspctrm(:,:,freqfourier.freq>=freqBands(k,1) & freqfourier.freq<=freqBands(k,2)),3); 
end

% GET DISTANCE RESIDUALS
Boolind = tril(repmat(true,numChan,numChan),-1);
Boolind = Boolind(:);
for k = 1:size(freqBands,1)
    datUse = epochRestingPre_coh(k,:,:);
    conVec = datUse(Boolind);
    dVec = d(Boolind);
    
    regY = log(conVec./(1-conVec));
    regX = 1./exp(dVec/10);
    mdl = fitlm(regX, regY);
    resVec = mdl.Residuals.Raw;
    conSqr = zeros(numChan);
    conSqr(Boolind) = resVec;
    conSqr = conSqr + conSqr';

    epochRestingPre_cohRes(k,:,:) = conSqr;
end

% calculate PLV
epochRestingPre_plv = zeros(length(freqBands), length(ChannelList),length(ChannelList));
epochRestingPre_plvRes = zeros(length(freqBands), length(ChannelList),length(ChannelList));
cfgi             = [];
cfgi.method      = 'plv';
plv_ft           = ft_connectivityanalysis(cfgi, freqfourier);
for k = 1:size(freqBands,1)
    epochRestingPre_plv(k,:,:) = nanmean(plv_ft.plvspctrm(:,:,freqfourier.freq>=freqBands(k,1) & freqfourier.freq<=freqBands(k,2)),3); 
end

% GET QEXP COMMUNICABILITY
epochRestingPre_qexp = zeros(length(freqBands), length(ChannelList),length(ChannelList));
for k = 1:size(freqBands,1)
    Wqexp = squeeze(epochRestingPre_plv(k,:,:));
    epochRestingPre_qexp(k,:,:) = getCommunicability(Wqexp,1,0);
end

% GET DISTANCE RESIDUALS
Boolind = tril(repmat(true,numChan,numChan),-1);
Boolind = Boolind(:);
for k = 1:size(freqBands,1)
    datUse = epochRestingPre_plv(k,:,:);
    conVec = datUse(Boolind);
    dVec = d(Boolind);
    
    regY = log(conVec./(1-conVec));
    regX = 1./exp(dVec/10);
    mdl = fitlm(regX, regY);
    resVec = mdl.Residuals.Raw;
    conSqr = zeros(numChan);
    conSqr(Boolind) = resVec;
    conSqr = conSqr + conSqr';

    epochRestingPre_plvRes(k,:,:) = conSqr;
end

%% PROCESS POST
datTimePost = 0:(1/FS):length(restAfter)/FS; datTimePost = datTimePost(1:end-1);
restingFT_post = [];
restingFT_post.trial{1} = restAfter;
restingFT_post.time{1} = datTimePost;
restingFT_post.label = ChannelList;
cfgir = [];
cfgir.trl = zeros(size(restingFT_post.trial,2), 3);
cfgir.trl(1,1) = 1;
cfgir.trl(1,2) = length(restingFT_post.time{1});
for i = 2:size(restingFT_post.trial,2)
    cfgir.trl(i,1) = 1 + length(restingFT_post.time{1})*(i-1);
    cfgir.trl(i,2) = length(restingFT_post.time{1}) + length(restingFT_post.time{1})*(i-1);
end
epochRestingPost = ft_redefinetrial(cfgir,restingFT_post);

% chunk into chunks
cfgi = [];
cfgi.length = 2; 
epochRestingPost_chunked = ft_redefinetrial(cfgi, epochRestingPost);
for t = 1:length(epochRestingPost_chunked.trial)
    epochRestingPost_chunked.time{t} = epochRestingPost_chunked.time{1};
end

% calculate measures
cfgi            = [];
cfgi.output     = 'fourier';
cfgi.method     = 'mtmfft';
cfgi.taper = 'dpss';
%cfgi.foi     = [1:30, 35:5:120];
cfgi.foilim = [0 120];
cfgi.tapsmofrq  = 4;
cfgi.keeptrials = 'yes';
cfgi.channel    = 'all';
freqfourier    = ft_freqanalysis(cfgi, epochRestingPost_chunked);

% frequency bands of interest
freqBands        = [1 4; 5 13;12 25;25 50;70 120];

% calculate coh
epochRestingPost_cohRes = zeros(length(freqBands), length(ChannelList),length(ChannelList));
epochRestingPost_coh = zeros(length(freqBands), length(ChannelList),length(ChannelList));
cfgi             = [];
cfgi.method      = 'coh';
coh_ft           = ft_connectivityanalysis(cfgi, freqfourier);
for k = 1:size(freqBands,1)
    epochRestingPost_coh(k,:,:) = nanmean(coh_ft.cohspctrm(:,:,freqfourier.freq>=freqBands(k,1) & freqfourier.freq<=freqBands(k,2)),3); 
end

% GET DISTANCE RESIDUALS
Boolind = tril(repmat(true,numChan,numChan),-1);
Boolind = Boolind(:);
for k = 1:size(freqBands,1)
    datUse = epochRestingPost_coh(k,:,:);
    conVec = datUse(Boolind);
    dVec = d(Boolind);
    
    regY = log(conVec./(1-conVec));
    regX = 1./exp(dVec/10);
    mdl = fitlm(regX, regY);
    resVec = mdl.Residuals.Raw;
    conSqr = zeros(numChan);
    conSqr(Boolind) = resVec;
    conSqr = conSqr + conSqr';

    epochRestingPost_cohRes(k,:,:) = conSqr;
end

% calculate PLV
epochRestingPost_plvRes = zeros(length(freqBands), length(ChannelList),length(ChannelList));
epochRestingPost_plv = zeros(length(freqBands), length(ChannelList),length(ChannelList));
cfgi             = [];
cfgi.method      = 'plv';
cfgi.trials      = t;
plv_ft           = ft_connectivityanalysis(cfgi, freqfourier);
for k = 1:size(freqBands,1)
    epochRestingPost_plv(k,:,:) = nanmean(plv_ft.plvspctrm(:,:,freqfourier.freq>=freqBands(k,1) & freqfourier.freq<=freqBands(k,2)),3); 
end

% GET DISTANCE RESIDUALS
Boolind = tril(repmat(true,numChan,numChan),-1);
Boolind = Boolind(:);
for k = 1:size(freqBands,1)
    datUse = epochRestingPost_plv(k,:,:);
    conVec = datUse(Boolind);
    dVec = d(Boolind);
    
    regY = log(conVec./(1-conVec));
    regX = 1./exp(dVec/10);
    mdl = fitlm(regX, regY);
    resVec = mdl.Residuals.Raw;
    conSqr = zeros(numChan);
    conSqr(Boolind) = resVec;
    conSqr = conSqr + conSqr';

    epochRestingPost_plvRes(k,:,:) = conSqr;
end

%% GET NETWORK MEASURES
% SEE http://complexity.es/school/neuroscience

% THRESHOLD NETWORK
numFeat = 2;
st = zeros(numChan,numFeat);
cw = zeros(numChan,numFeat);
ec = zeros(numChan,numFeat);
deg = zeros(numChan,numFeat);
bw = zeros(numChan,numFeat);
wdz = zeros(numChan,numFeat);
part = zeros(numChan,numFeat);

for feat = 1:numFeat
    freqUse = 2;
    if feat == 1
        featUse = squeeze(epochRestingPre_coh(freqUse,:,:)); 
    elseif feat ==2
        featUse = squeeze(epochRestingPre_plv(freqUse,:,:)); 
    end
    
    W = threshold_proportional(featUse, 0.6);   % look at this! "W" is now our matrix with threshold
    
    % GET NODE STRENGTH
    st(:,feat) = sum(W,2);
    
    % GET WEIGHTED CLUSTERING
    cw(:,feat) = clustering_coef_wu(W);
    
    % GET EIGENVECTOR CENTRALITY
    ec(:,feat) = eigenvector_centrality_und(W);

    % BETWEENNESS CENTRALITY
    bw(:,feat) = betweenness_bin(W);

    % WITHIN MODULE Z & PARTICIPATION
    [M,~] = community_louvain(W);
    wdz(:,feat) = module_degree_zscore(W,M);
    part(:,feat) = participation_coef(W,M);
    
    % GET DEGREES
    deg(:,feat) = degrees_und(W);
end

restPre.st = st;
restPre.cw = cw;
restPre.ec = ec;
restPre.deg = deg;
restPre.dist = d;
restPre.part = part;
restPre.bw = bw;
restPre.wdz = wdz;
restPre.epochRestingPre_coh =epochRestingPre_coh;
restPre.epochRestingPre_plv =epochRestingPre_plv;
restPre.epochRestingPre_cohRes =epochRestingPre_cohRes;
restPre.epochRestingPre_plvRes =epochRestingPre_plvRes;
restPre.qexp = epochRestingPre_qexp;
restPre.volPea = volPeaPre;
restPre.plv = plvPre;
restPre.pea = peaPre;
restPre.plvRes = plvPreRes;
restPre.peaRes = peaPreRes;


% POST
numFeat = 2;
st = zeros(numChan,numFeat);
cw = zeros(numChan,numFeat);
ec = zeros(numChan,numFeat);
deg = zeros(numChan,numFeat);
bw = zeros(numChan,numFeat);
wdz = zeros(numChan,numFeat);
part = zeros(numChan,numFeat);

for feat = 1:numFeat
    freqUse = 2;
    if feat == 1
        featUse = squeeze(epochRestingPost_coh(freqUse,:,:)); 
    elseif feat ==2
        featUse = squeeze(epochRestingPost_plv(freqUse,:,:)); 
    end
    
    W = threshold_proportional(featUse, 0.6);   % look at this! "W" is now our matrix with threshold
    
    % GET NODE STRENGTH
    st(:,feat) = sum(W,2);
    
    % GET WEIGHTED CLUSTERING
    cw(:,feat) = clustering_coef_wu(W);
    
    % GET EIGENVECTOR CENTRALITY
    ec(:,feat) = eigenvector_centrality_und(W);

    % BETWEENNESS CENTRALITY
    bw(:,feat) = betweenness_bin(W);

    % WITHIN MODULE Z & PARTICIPATION
    [M,~] = community_louvain(W);
    wdz(:,feat) = module_degree_zscore(W,M);
    part(:,feat) = participation_coef(W,M);
    
    % GET DEGREES
    deg(:,feat) = degrees_und(W);
end

restPost.st = st;
restPost.cw = cw;
restPost.ec = ec;
restPost.deg = deg;
restPost.dist = d;
restPost.part = part;
restPost.bw = bw;
restPost.wdz = wdz;
restPost.epochRestingPost_coh =epochRestingPost_coh;
restPost.epochRestingPost_plv =epochRestingPost_plv;
restPost.epochRestingPost_cohRes =epochRestingPost_cohRes;
restPost.epochRestingPost_plvRes =epochRestingPost_plvRes;
restPost.plv = plvPost;
restPost.pea = peaPost;
restPost.plvRes = plvPostRes;
restPost.peaRes = peaPostRes;

outdata.restPre = restPre;
outdata.restPost = restPost;

save([dtype], '-struct', 'outdata', '-v7.3');

