
%% 1. perfectly correlated


allData = randn(1,1e4);
allData = repmat(allData,2048,1);

%If interested to explore adding noise
%allData = allData+0.1.*randn(2048,1e4);

%[activityICG,outPairID] =  ICG(allData);
%[VarAct,KurAct,TscAct, VarK] = calcVKT(activityICG,fs);

%% 2. Independent/uncorrelated

allData = randn(2048,1e5);

%[activityICG,outPairID] =  ICG(allData);
%[VarAct,KurAct,TscAct, VarK] = calcVKT(activityICG,fs);



%% Data matched

%load in some Data to match (
CellAct = dat;
%How much do you want to match
allData = CellAct(1:2048,:);


%% Time shuffling within in modules 

[activityICG,outOGID] =  ICG(allData);
[VarAct,KurAct,TscAct, VarK] = calcVKT(activityICG,fs);


figure(1)
subplot(2,2,3)
semilogx(VarK(1:end-1),VarAct,'r.-','MarkerSize',5)
hold on


% ICG level at which to shuffle within e.g., (2/8/32)

for msiz = [2 4 6]
timeShufact = allData;

modLVL  = msiz;

for mm = 1:size(outOGID{modLVL},1)
    
    permT = randi(size(allData,2),1);
    
    
    for nn = 1:size(outOGID{modLVL},2)  
        nID = outOGID{modLVL}(mm,nn);
        timeShufact(nID,:) = [timeShufact(nID,permT:end)  timeShufact(nID,1:permT-1)];
        
    end
    
end

[activityICG,outOGID] =  ICG(allData);
[VarAct,KurAct,TscAct, VarK] = calcVKT(activityICG,fs);


figure(1)
subplot(2,2,3)
semilogx(VarK(1:end-1),VarAct,'r.-','MarkerSize',5)
hold on

end


%% Variance matched (timescale destroyed)

%covariance matrix to match
Sigma = cov(allData');

% Matched modular static
% 4 /32/128

for msiz = [3 6 8]


modLVL  = msiz;

%build the cov matrix (e.g., 2048 by 2048)
covMod_Binary = zeros(2048);


for mm = 1:size(outOGID{modLVL},1)
    
    C = nchoosek(outOGID{modLVL}(mm,:),2);
    
    for pp = 1:size(C,1)
    covMod_Binary(C(pp,1),C(pp,2)) = 1;
    end
    
    
end

%add transpose for opposite entries
covMod_Binary = covMod_Binary+covMod_Binary';
%add eyes for self terms
covMod_Binary = covMod_Binary+eye(2048);


A_PD_Mod = Sigma.*covMod_Binary;

[V,D] = eig(A_PD_Mod);   % Calculate the eigendecomposition of your matrix (A = V*D*V')    
% where "D" is a diagonal matrix holding the eigenvalues of your matrix "A"
d= diag(D);      % Get the eigenvalues in a vector "d" 
d(d <= 1e-7) = 1e-7;  % Set any eigenvalues that are lower than threshold 1e-7 to a fixed non-zero small value 
D_c = diag(d);        % Build the updated diagonal matrix "D_c"
A_PD_Mod = V*D_c*V';

X = randn(5e3,2048);
X = bsxfun(@minus, X, mean(X));
X = X * inv(chol(cov(X)));
X = X * chol(A_PD_Mod);


% [activityICG,outOGID] =  ICG(allData);
% [VarAct,KurAct,TscAct, VarK] = calcVKT(activityICG,fs);


end

yline(3)


%% Poisson null independent


dt = 1/1000;
t = 0:dt:100;
nBins = numel(t);


srate = 30; 
spikeTimeStep = t;
spikeTimeStepHist = spikeTimeStep(1)-1/(2*srate):1/srate:spikeTimeStep(end)+1/(2*srate);
spikeTimeSubStep = spikeTimeStep(1):1/srate:spikeTimeStep(end);

nN = 8*1024;
%
% Independent simulations

spikeRate = nan(nN,numel(spikeTimeSubStep));

for nn = 1:nN
    
    if ~mod(nn,200)
        nn/nN
    end
    
   
    GW = 1.*randn(1,nBins)+7;
    
    %Generate the spikes
    spike2 = rand(nBins,1) < GW(:).*dt;
    spike2 = find(spike2);
    spkTim = spike2.*dt;
    
   
    %Now go through and calc rates
    Ncnt = histcounts(spkTim,spikeTimeStepHist);
    
    spikeRate(nn,:) = Ncnt;
    
end






%% Poisson null correlated

%How do you want to generate a global signal
smoothBM = @(x,n) conv(x,gausswin(n)/sum(gausswin(n)),'same');

hSlp = 1;
G = abs(wfbm(hSlp,nBins));

% 
G35 = 3.5.*G./mean(G);

spikeRateCorrSparse = nan(nN,numel(spikeTimeSubStep));

for nn = 1:nN
        
    %Generate the spikes
    spike2 = rand(nBins, 1) < G35(:).*dt;
    spike2 = find(spike2);
    spkTim = spike2.*dt;
    
   
    %Now go through and calc rates
    Ncnt = histcounts(spkTim,spikeTimeStepHist);
    
    spikeRateCorrSparse(nn,:) = Ncnt;
    
end

