function [VarAct,KurAct,TscAct, VarK] = calcVKT(activityICG,fs)

%Input:
%activityICG - ICG ouput
%fs - Sampling frequency

%Output:
%Variance/Timescale/Kurtosis across scales

%% Variance
VarAct = zeros(1,numel(activityICG));

for rr = 1:numel(activityICG)
    
    temp = activityICG{rr};
    
    %vtemp = mean(temp.^2,2)- mean(temp,2).^2;
    vtemp = var(temp,0,2);

    VarAct(rr) = mean(vtemp);  

end

numClus = size(VarAct,2)-1;
VarK = 2.^(0:numClus);

%% Timescale 

%In seconds
maxT = 3; 

nLags = ceil(maxT*(1/fs)); %fs is sampling rate of recording

TscAct = zeros(1,numel(activityICG));

for rr = 1:numel(activityICG)
    
    temp = activityICG{rr};

    autocorrData = nan(size(temp,1),nLags+1);

    for nN = 1:size(temp,1)

        lvlAct = temp(nN,:);
        [xCorrAct,~] = crosscorr(lvlAct,lvlAct,nLags);

        autocorrData(nn,:) = xCorrAct(nLags+1:end);
    end


    ogT = (0:nLags)./fs;

    newX = 0:0.01:maxT;
    newAC = interp1(ogT,mean(autocorrData),newX);

    %AUC, also works fitting an exp
    t_AUC = trapz(newX,newAC);

    TscAct(rr) = t_AUC;

end

%% Kurtosis

KurAct = zeros(1,numel(activityICG)-1);

for rr = 1:numel(activityICG)-1
    [numel(activityICG)-1 rr]

    activity = activityICG{rr};
  
    %Calculate pairwise correlations
    rho = corr(activity');

    C = triu(rho,1);
    clearvars rho

    %Grab just the upper triangle
    %Find the indices of each element
    upTriIndx = find(triu(true(size(C)),1)>0);
    %[allRIndx, allCIndx] = ind2sub(size(C),upTriIndx);

    %now it is a vector of correlation values
    C = C(upTriIndx);
    clearvars upTriIndx

    C = C(~isnan(C));

    KurAct(rr) = kurtosis(C,0);
    clearvars C

end

end