function [Syn,Redncy] = discreteEntropySynRed(S1,S2)

%Net redundancy estimator across pairwise binarised signals

%%
%Input
%Binarised signals (S1/S2)

%Ouput:
%Redncy - Redundancy
%Syn - Synergy

%% 

%S1
[~,ia,ic] = unique(S1(:),'rows');
stateCount = histcounts(ic,numel(ia));
prob_S1 = stateCount./sum(stateCount);
localEntrpy_S1 = -1.*log2(prob_S1(ic)); %[bits]

%S2
[~,ia,ic] = unique(S2(:),'rows');
stateCount = histcounts(ic,numel(ia));
prob_S2 = stateCount./sum(stateCount);
localEntrpy_S2 = -1.*log2(prob_S2(ic)); %[bits]

%S1&S2
[~,ia,ic] = unique([S1(:) S2(:)],'rows');
stateCount = histcounts(ic,numel(ia));
prob_S1S2 = stateCount./sum(stateCount);
localEntrpy_S1S2 = -1.*log2(prob_S1S2(ic)); %natural log [bits]

%Red 
localRedncy = min([localEntrpy_S1(:)'; localEntrpy_S2(:)']);
Redncy = mean(localRedncy);

%Individual information
%U1
U1 = mean(localEntrpy_S1) - Redncy;
%U2
U2 = mean(localEntrpy_S2) - Redncy;
%Synergy
Syn = mean(localEntrpy_S1S2) - U1 - U2 - Redncy;

end

