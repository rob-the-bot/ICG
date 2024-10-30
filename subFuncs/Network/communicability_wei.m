function   F = communicability_wei(CIJ)

%inputs TEST
%           CIJ    weighted connection matrix
%           i      row
%           j      column
%
%outputs
%           F      communicability
%=================================================

N = size(CIJ,1);
B = sum(abs(CIJ));
C = diag(B); 
D = C^(-(1/2)); % Normalization factor for each region
E = D * CIJ * D;
F = expm(E).*~eye(N);