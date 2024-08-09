

%% power law connectivity network

gamvals = 1.5:0.25:4.5;

gamma = gamvals(gg);
N = 2^12;
x = randht(20*N,'xmin',1,'powerlaw',gamma);
x= floor(x);
w= x(x<N&x>1);

A = CL_generator(w);
CIJ = A;


CIJ = CIJ - diag(diag(CIJ));
bdIndx = find(sum(CIJ,1)==0 & (sum(CIJ,2)==0)');
CIJ(bdIndx,:) =[];
CIJ(:,bdIndx) =[];

CIJ = CIJ(1:N,1:N);
CIJ = full(CIJ);


%% Watts and Strogatz

grph = WattsStrogatz(N,k,beta);
prob = full(adjacency(grph));

%% make modular small world OS 

modvals = [0 1 2 3 4 5 6 7];
mx_lvl = 11;
N = 2^11;
E = 100;
[prob,CIJ,~] = makefractalCIJ(mx_lvl,E,sz_cl,1);


%% Spatial connectivity
nSqr = 64;

nsamp = 2;

x = 1:nsamp:nsamp*nSqr;
y = 1:nsamp:nsamp*nSqr;
[X,Y] = meshgrid(x,y);

%D = pdist([X(:) Y(:)]);
DX = pdist(X(:));
ZX = squareform(DX);
ZX(ZX>nSqr/2) = nSqr-ZX(ZX>nSqr/2);

DY = pdist(Y(:));
ZY = squareform(DY);
ZY(ZY>nSqr/2) = nSqr-ZY(ZY>nSqr/2);

Deuc = sqrt(ZX.^2+ZY.^2);


%Sweep slope param!
prob = exp(-1*slope*Deuc);



%%
function A = CL_generator(w)
% INPUT:    w = expected degrees vector
% OUTPUT:   A = binary adjacency matrix of a graph G in G(w)
    n = length(w);
    m = (dot(w,w)/sum(w))^2 + sum(w); 
    m = ceil(m/2);
    wsum = [0 ; cumsum(w(:))]; 
    wsum = wsum/wsum(end);
    I = discretize(rand(m,1),wsum);
    J = discretize(rand(m,1),wsum);
    A = sparse([I;J],[J;I],1,n,n);
    A = spones(A);
end