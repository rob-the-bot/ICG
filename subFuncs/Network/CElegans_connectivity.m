
load('chemAlligned_male_herm.mat')
load('gapJuncAlligned_male_herm.mat')

%%
hChem = double(hChem>0);
mChem = double(mChem>0);
hGJ = double(hGJ>0);
mGJ = double(mGJ>0);
%%

nonEmptyCellIDS = nan(1,numel(cellIDs));

for cc = 1:numel(cellIDs)
    
    nonEmptyCellIDS(cc) = ~isempty(cellIDs{cc});
end

nonEmptyCell = find(nonEmptyCellIDS);




%% Contrast connectivity of functional ICG pairings

connectivityhGJ = cell(1,nN);
connectivitymGJ = cell(1,nN);
connectivityhChem = cell(1,nN);
connectivitymChem = cell(1,nN);


for ICGlevel = 2:numel(activityICG)-1
    ICGlevel
    
    allhgj =[];
    allmgj =[];
    
    allhChem =[];
    allmChem =[];
    
    for pp = 1:size(outPairOGID{ICGlevel},1)
        
        CGpairs = outPairOGID{ICGlevel}(pp,:);
        
        
        possPairs = nchoosek(CGpairs,2);
        
        for ll = 1:size(possPairs,1)
            
            %GJ are symmetric so only need one pair
            preCell = cellIDs{possPairs(ll,1)};
            postCell = cellIDs{possPairs(ll,2)};
                 
            for prC = 1:numel(preCell)
                for poC = 1:numel(postCell)
                    
                    prCTemp = find(preCell{prC} == cGJ);
                    poCTemp = find(postCell{poC} == cGJ);
                    
                    if isempty(prCTemp)
                        continue
                    elseif isempty(poCTemp)
                        continue
                    end
                
                    
                    allhgj = [allhgj hGJ(prCTemp,poCTemp)];
                    allmgj = [allmgj mGJ(prCTemp,poCTemp)];
                    
                    
                    
                end
                
            end
            
            %Chemical pre/post            
            for prC = 1:numel(preCell)
                for poC = 1:numel(postCell)
                    
                    prCTemp = find(preCell{prC} == cpre);
                    poCTemp = find(postCell{poC} == cpost);
                    
                    if isempty(prCTemp)
                        continue
                    elseif isempty(poCTemp)
                        continue
                    end
                
                    allhChem = [allhChem hChem(prCTemp,poCTemp)];
                    allmChem = [allmChem mChem(prCTemp,poCTemp)];
                    
                end
            end
            
            %Chemical pre/post swapped
            preCell = cellIDs{possPairs(ll,2)};
            postCell = cellIDs{possPairs(ll,1)};
            
            for prC = 1:numel(preCell)
                for poC = 1:numel(postCell)
                    
                    prCTemp = find(preCell{prC} == cpre);
                    poCTemp = find(postCell{poC} == cpost);
                    
                    if isempty(prCTemp)
                        continue
                    elseif isempty(poCTemp)
                        continue
                    end
                 
                    allhChem = [allhChem hChem(prCTemp,poCTemp)];
                    allmChem = [allmChem mChem(prCTemp,poCTemp)];
                    
                end
            end
            
         
        end
        
    end
    
    connectivityhGJ{ICGlevel} = allhgj;
    connectivitymGJ{ICGlevel} = allmgj;
    
    connectivityhChem{ICGlevel} = allhChem;
    connectivitymChem{ICGlevel} = allmChem;
end

%

hgj_ICG = [];
mgj_ICG = [];
hChem_ICG = [];
mChem_ICG = [];


for ICGlevel = 2:numel(activityICG)-1
    
    allhgj = connectivityhGJ{ICGlevel};
    allmgj = connectivitymGJ{ICGlevel} ;
    
    allhChem = connectivityhChem{ICGlevel};
    allmChem = connectivitymChem{ICGlevel};
    
    hgj_ICG = [hgj_ICG mean(allhgj)];
    mgj_ICG = [mgj_ICG mean(allmgj)];
    hChem_ICG = [hChem_ICG mean(allhChem)];
    mChem_ICG = [mChem_ICG mean(allmChem)];
    
end

