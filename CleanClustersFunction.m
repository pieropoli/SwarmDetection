%function CleanClustersFunction(cluster,outdir)
clear
clc
close
addpath ~/Seismic_Matlab_Functions/

%% paramters
inputdirectory = 'TestsInfSize';
outdirectory = 'SwamrsCleanSize300';
maxsize = 300;

mkdir(outdirectory)

outI = loaddatafromdir(inputdirectory);

for ii1 = [11 17]
    
    
    
    
    out = loaddatafromdir([inputdirectory '/' outI(ii1).name]);
    
    count = 1;
    for i1d = 1 : length(out)
        
        load([inputdirectory '/' outI(ii1).name '/' out(i1d).name])
        
        if  cluster.size< maxsize
            clusterout{count} = cluster;
            count = count + 1;
            
        end
        
        
        
    end
    
    clear cluster
    cluster = clusterout;
    test = ones(1,length(cluster));
    
    for i1 = 1 : length(cluster)
        
        clusref = cluster{i1};
        
        
        
        %% is taget a daughter in other groups?
        for i2 = 1 : length(cluster)
            if i2 ~= i1
                clustest = cluster{i2};
                
                F = intersect(clustest.events_id,clusref.ref_event_id);
                
                if isempty(F) == 0 && clustest.strength > clusref.strength
                    test(i1) = 0;
                    
                    disp('.')
                end
                
            end
        end
        
    end
    
    
    for i1 = 1 : length(cluster)
        if test(i1)>0
            
            
            clusref = cluster{i1};
            dauthers = clusref.events_id;
            Q(1) = clusref.strength;
            id(1) = i1;
            %% is taget a daughter in other groups?
            cnt = 2;
            for i2 = 1 : length(cluster)
                if i2 ~= i1 && test(i2)>0 && test(i1)>0
                    clustest = cluster{i2};
                    
                    F = intersect(dauthers,clustest.ref_event_id);
                    
                    if isempty(F) == 0 %&& clustest.strength > clusref.strength
                        Q(cnt)=clustest.strength;
                        id(cnt) = i2;
                        
                        disp('.')
                    end
                    
                end
            end
            
            
            
            if length(Q)>1
                
                IX = find(Q<max(Q));
                if isempty(IX) == 1;
                    IX = 2;
                end
                test(id(IX))=0;
                id(IX)
                
            end
            clear Q id
        end
        
    end
    
    
    
    %% retain only non duplicate clusters
    
    
    F = find(test==1);
    
    for iF = 1 : length(F)
        
        clusterselected{iF} = cluster{F(iF)};
        
        
    end
    
    
    %% SAVE
    name = [outdirectory '/' outI(ii1).name '_Clean_maxsize_' num2str(maxsize)];
    save(name,'clusterselected');
    
    clear clusterselected test id Q id cluster clusterout
end
