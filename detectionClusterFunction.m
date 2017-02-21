function detectionClusterFunction(N,cat,REF,time,dirname,eventindex,DNmax,TNmax,TNmax_end,maxsizeClust,maxTimeClus,fracD,Qlim,Nlim)

d = deg2km(distance(cat.lat(REF),cat.lon(REF),cat.lat,cat.lon));
timei = time;
% fucntion to detect cluster of seismicity after Chen and Shearer 2016 A new method to identify earthquake swarms applied to seismicity near the San Jacinto Fault, California
%% input
% N, max number of neighbour events analyzed
% cat, seismic catalog in structure format
% REF, reference event
% time, time of events in the catalog in matlab datenum
% dirname, otuput directory to save data
% eventindex, index of the events analyzed [1:k]
% DNmax, max sise of space window
% TNmax, max size of positive time window
% TNmax_end, max size of negative time window
% maxsizeClust, max size of cluster in space, after this size in the code
% stop and does not calculate the clustering values anymore
% maxTimeClus, as above but for maximum time
% fracD, fractal dimension fault
% Qlim, Cluster strenght limit used to be saved
% Nlim, minimum number of events to save a cluster
% Piero Poli ppoli@mit.edu

%% example:
% N=100
% TNmax = 5;
% DNmax=3;
% TNmax_end = 2;
% Nstart = 2;
% Qlim=2;
% fracD = 1.6;
% maxsizeClust=Inf;
% maxTimeClus=Inf;
% Nlim=5;
% %% read catalog
% catname='SA19902017USGSM4';
% [timeevent,late,lone,depth,mag,magtype] = ParseUSGS_Catalog([catname '.csv']);
% cat.matlabtime = datenum(timeevent);
% cat.mag = mag;
% cat.lat=late;
% cat.lon=lone;
% cat.depth = depth;
% dirname = output
% %% sort
% [cat.matlabtime,ix] = sort(cat.matlabtime);
% cat.mag = cat.mag(ix);
% cat.lon=cat.lon(ix);
% cat.lat=cat.lat(ix);
% cat.depth=cat.depth(ix);
% 
% time = datevec(cat.matlabtime);
% CNT = 1;
% eventindex = 1 : length(cat.matlabtime);
% NOW RUN THE CODE
% detectionClusterFunctionShearer(N,cat,REF,time,dirname,eventindex,DNmax,TNmax,TNmax_end,maxsizeClust,maxTimeClus,fracD,Qlim,Nlim)


%% CODE!!!

cluster.maxTimeClus=maxTimeClus;
cluster.maxsizeClust=maxsizeClust;
cluster.TNmax_end=TNmax_end;
cluster.TNmax=TNmax;
cluster.DNmax = DNmax;
cluster.Nlim=Nlim;
%% get time
dt = zeros(size(cat.lon));
for ic = 1 : length(cat.lon)
    
    dt(ic) =  etime(timei(ic,:),datevec(cat.matlabtime(REF,:)))./86400;
    
end
clear timei
ni = dt.*d.^fracD;
%% for catalog
% sort
id = eventindex;
[~,indx] = sort(ni);
mag = cat.mag(indx);
lon = cat.lon(indx);
lat = cat.lat(indx);
depth = cat.depth(indx);
abstime = cat.matlabtime(indx);
d = d(indx);
dt = dt(indx);
id = id(indx);
%remove
F = find(ni>=0);
dii = d(F);
dti=dt(F);

Q = zeros(N-2,1);Tmaxi=Q; Dmaxi=Q;indN=Q;Nev=Q;
Tev=Q;
Dev=Q;
cnt=1;

for in =  3 : N
    if length(dii) > in
        di = dii(1:in);
        ti = dti(1:in);
        %% estimate duration and size of cluster
        Dmax = (max(di));
        Tmax = max(ti)-min(ti);
        
        if Dmax > maxsizeClust
            
            break
        end
        
        if Tmax > maxTimeClus
            
            break
        end
        
        %% define out polygon
        xout = [-Tmax*TNmax -Tmax*TNmax Tmax*TNmax_end Tmax*TNmax_end];
        yout = [0 DNmax*Dmax DNmax*Dmax 0];
        [outside,~] = inpolygon(dt,d,xout,yout);
        %% define in polygon
        xin = [0 0 Tmax Tmax];
        yin = [0 Dmax Dmax 0];
        [inside,~] = inpolygon(dt,d,xin,yin);
        %% make estimate strength
        Q(in-2)=sum(inside)/((sum(outside)-sum(inside))+1);
        Dmaxi(in-2) = Dmax;
        Tmaxi(in-2) = Tmax;
        Nev(in-2) = sum(inside);
        Tev(in-2)=Tmax;
        Dev(in-2)=Dmax;
        indN(in-2)=in;
        %cnt = cnt+1;

    end
    clear di ti
end

%% get max Q
[maxQ,indQ] = max(Q);

%% check if Q is good
if maxQ > Qlim && Nev(indQ) > Nlim
    cluster.Q = Q;
    cluster.Qdist = Dev;
    cluster.Qtime = Tev;clear Tev Dev
    cluster.maxQind = indQ;
    Tmax = Tmaxi(indQ);
    Dmax = Dmaxi(indQ);
    
    
    xout = [-Tmax*TNmax -Tmax*TNmax Tmax*TNmax_end Tmax*TNmax_end];
    yout = [0 DNmax*Dmax DNmax*Dmax 0];
    [outside,~] = inpolygon(dt,d,xout,yout);
    
    xin = [0 0 Tmax Tmax];
    yin = [0 Dmax Dmax 0];
    [inside,~] = inpolygon(dt,d,xin,yin);

    %% save cluster information
    % outside
    cluster.dist = d(inside);
    cluster.relativetime = dt(inside);
    cluster.elat = lat(inside);
    cluster.elon = lon(inside);
    cluster.time = abstime(inside);
    cluster.mag = mag(inside);
    cluster.ref_event_id = (REF);
    cluster.events_id = id(inside);
    cluster.duration = Tmax;
    cluster.size = Dmax;
    cluster.strength = maxQ;
    cluster.depth = depth(inside);
    % outside
    cluster.dist_out = d(outside);
    cluster.relativetime_out = dt(outside);
    cluster.elat_out = lat(outside);
    cluster.elon_out = lon(outside);
    cluster.time_out = abstime(outside);
    cluster.mag_out = mag(outside);
    name = [dirname '/cluster' num2str(REF)];
    disp(name)
    save(name,'cluster')
    
end
end
