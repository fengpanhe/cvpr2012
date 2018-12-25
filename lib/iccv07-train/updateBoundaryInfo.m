function [bndinfo2, err] = updateBoundaryInfo(bndinfo, result)

%% Initialize

ne = bndinfo.ne;
nde = ne*2;
regions = result.regions;
nr = numel(regions);
chains = result.edgeletChain;
nc = numel(chains);
spLR = result.chainLR;

e2chain = zeros(nde, 1);
for k = 1:numel(chains)
    e2chain(chains{k}) = k;
end


%% Order chains and split when disconnected

enext = bndinfo.edges.adjacency;
eprev = cell(nde, 1);
for k1 = 1:numel(enext)
    enext{k1} = enext{k1}(:)';
    for k2 = enext{k1}(:)'
        eprev{k2}(end+1) = k1;
    end
end

ndc = numel(chains);  % number of directed chains
nc = ndc / 2;
nc2 = 0;
firste = cell(nc, 1);
for k = 1:nc
    if numel(chains{k})==1
        firste{k} = 1;
    else    
        firste{k} = [];    
        for k2 = 1:numel(chains{k})
            e = chains{k}(k2);
            if ~any(e2chain(eprev{e})==k)
                firste{k}(end+1) = k2;
            end        
        end
    end    
    if isempty(firste{k})
        firste{k} = 1;
    end         
   	nc2 = nc2 + numel(firste{k});
end
chains2 = cell(1, nc2*2);
spLR2 = zeros(nc2, 2);
e2chain2 = zeros(bndinfo.ne*2, 1);

nc2max = nc2;
nc2 = 0;
for k = 1:nc        
    oldchain = chains{k};            
    for k2 = 1:numel(firste{k})
        subchain = zeros(1, numel(oldchain));
        clen = 1;
        subchain(clen) = firste{k}(k2);
        
        nc2 = nc2+1;        
        done = false;
        while ~done
            done = true;
            curre = oldchain(subchain(clen));            
            for k3 = enext{curre}  
                if e2chain(k3)==k
                    clen = clen+1;
                    subchain(clen) = find(oldchain==k3);
                    e2chain(k3) = 0;
                    done = false;
                    break;
                end
            end
        end
        subchain = subchain(1:clen);
        
        chains2{nc2} = oldchain(subchain);
        e2chain2(chains2{nc2}) = nc2;
        oldchain(subchain) = 0;
        spLR2(nc2, :) = spLR(k, :);
        
        chains2{nc2+nc2max} = ...
            mod(chains2{nc2}(end:-1:1)+bndinfo.ne-1, bndinfo.ne*2)+1;
        e2chain2(chains2{nc2+nc2max}) = nc2+nc2max;
        %spLR2(nc2+nc2max, :) = spLR(k, [2 1]);
        
    end
end
nc2 = nc2*2;        

%% Get edglet adjacency and other information

nc = nc2/2;
ejunctions = bndinfo.edges.junctions;
ejunctions = [ejunctions ; ejunctions(:, [2 1])];
ejunctions2 = zeros(nc2, 2);
thetaDirected = bndinfo.edges.thetaDirected;
thetaDirected = [thetaDirected ; thetaDirected+pi];

enext2 = cell(nc2, 1);
theta2 = zeros(nc2, 1);
thetaStart = zeros(nc2, 1);
thetaEnd = zeros(nc2, 1);
for k = 1:nc2
    enext2{k} = e2chain2(enext{chains2{k}(end)});
    enext2{k} = enext2{k}(enext2{k}>0);
    ejunctions2(k, :) = [ejunctions(chains2{k}(1), 1) ejunctions(chains2{k}(end), 2)];
    jpos1 = bndinfo.junctions.position(ejunctions(k, 1), :);
    jpos2 = bndinfo.junctions.position(ejunctions(k, 2), :);
    theta2(k) = atan2(-(jpos2(2)-jpos1(2)), jpos2(1)-jpos1(1));
    if isfield(bndinfo.edges, 'thetaStart')
        thetaStart(k) = bndinfo.edges.thetaStart(chains2{k}(1));
        thetaEnd(k) = bndinfo.edges.thetaEnd(chains2{k}(end));
    else
        thetaStart(k) = thetaDirected(chains2{k}(1));
        thetaEnd(k) = thetaDirected(chains2{k}(end));
    end
end


%% Create data structure

bndinfo2.nseg = nr;
bndinfo2.edges.chains = chains2;
bndinfo2.edges.adjacency = enext2;
bndinfo2.edges.junctions = ejunctions2;
bndinfo2.edges.spLR = spLR2;
bndinfo2.edges.thetaDirected = theta2;
bndinfo2.edges.thetaStart = thetaStart;
bndinfo2.edges.thetaEnd =thetaEnd;
bndinfo2.junctions.position = bndinfo.junctions.position;
bndinfo2.ne = nc;
bndinfo2.imname = bndinfo.imname;
bndinfo2.imsize = bndinfo.imsize;
bndinfo2.names = bndinfo.names;
bndinfo2.nj = bndinfo.nj;
bndinfo2.type = bndinfo.type;
bndinfo2.regions2sp = regions;
bndinfo2.origne = bndinfo.ne;
bndinfo2.orignseg = bndinfo.nseg;

splab = zeros(bndinfo.nseg, 1);
for k = 1:numel(regions)
    splab(regions{k}) = k;
end
bndinfo2.wseg = uint16(splab(bndinfo.wseg));

stats = regionprops(bndinfo.wseg, 'Area');
area = cat(1, stats(:).Area);
bndinfo2.spArea = area;
[bndinfo2.labels, err] = transferLabels(bndinfo.labels, regions, area);
bndinfo2 = processGtBoundaryLabels(bndinfo2);

bndinfo2 = orderfields(bndinfo2);
