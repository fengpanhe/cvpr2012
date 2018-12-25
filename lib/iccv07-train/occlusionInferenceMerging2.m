function [lab, energy, result, valdata, dispim] = ...
    occlusionInferenceMerging2(pB, bndinfo, delta, DO_VAL)
%
% [lab, energy, edgeletChain, chainLR, regions] = ...
%   occlusionInferenceMerging(pB, pC, edgeind, bndinfo))
%
% Solves for y given boundary and continuity likelihoods and an adjacency
% list for each edge.  The MRF is constructed as a directed graph (that is,
% Eij(0,1)~=Eij(1,0)).  Each edgelet in the image forms two nodes in the
% graph (one for occlusions in each direction).  A link between these two
% nodes disallows them to both be 1.  Links between adjacent nodes in the
% directed graph encode the conditional term P(ej | ei).  Auxialliary
% junction nodes disallow edglets to begin from nowhere or end nowhere.
%
% Inputs:
%  pB(nedglets, 2): unary likelihoods
%  pC(npairs): p(ei = 1 | ej=1)
%  edgind(npairs, 2): edgelet indices
%  eadj{nedglets}: edge adjacency 
%

% Initial:
%   E*_i = [ -log[P(ei=0 | x)]       0 ]
%   E*_ji = [ 0                         0       
%            -log[P(ei=0|ej=1, x)]      -log[P(ei=1|ej=1, x)] ]
%
% Convert to below if not on a border junction:
%   K_i = min_j(E*_ji(2,2))
%   E_i = [E*_i(0)      K_i] 
%   E_ji = [0                       0       
%           E*_ji(2,1)-E*_i(0)      E*_ji(2,2)-K_i ]
%
% Junctions:
%   Infinite penalty for contour end and start, except at border junctions
%   Small penalty for T-junctions


% Inference
%
%  Edglet "off" energy = -log[P(ei=0 | x)]
%  Edglet "on" energy = min_j(-log[P(ei=1|ej=1, x)]) for j edglets on
%  Merging cost = sum_i(E(ei=1)-E(ei=0)) with sum over boundary between two  
%  regions.
%
%  Until merging cost is positive:
%    1) Merge two regions (Ri,Rj) with lowest cost (greatest drop in energy)
%    2) Update merging costs
%       a) if adj(Ri,Rk) and adj(Rj,Rk), C(Ri+Rj,Rk) = C(Ri,Rk)+C(Rj,Rk)
%       b) ek is adjacent to a removed edgelet, update its energy
%
%  Note: each edgelet chain could potentially face two directions if "on".
%  Set the cost based on the lowest energy direction.


%% Initialize edge and region information
% for each edgelet store:
%     adjacent edgelets, pairwise energy with adjacent edgelets, whether
%     the edgelet is on (all on initially), minimum "on" energy, "off"
%     energy
%
% for each region store:
%     which other regions it is adjacent to
%     the set of superpixels in each region
%
% for each pair of adjacent regions store:
%     edgelet list connecting them (edgelet ids forward and backwards)
%     current cost of merging

global DO_DISPLAY;

if ~exist('DO_VAL', 'var')
    DO_VAL = 0;    
end

if DO_VAL
    nval = 0;
    valstep = floor(bndinfo.nseg / 100);
    stats = regionprops(bndinfo.wseg, 'Area');
    areas = cat(1, stats(:).Area);
    segError = zeros(bndinfo.nseg, 1);
    mergeCost = zeros(bndinfo.nseg, 1);    
else
    valdata = [];
end

ne = bndinfo.ne;
nde = bndinfo.ne*2;
nsp = bndinfo.nseg;
wseg = bndinfo.wseg;

% spLR(i, [1 2]) gives the left and right, resp., 
spLR = bndinfo.edges.spLR;
spLR = [spLR ; spLR(:, [2 1])];

% edge adjacency: i-->enext{i}(j), eprev{i}(j)-->i 
enext = bndinfo.edges.adjacency;
eprev = cell(nde, 1);
for k1 = 1:numel(enext)
    enext{k1} = enext{k1}(:)';
    for k2 = enext{k1}(:)'
        eprev{k2}(end+1) = k1;
    end
end

% region adjacency (sparse): radj(i,j) gives index of edgelet chain for Ri
% occludes Rj
echain = num2cell(1:nde);  % edgelet chain between each pair of regions
rsp = num2cell(1:nsp);
validr = true(nsp, 1);
edgeOn = true(ne, 1);
evalid = repmat(edgeOn, [2 1]);
e2chain = (1:nde);

% unite borders between two superpixels that are split
reme = [];
spLRm = zeros([nsp nsp], 'uint16');
for k = 1:ne
    s1 = spLR(k, 1);  s2 = spLR(k,2);
    if spLRm(s1,s2)==0
        spLRm(s1,s2) = k;
    else
        reme(end+(1:2)) = [k k+ne];
        k2 = spLRm(s1, s2);
        echain{k2} = [echain{k2} echain{k}];
        echain{k2+ne} = [echain{k2+ne} echain{k+ne}];
        e2chain(echain{k2}) = k2;
        e2chain(echain{k2+ne}) = k2+ne;
    end
end
spLR(reme, :) = [];
echain(reme) = [];
e2chainshift = ones(numel(e2chain), 1);
e2chainshift(reme) = 0;
e2chainshift = cumsum(e2chainshift);
e2chain = e2chainshift(e2chain);
clear spLRm

% set initial energies
unary = -log(pB(:, 2));
Eon = -log(1-pB(:,1)+1E-10);
Eoff = -log(pB(:, 1)+1E-10); % energy for edge being off


%% Iteratively merge regions

mincost = -Inf;

Ecost = zeros(numel(echain), 1);
for k = 1:numel(echain)
    %Ecost(k) = sum(Eoff(echain{k}))-sum(Eon(echain{k}));
    Ecost(k) = max(log((1-pB(echain{k},1)) ./ pB(echain{k},1))) + delta;
end

iter = 0;
while (mincost < 0) %&& sum(validr)>100
    
    iter = iter + 1;    
        
    npair = numel(echain)/2;               
    
    % computeEnergy(echain, minlab+1)
    maxEcost = max([Ecost(1:npair) Ecost(npair+1:end)], [], 2);
    [mincost, minind] = min(maxEcost);    
    
    if mincost >= 0
        break;
    end
    
    if DO_VAL && mod(iter-1, valstep)==0        
        nval = nval + 1;
        nsp = max(wseg(:));
        segError(nval) = ...
            getSegmentationErrorFast(bndinfo.labels, nsp, rsp(validr), areas);
        mergeCost(nval) = mincost;
    end    
    
    if DO_DISPLAY && mod(iter, 500)==0 
%         Etotal = sum(Eoff(~edgeOn));  %0;
%         for c = 1:npair
%             Etotal = Etotal + min(sum(Eon(echain{c})), sum(Eon(echain{c+npair})));
%         end        
                
        %disp(num2str([iter Etotal sum(validr) sum(edgeOn) mincost]));
        
        displayRegions(rsp(validr), bndinfo.wseg, 1);        
        
        if 0
        nstr = num2str(100000 + iter);
        imwrite(tmpim, ['./tmp/merging' nstr(2:end) '.jpg']);        
        end
        
        %movframe(iter/25) = im2frame(tmpim);        
%        movframe(iter/25) = getframe(1);
    end    
    
    % remove influence from other edglets
    thischain = [echain{minind} echain{minind+npair}];
    for k = thischain        
        
        for k2 = enext{k}
            if ~(e2chain(k2)==minind || e2chain(k2)==minind+npair) && evalid(k2)
                Eij{k2}(eprev{k2}==k) = Inf;                                           
                if numel(eprev{k2})>1
                    Eon(k2) = min(Eij{k2});
                    if Eon(k2)==Inf
                        % this should only happen for a border 4-junction
                        % in which one edgelet is missing (due to ignoring
                        % edglets on the border)
                        Eon(k2) = unary(k2);  
                    end
                else
                    Eon(k2) = unary(k2);
                end
                c = e2chain(k2); 
%                Ecost(c) = sum(Eoff(echain{c}))-sum(Eon(echain{c}));
                Ecost(c) = max(log((1-pB(echain{c},1)) ./ pB(echain{c},1))) + delta;
            end
        end
    end    

    keep  =true(npair, 1);
    reme = [minind minind+npair];
    keep(minind) = false;
    
    r1 = spLR(minind, 1);
    r2 = spLR(minind, 2);
    
    newr = r1;
    
    % get neighbors to left and right of r1 and r2
    ind1L = find(spLR(:, 1)==r1);    
    ind2L = find(spLR(:, 1)==r2);  
    left1 = spLR(ind1L, 2);
    left2 = spLR(ind2L, 2);    
    ind1R = find(spLR(:, 2)==r1);
    ind2R = find(spLR(:, 2)==r2);

    spLR([ind1L ; ind2L], 1) = newr;    
    spLR([ind1R ; ind2R], 2) = newr;
    
    % unite any split borders that may arise
    indL = [ind1L ; ind2L]';  
    s1 = [left1 ; left2]';  
    for k1 = 1:numel(s1)
        for k2 = k1+1:numel(s1)
            if s1(k1)==s1(k2) && ~any(reme==indL(k2))
                
                %reme(end+1) = indL(k2);
                keep(mod(indL(k2)-1,npair)+1) = false;
                i1 = indL(k1);
                echain{i1} = [echain{i1} echain{indL(k2)}];
                e2chain(echain{i1}) = i1;
%                Ecost(i1) = sum(Eoff(echain{i1}))-sum(Eon(echain{i1}));
                Ecost(i1) = max(log((1-pB(echain{i1},1)) ./ pB(echain{i1},1))) + delta;
                
                i1 = mod(indL(k1)+npair-1, 2*npair)+1;
                i2 = mod(indL(k2)+npair-1, 2*npair)+1;                
                %reme(end+1) = i2;                
                echain{i1} = [echain{i1} echain{i2}];
                e2chain(echain{i1}) = i1;
%                Ecost(i1) = sum(Eoff(echain{i1}))-sum(Eon(echain{i1}));                
                Ecost(i1) = max(log((1-pB(echain{i1},1)) ./ pB(echain{i1},1))) + delta;
                
            end
        end
    end
        
    keep = [keep ; keep];   
    
    edgeOn(mod(echain{minind}-1,ne)+1) = false; 
    evalid = [edgeOn ; edgeOn];    
    
    % remove extra edges, regions
    if 1 || r1~=r2
        validr(r2) = false;
    end
    rsp{newr} = [rsp{r1} rsp{r2}];    
    echain = echain(keep);
    spLR = spLR(keep, :);
    Ecost = Ecost(keep);
    
    e2chainshift = max(cumsum(keep),1);
    e2chain = e2chainshift(e2chain);
    
end

%movie2avi(movframe, './tmp/merging.avi')

% if DO_DISPLAY
%     disp(num2str(sum(validr)))
% end

lab = ones(ne, 1);
energy = sum(Eoff(~edgeOn));  %0;
for c = 1:npair
    [minval, minlab] = min([sum(Eon(echain{c})) sum(Eon(echain{c+npair}))]);
    energy = energy + minval;
    lab(mod(echain{c}-1, ne)+1) = minlab+1;
end 

rshift = cumsum(validr);
spLR= rshift(spLR);

% lab(edgeOn) = minlab+1;
% lab(~edgeOn) = 1;
result.edgeletChain = echain;
result.chainLR = spLR;
result.regions = rsp(validr);
result.edgeLabels = lab;
result.approxEnergy = energy;

if DO_VAL
    valdata.mergeCost = mergeCost(1:nval);
    valdata.segError = segError(1:nval);
end

dispim = displayRegions(result.regions, bndinfo.wseg, DO_DISPLAY);




%% Display function
function tmpim = displayRegions(rsp, wseg, fignum)

nsp = max(wseg(:));

elab = zeros(nsp, 1);
for k = 1:numel(rsp)
    elab(rsp{k}) = k;
end  

tmpim = label2rgb(elab(wseg));
[tx, ty] = gradient(elab(wseg));
te = tx~=0 | ty~=0;
tmpim(repmat(te, [1 1 3])) = 0;

if fignum > 0
    figure(fignum), hold off, imagesc(tmpim), axis image, drawnow;
end