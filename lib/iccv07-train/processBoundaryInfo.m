function bndinfo = processBoundaryInfo(bndinfo)

wseg = bndinfo.wseg;

[imh, imw] = size(wseg);

nseg = double(max(wseg(:)));

%% Get adjacent superpixels
edgeim = false(size(wseg));
edgeim(2:end-1, 2:end-1) = (wseg(2:end-1,2:end-1)==0);
nedge = sum(sum(edgeim));

eind = find(edgeim);

ecount = zeros(nseg+1, nseg+1, 'uint16');
[ey, ex] = ind2sub([imh imw], eind); 
  

% for each edge pixel, get which superpixels are adjacent
spcomb = zeros(nedge, 4, 'uint32');
cvalid = false(nedge, 4);

% diag 1, diag 2, vert, horz
[spcomb(:, 1), cvalid(:, 1)] = ...
    getSpVals(wseg, edgeim, [-1 1], [-1 1], nedge, nseg);
[spcomb(:, 2), cvalid(:, 2)] = ...
    getSpVals(wseg, edgeim, [-1 1], [1 -1], nedge, nseg);
[spcomb(:, 3), cvalid(:, 3)] = ...
    getSpVals(wseg, edgeim, [0 0], [-1 1], nedge, nseg);
[spcomb(:, 4), cvalid(:, 4)] = ...
    getSpVals(wseg, edgeim, [-1 1], [0 0], nedge, nseg);

% get unique superpixel pairs
ucomb = sort(spcomb(cvalid(:)));
ucomb = ucomb([true(1) ; ucomb(1:end-1)~=ucomb(2:end)]);
nadj = numel(ucomb);

% create list of adjacent pairs of superpixels
adjlist = zeros(nadj, 2);
[adjlist(:, 1), adjlist(:, 2)] = ind2sub([nseg nseg], ucomb);

% create lookup table from pair of superpixels to index in adjlist
adjLookup = spalloc(nseg, nseg, nadj);
adjLookup(ucomb(:)) = (1:numel(ucomb))';

%% Get edglets between pairs of superpixels

% create cell array of edgelets
edges = cell(nadj, 1);
ecount = zeros(nadj, 1);

eim = zeros(size(wseg), 'uint32');

for k = 1:size(spcomb(:, 1))
    eadj = spcomb(k, cvalid(k, :));    
%     if numel(eadj)>1
%         eadj = sort(eadj);
%         eadj2 = eadj(2:end);
%         eadj = [eadj(1) eadj2(eadj(2:end)~=eadj(1:end-1))];
%     end
    if numel(eadj)>=1 
        eadj = eadj(1);
        eadj = full(adjLookup(eadj));
        ecount(eadj) = ecount(eadj) + 1;
    end
end
for k = 1:numel(ecount)
    edges{k} = zeros([ecount(k) 1], 'uint32');
end
ecount(:) = 0;
for k = 1:size(spcomb(:, 1))
    eadj = spcomb(k, cvalid(k, :));
%     if numel(eadj)>1
%         eadj = sort(eadj);
%         eadj2 = eadj(2:end);
%         eadj = [eadj(1) eadj2(eadj(2:end)~=eadj(1:end-1))];
%     end
    if numel(eadj)>=1
        eadj = eadj(1);
        eadj = full(adjLookup(eadj));
        ecount(eadj) = ecount(eadj) + 1;
        edges{eadj}(ecount(eadj)) = eind(k);
        eim(eind(k)) = eadj;
    end
%     for k2 = 1:numel(eadj)
%         edges{eadj(k2)}(ecount(eadj(k2))) = eind(k);
%     end
end

edgeAdjacency = spalloc(nadj, nadj, round(nadj*2.5));

ejunctions = cell(nadj, 1);

%% Get junctions and adjacency of edglets
junctions = edgeim2 & imfilter(single(edgeim2), [0 1 0; 1 0 1; 0 1 0])>2;
[yj, xj] = find(junctions);
jind = yj + (xj-1)*imh;
nj = numel(yj);
jim = zeros(imh, imw);
jim(yj + (xj-1)*imh) = (1:nj);
jlist = cell(nj, 1);
for k = 1:nj    
    tmpe = reshape(eim(yj(k)-1:yj(k)+1, xj(k)-1:xj(k)+1), [9 1]);
    tmpe = sort(tmpe(tmpe>0));
    if numel(tmpe)>1
        jlist{k} = tmpe([true(1) ; (tmpe(2:end)~=tmpe(1:end-1))]);          
    else
        jlist{k} = tmpe;
    end
end
for k = 1:nj
    tmpj = [jim([yj(k)-1 yj(k)+1], xj(k)) jim(yj(k), [xj(k)-1 xj(k)+1])'];
    tmpj = tmpj(tmpj>0);
    for k2 = 1:numel(tmpj)
        jlist{k} = unique([jlist{k} ; jlist{tmpj(k2)}]);
    end    
    for k1 = 1:numel(jlist{k})
        ejunctions{jlist{k}(k1)}(end+1) = jind(k); 
        for k2 = k1+1:numel(jlist{k})
            edgeAdjacency(jlist{k}(k1), jlist{k}(k2)) = true;            
        end
    end
end
    
bndinfo.nseg = nseg;
bndinfo.edges = edges;
bndinfo.adjlist = adjlist;
bndinfo.adjLookup = adjLookup;

[e1, e2] = find(edgeAdjacency);
edgeAdjList = [e1(:) e2(:)];
bndinfo.edgeAdjacency = edgeAdjacency;
bndinfo.edgeAdjList = edgeAdjList;
bndinfo.ejunctions = ejunctions;



%% Helper function to get indices of adjacent pairs of superpixels
function [spcomb, cvalid] = getSpVals(wseg, edgeim, sx, sy, nedge, nseg)

% get superpixel indices for each (sx,sy)
spvals = zeros([nedge 2], 'uint32');
for k = 1:2
    tmpedgeim = false(size(wseg));
    tmpedgeim(sy(k)+2:end-1+sy(k), sx(k)+2:end-1+sx(k)) = edgeim;
    spvals(:,k) = wseg(tmpedgeim);
end

% ensure that spvals(:,1)<=spvals(:, 2)
ind = spvals(:, 1)>spvals(:, 2);
spvals(ind, :) = [spvals(ind, 2) spvals(ind, 1)];

% check that pair are non-zero and non-equal
cvalid = (spvals(:, 1)>0) & (spvals(:, 2)>0) & (spvals(:, 1)~=spvals(:, 2));

% encode spvals pair as single value
spcomb = zeros(nedge, 1, 'uint32');
spcomb(cvalid) = spvals(cvalid, 1) + (spvals(cvalid, 2)-1)*uint32(nseg);



