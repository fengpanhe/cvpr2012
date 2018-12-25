function [X, Y] = getBoundaryFeatures3(bndinfo, pbim, gconf)
% [X, Y] = getBoundaryFeatures(bndinfo, pbim, gconf)
% 
% Computes features for each edgelet based on boundary and geometry
% confidence images.  This version (3) computes geometry confidences as
% means of superpixel values on either side of the edglet.
%
% Input:
%   bndinfo - structure of superpixel boundaries
%   pbim(imh, imw, norient) - probability of boundary image at each orientation
%   gconf(imh, imw, ngeoms) - confidence in each geometric label
%
% Output:
%   X(nboundaries, :) - features for the boundaries
%   Y(nboundaries, 2) - geometry on each size or 0 for no edge


edges = bndinfo.edges;
nbnd = bndinfo.ne*2; % number of directed edges
ne = bndinfo.ne;

X = zeros(ne, 19);
Y = zeros(nbnd, 1);

if isfield(bndinfo.edges, 'boundaryType')
    Y = bndinfo.edges.boundaryType;
end

% make gconf(:, :, [gnd planar porous solid sky])
gconf = cat(3, gconf(:, :, 1), sum(gconf(:, :, 2:4), 3), gconf(:, :, 5:7));

[imh, imw, ng] = size(gconf);

% get mean confidence in each superpixel
gcabs = zeros(bndinfo.nseg, ng);
stats = regionprops(bndinfo.wseg, 'PixelIdxList');
idx = {stats(:).PixelIdxList};
for g = 1:ng
    gim = reshape(gconf(:, :, g), [imh*imw 1]);    
    for k = 1:numel(idx)
        gcabs(k, g) = sum(gim(idx{k}))/numel(idx{k});
    end
end

% get discretize orientation into 1=1|2, 2=1/2, 3=1_2, 4=2\1
theta = bndinfo.edges.thetaUndirected;
rels = (theta < -3*pi/8) + (theta < -pi/8) + (theta < pi/8) +  (theta < 3*pi/8); 
rels = mod(rels, 4) + 1;

% need to switch some values if thetaDirected~=thetaUndirected
theta2 = bndinfo.edges.thetaDirected;
doSwitch = (theta2 <= -pi/2) | (theta2 > pi/2);
doSwitch(ne+1:nbnd) = ~doSwitch;

thetaDirected = repmat(theta, [2 1]);
thetaDirected(doSwitch) = thetaDirected(doSwitch) + pi;
nAngleBins = 16;
thetaDirected = mod(thetaDirected, 2*pi);
thetaDirected = max(ceil(thetaDirected / (2*pi) * nAngleBins), 1);

% pbim(imh, imw, [ -, \, |, / ]) (direction of edge)
% rels(k): 1|2, 1/2, 1_2, 2\1 (relationship)

pbmap = [3 4 1 2];
spLR = bndinfo.edges.spLR;
spLR = [spLR ; spLR(:, [2 1])];

juncts = bndinfo.edges.junctions;
jpos = bndinfo.junctions.position;
jpos = jpos ./ repmat([imh imw], [size(jpos, 1) 1]);
jpos1 = jpos(juncts(:, 1), :);
jpos2 = jpos(juncts(:, 2), :);
jposMean = (jpos1 + jpos2)./2;

% compute features
for k = 1:ne 
    rel = rels(k);
    eind = edges.indices{k};
    pbsubi = pbmap(rel);
    
    s1 = spLR(k, 1);
    s2 = spLR(k, 2);
    
    if ~isempty(eind)

        X(k, (1:ng)) = gcabs(s1, :); % f1-f5 (geom mean left)
        X(k, 2*ng+(1:ng)) = gcabs(s2, :); % f11-f15  (geom mean right)
        X(k, ng+(1:ng)) = X(k, (1:ng))-X(k, 2*ng+(1:ng)); % f6-f10 (geom diff left-right)           
        
        ind = eind + (pbsubi-1)*imh*imw;
        X(k, 3*ng+1) = sum(pbim(ind))/numel(ind); % f16 pb
        X(k, 3*ng+2) = thetaDirected(k); % f17 directed angle          
        
        % edge position
        X(k, 3*ng+3) = jposMean(k, 1);  % f18 x
        X(k, 3*ng+4) = jposMean(k, 2);  % f19 y           
    end
end
   
