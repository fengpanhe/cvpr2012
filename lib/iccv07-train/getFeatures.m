function [X, Y] = getFeatures(bndinfo, im, pbim, gconf)
% [X, Y] = getBoundaryFeatures(bndinfo, pbim, gconf)
% 
% Computes features for each edgelet based on boundary and geometry
% confidence images.  This version (3) computes geometry confidences as
% means of superpixel values on either side of the edglet.
%
% Input:
%   bndinfo - structure of superpixel boundaries
%   in(imh, imw, 3) - color image
%   pbim(imh, imw, norient) - probability of boundary image at each orientation
%   gconf(imh, imw, ngeoms) - confidence in each geometric label
%
% Output:
%   X(nboundaries, :) - features for the boundaries
%   Y(nboundaries, 2) - geometry on each size or 0 for no edge


%% Initialization

nsp = bndinfo.nseg;
edges = bndinfo.edges;
nbnd = bndinfo.ne*2; % number of directed edges
ne = bndinfo.ne;
[imh, imw] = size(bndinfo.wseg);

% set edge labels if ground truth is available
Y = zeros(nbnd, 1);
if isfield(bndinfo.edges, 'boundaryType')
    Y = bndinfo.edges.boundaryType;
end

X.edge.pb = zeros(ne, 1);
X.edge.theta = zeros(ne, 1);
X.edge.length = zeros(ne, 1);

X.region.colorMean = zeros(nsp, 3); % Lab color mean
X.region.colorHist = [];  % set number of bins later (max 512)
%X.region.gradHist = []; % set number of bins later (max 64)
X.region.x = zeros(nsp, 3); % min, mean, max
X.region.y = zeros(nsp, 3); % min, mean, max
X.region.area = zeros(nsp, 1);
X.region.geomContext = zeros(nsp, 5);


%% Edge statistics

% get discretize orientation into 1=1|2, 2=1/2, 3=1_2, 4=2\1
theta = bndinfo.edges.thetaUndirected;
rels = (theta < -3*pi/8) + (theta < -pi/8) + (theta < pi/8) +  (theta < 3*pi/8); 
rels = mod(rels, 4) + 1;

X.edge.theta = bndinfo.edges.thetaDirected;  % directed edge angle

% pbim(imh, imw, [ -, \, |, / ]) (direction of edge)
pbmap = [3 4 1 2];

%X.edge.pbOrient = zeros(ne,4);
% compute features
for k = 1:ne 

    eind = edges.indices{k};
    
    X.edge.length(k) = numel(eind); % edge length
            
    pbsubi = pbmap(rels(k));    
    ind = eind + (pbsubi-1)*imh*imw;
    X.edge.pb(k) = sum(pbim(ind))/numel(ind);  % mean pb               
    
    %X.edge.pbOrient(k, pbsubi) = X.edge.pb(k);
    
end


%% Region statistics

% get area and position stats
stats = regionprops(bndinfo.wseg, 'PixelIdx', 'Area', 'Centroid', 'BoundingBox');

area = cat(1, stats(:).Area);
X.region.area = area / (imh*imw);

bbox = cat(1, stats(:).BoundingBox);
centroid = cat(1, stats(:).Centroid);
minx = bbox(:,1);
meanx = centroid(:, 1);
maxx = minx + bbox(:,3);
miny = bbox(:,2);
meany = centroid(:, 2);
maxy = miny + bbox(:,4);
X.region.x = [minx meanx maxx]; % min, mean, max
X.region.y = [miny meany maxy]; % min, mean, max

% get lab color image
im = RGB2Lab(im);

% get discrete image with nb values per channel
nb = 8;
mincols = repmat(min(min(im)), [imh imw]);
maxcols = repmat(max(max(im))+1E-10, [imh imw]);
imd = floor((im-mincols)./(maxcols-mincols)*nb);
imd = imd(:, :, 1) + imd(:, :, 2)*nb + imd(:,:,3)*nb*nb + 1;

% get discrete texture image (does not help much)
% [gx, gy] = gradient(im(:, :, 1));
% gx = log(abs(gx)+1);
% gy = log(abs(gy)+1);
% gradim = [log(abs(gx(:))+1)  log(abs(gy(:))+1)];
% mingrad = repmat(min(gradim), [imh*imw 1]);
% maxgrad = repmat(max(gradim)+1E-10, [imh*imw 1]);
% gradim = floor((gradim-mingrad)./(maxgrad-mingrad)*nb);
% gradim = gradim(:, 1) + nb*gradim(:, 2) + 1;

% make gconf(:, :, [gnd planar porous solid sky])
gconf = cat(3, gconf(:, :, 1), sum(gconf(:, :, 2:4), 3), gconf(:, :, 5:7));

% compute mean color
idx = {stats(:).PixelIdxList};
im = reshape(im, [imh*imw 3]);
for c = 1:3
    cim = im(:, c);
    for k = 1:nsp    
        X.region.colorMean(k,c) = sum(cim(idx{k})) / area(k);
    end
end

% compute mean geometric context
gconf = reshape(gconf, [imh*imw 5]);
for g = 1:5
    gim = gconf(:, g);
    for k = 1:nsp    
        X.region.geomContext(k,g) = sum(gim(idx{k})) / area(k);
    end
end

% compute histograms of color and texture
imd =reshape(imd, [imh*imw 1]);
wseg = bndinfo.wseg;
colorHist = zeros(nsp, nb*nb*nb);
%gradHist = zeros(nsp, nb*nb);
for k = 1:imh*imw
    s = wseg(k);
    colorHist(s, imd(k)) = colorHist(s, imd(k)) + 1;
    %gradHist(s, gradim(k)) = gradHist(s, gradim(k)) + 1;
end

keep = sum(colorHist, 1)>0;  % only keep bins that have at least one value
colorHist = colorHist(:, keep);
X.region.colorHist = single(colorHist ./ repmat(area, [1 sum(keep)]));

% keep = sum(gradHist, 1)>0;  % only keep bins that have at least one value
% gradHist = gradHist(:, keep);
% X.region.gradHist = single(gradHist ./ repmat(area, [1 sum(keep)])); 
