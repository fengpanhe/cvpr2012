function [X2, Y2] = updateFeatures(bndinfo, X)
% [X2, Y2] = updateFeatures(bndinfo, X, result)
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

ne = bndinfo.ne;
nr = bndinfo.nseg;

X2.edge.pb = zeros(ne, 1);
X2.edge.theta = zeros(ne, 1);
X2.edge.length = zeros(ne, 1);

X2.region.colorMean = zeros(nr, 3); % Lab color mean
X2.region.colorHist = zeros([nr size(X.region.colorHist,2)], 'single'); 
X2.region.gradHist = zeros([nr size(X.region.gradHist,2)], 'single');
X2.region.x = zeros(nr, 3); % min, mean, max
X2.region.y = zeros(nr, 3); % min, mean, max
X2.region.area = zeros(nr, 1);
X2.region.geomContext = zeros(nr, 5);


%% Assign labels (if required)

if nargout==2
    Y2 = bndinfo.edges.boundaryType;
end


%% Edge statistics

chains = bndinfo.edges.chains;
elength = X.edge.length;

X2.edge.theta = bndinfo.edges.thetaDirected;  % directed edge angle

% compute features
for k = 1:ne 
    ind = mod(chains{k}-1, bndinfo.origne)+1;
    X2.edge.length(k) = sum(elength(ind));
    X2.edge.pb(k) = sum(X.edge.pb(ind).*elength(ind)) / X2.edge.length(k);                      
end


%% Region statistics

area = bndinfo.spArea(:);
regions = bndinfo.regions2sp;

for k = 1:numel(regions)
    
    ind = regions{k};    
    spArea = area(ind);
    totalArea = sum(spArea);
    
    % get area and position stats
    X2.region.area(k) = totalArea;
    
    minx = min(X.region.x(ind, 1));
    maxx = max(X.region.x(ind, 3));
    meanx = sum(spArea .* X.region.x(ind, 2)) / totalArea;
    X2.region.x(k, :) = [minx meanx maxx];
    
    miny = min(X.region.y(ind, 1));
    maxy = max(X.region.y(ind, 3));
    meany = sum(spArea .* X.region.y(ind, 2)) / totalArea;
    X2.region.y(k, :) = [miny meany maxy];    


    X2.region.colorMean(k, :) = sum(repmat(spArea, [1 3]) .* ...
        X.region.colorMean(ind, :)) / totalArea;
    
    X2.region.geomContext(k,:) = sum(repmat(spArea, [1 5]) .* ...
        X.region.geomContext(ind, :)) / totalArea;
        
    sh = size(X.region.colorHist, 2);
    X2.region.colorHist(k, :) = sum(repmat(spArea, [1 sh]) .* ...
        X.region.colorHist(ind, :)) / totalArea;

    sh = size(X.region.gradHist, 2);
    X2.region.gradHist(k, :) = sum(repmat(spArea, [1 sh]) .* ...
        X.region.gradHist(ind, :)) / totalArea;

end    
