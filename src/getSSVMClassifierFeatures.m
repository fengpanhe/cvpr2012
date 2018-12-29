function [features, lables] = getSSVMClassifierFeatures(bndinfo, combinedFeatures)
% X is the raw data
% ind is the set of indices for which the edgeFeatures should be computed
%
% edgeFeatures:
%   Edge edgeFeatures (1-6)
%        1:  Pb
%        2:  Length / Perimeter
%        3:  Smoothness
%        4:  Angle
%      5-6:  Continuity
%      7-8:  Convexity (area and ratio) - not used
%        9:  Chain length


boundarylabs =  bndinfo.edges.boundaryType;
boundarylabs = (boundarylabs(1:end/2)>0) + 2*(boundarylabs(end/2+1:end)>0);

ind = (1:bndinfo.ne);
edgeFeatures = getEdgeFeatures(bndinfo, combinedFeatures.edgeInfo, ind);

TJInfos =  combinedFeatures.TJInfo;
TJnum = numel(TJInfos);

features = zeros([TJnum * 3, 100], 'single');
lables = zeros([TJnum * 3, 1], 'single');

for k = 1:TJnum
    col = 1;
    row = k * 3 - 2 : k * 3;

    tjinfo = TJInfos{k};

    features(row, col) = k;
    col = col + 1;

    features(row, col : col + 1) = cell2mat(tjinfo.angles);
    col = col + 2;

    features(row, col : col + 35) = cell2mat(tjinfo.edgeConvexityFeature);
    col = col + 36;

%     features(row, col : col + 8) = edgeFeatures(tjinfo.edgeId, 1 : 9);
    col = col + 9;
end



function edgeFeatures = getEdgeFeatures(bndinfo, edgeInfo, ind)
%% Edge Features
ndata = numel(ind);

edgeFeatures = zeros([ndata 9], 'single');

if isempty(ind)
    return;
end

spLR = bndinfo.edges.spLR;
s1 = spLR(ind, 1);
s2 = spLR(ind, 2);

f = 0;


edgeFeatures(:, f+1) = edgeInfo.pb(ind);

perim = zeros(bndinfo.nseg, 1);
for k = 1:numel(edgeInfo.length)
    perim(spLR(k, 1)) = perim(spLR(k, 1)) + edgeInfo.length(k);
    perim(spLR(k, 2)) = perim(spLR(k, 2)) + edgeInfo.length(k);
end
minperim = min([perim(s1) perim(s2)], [], 2);
edgeFeatures(:, f+2) = edgeInfo.length(ind) ./ minperim; % edge length / perim of smaller region

edgeFeatures(:, f+3) = edgeInfo.smoothness(ind); % measure of smoothess
    
theta = edgeInfo.theta;
% discrete angle
edgeFeatures(:, f+4) = max(ceil((mod(theta(ind),2*pi) / (pi*2) * 16 - 1E-10)),1); 

% relative angle (continuity)
theta1 = mod(edgeInfo.thetaStart*180/pi, 360);
theta2 = mod(edgeInfo.thetaEnd*180/pi, 360);
maxc = zeros(ndata, 2);
eadj = bndinfo.edges.adjacency;
ne = bndinfo.ne;
for k = 1:ndata
    ki = ind(k);
    ra = abs(theta2(ki)-theta1(eadj{ki}));
    ra = ra - 180*(ra>180);
    if isempty(ra), maxc(k,1) = 0;
    else maxc(k,1) = min(ra);
    end    
    ra = mod(abs(theta2(ne+ki)-theta1(eadj{ne+ki})), 180+1E-5);         
    if isempty(ra), maxc(k,2) = 0;
    else maxc(k,2) = min(ra);
    end    
end
edgeFeatures(:, f+(5:6)) = [min(maxc, [], 2) max(maxc, [], 2)];

edgeFeatures(:, f+8) = edgeInfo.convRatio;

edgeFeatures(:, f+9) = edgeInfo.edge2chain(ind);

f = f + 9;