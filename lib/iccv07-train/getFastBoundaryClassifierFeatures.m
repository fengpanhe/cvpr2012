function [tx, categoryFeatures] = getFastBoundaryClassifierFeatures(bndinfo, X, ind)
% X is the raw data
% ind is the set of indices for which the features should be computed
%
% tx:
%   Edge features (1-6)
%        1:  Pb
%        4:  Angle
%      5-6:  Continuity
%      7-8:  Convexity (area and ratio) - not used
%        9:  Chain length
%   Region features (7-17)+3
%    10-11:  Area
%       12:  Color Mean Difference
%       13:  Color Entropy Difference
%       14:  Gradient Entropy Difference - not used
%    15-16:  Position (x,y)
%    17-18:  Extent Overlap (x, y)
%    an additional 10 features of position/overlap
%   Geometry features (16-39)+3
%    19-28:  Geometric Context Mean
%    29-33:  Geometric Context Difference
%    34   :  Geometric Context Sum Abs Difference
%    35-36:  Geometric Context Most Likely Label (G V or S)
%    37-40:  Depth under- and over-estimates for each side
%    31-43:  Depth, min1-min2, max1-max2, min(max12) - max(min12)
%    44-47:  Depthcol, each sp, diff, abs diff




ng = 5; % five geometric classes

ndata = numel(ind);

[imh, imw] = size(bndinfo.wseg);

tx = zeros([ndata 33], 'single');

spLR = bndinfo.edges.spLR;
s1 = spLR(ind, 1);
s2 = spLR(ind, 2);

categoryFeatures = [];
f = 0;


%% Edge features
tx(:, f+1) = X.edge.pb(ind);
    
area1 = X.region.area(s1);
area2 = X.region.area(s2);

f = f+1;

%% Region features

% area

tx(:, f+(1:2)) = [min([area1 area2], [], 2) max([area1 area2], [], 2)];

% color
tx(:, f+3) = sqrt(sum((X.region.colorMean(s1, :)-X.region.colorMean(s2, :)).^2, 2));

% position
tx(:, f+4) = (X.region.y(s1, 3))-(X.region.y(s2, 3)); % difference of tops
tx(:, f+5) = (X.region.y(s1, 1))-(X.region.y(s2, 1)); % difference of bottoms
tx(:, f+6) = (X.region.y(s1, 3))-(X.region.y(s2, 1)); % top1 - bottom2
tx(:, f+7) = (X.region.y(s1, 1))-(X.region.y(s2, 3)); % bottom1 - top2
tx(:, f+8) = X.region.y(s1, 3) - X.region.y(s1, 1); % top1 - bottom1
tx(:, f+9) = X.region.y(s2, 3) - X.region.y(s2, 1); % top2 - bottom2
tx(:, f+10) = (X.region.x(s1, 1))-(X.region.x(s2, 1)); % left1 - left2
tx(:, f+11) = (X.region.x(s1, 3))-(X.region.x(s2, 3)); % right1 - right2
tx(:, f+12) = X.region.x(s1, 3) - X.region.x(s1, 1); % right1 - left1
tx(:, f+13) = X.region.x(s2, 3) - X.region.x(s2, 1); % right2 - left2

% x alignment
x1 = X.region.x(s1, [1 3]);
x2 = X.region.x(s2, [1 3]);
tx(:, f+14) = (min([x1(:, 2) x2(:, 2)], [], 2)-max([x1(:, 1) x2(:, 1)], [], 2)) ./ ...
    (max([x1(:, 2) x2(:, 2)], [], 2)-min([x1(:, 1) x2(:, 1)], [], 2));

% y overlap
y1 = X.region.y(s1, [1 3]);
y2 = X.region.y(s2, [1 3]);
tx(:, f+15) = (min([y1(:, 2) y2(:, 2)], [], 2)-max([y1(:, 1) y2(:, 1)], [], 2)) ./ ...
    (max([y1(:, 2) y2(:, 2)], [], 2)-min([y1(:, 1) y2(:, 1)], [], 2));

f = f + 15;


%% 3D Geometry features

% geometric context features
gc = X.region.geomContext;

tx(:, f+(1:ng)) = gc(s1, :);
tx(:, f+ng+(1:ng)) = gc(s2, :);
tx(:, f+2*ng+(1:ng)) = tx(:, f+(1:ng))-tx(:, f+ng+(1:ng));
tx(:, f+3*ng+1) = sum(abs(tx(:, f+2*ng+(1:ng))), 2)/2;

[maxval, maxlab] = max([gc(:, 1) sum(gc(:, 2:4), 2) gc(:, 5)], [], 2);
tx(:, f+3*ng+2) = (maxlab(s1)-1)*3+ maxlab(s2);
categoryFeatures(end+1) = f+3*ng+2;

f = f + 17;
