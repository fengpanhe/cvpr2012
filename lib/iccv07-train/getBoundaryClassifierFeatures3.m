function [tx, categoryFeatures] = getBoundaryClassifierFeatures3(bndinfo, X, ind)
% X is the raw data
% ind is the set of indices for which the features should be computed
%
% tx:
%   Edge features (1-6)
%        1:  Pb
%        2:  Length / Perimeter
%        3:  Smoothness
%        4:  Angle
%      5-6:  Continuity
%   Region features (7-17)
%      7-8:  Area
%        9:  Color Mean Difference
%       10:  Color Entropy Difference
%       11:  Gradient Entropy Difference
%    12-13:  Position (x,y)
%    14-15:  Extent Overlap (x, y)
%   Geometry features (16-39)
%    16-25:  Geometric Context Mean
%    26-30:  Geometric Context Difference
%    31   :  Geometric Context Sum Abs Difference
%    32-33:  Geometric Context Most Likely Label (G V or S)
%    34-35:  Geometric T-junctions




ng = 5; % five geometric classes

ndata = numel(ind);

[imh, imw] = size(bndinfo.wseg);

tx = zeros([ndata 35], 'single');

spLR = bndinfo.edges.spLR;
s1 = spLR(ind, 1);
s2 = spLR(ind, 2);

categoryFeatures = [];
f = 0;


%% Edge features
tx(:, f+1) = X.edge.pb(ind);

perim = zeros(bndinfo.nseg, 1);
for k = 1:numel(X.edge.length)
    perim(spLR(k, 1)) = perim(spLR(k, 1)) + X.edge.length(k);
    perim(spLR(k, 2)) = perim(spLR(k, 2)) + X.edge.length(k);
end
minperim = min([perim(s1) perim(s2)], [], 2);
tx(:, f+2) = X.edge.length(ind) ./ minperim; % edge length / perim of smaller region

juncts = bndinfo.edges.junctions(ind, :);
jpos1 = bndinfo.junctions.position(juncts(:, 1), :);
jpos2 = bndinfo.junctions.position(juncts(:, 2), :);
directLength = abs(jpos2(:, 1)-jpos1(:,1)) + abs(jpos2(:, 2)-jpos1(:,2));
tx(:, f+3) = directLength ./ X.edge.length(ind); % measure of smoothess
    
theta = X.edge.theta;
% discrete angle
tx(:, f+4) = max(ceil((mod(theta(ind),2*pi) / (pi*2) * 16 - 1E-10)),1); 
categoryFeatures(end+1) = f+4;

% relative angle (continuity)
theta = mod([theta ; theta+pi]/pi*180, 360);
maxc = zeros(ndata, 2);
eadj = bndinfo.edges.adjacency;
ne = bndinfo.ne;
for k = 1:ndata
    ki = ind(k);
    ra = mod(abs(theta(ki)-theta(eadj{ki})), 180+1E-5);
    if isempty(ra), maxc(k,1) = 0;
    else maxc(k,1) = min(ra);
    end    
    ra = mod(abs(theta(ne+ki)-theta(eadj{ne+ki})), 180+1E-5);         
    if isempty(ra), maxc(k,2) = 0;
    else maxc(k,2) = min(ra);
    end    
end
tx(:, f+(5:6)) = [min(maxc, [], 2) max(maxc, [], 2)];

f = f + 6;


%% Region features

% area
area1 = X.region.area(s1);
area2 = X.region.area(s2);
tx(:, f+(1:2)) = [min([area1 area2], [], 2) max([area1 area2], [], 2)];

% color
tx(:, f+3) = sqrt(sum((X.region.colorMean(s1, :)-X.region.colorMean(s2, :)).^2, 2));

ch = X.region.colorHist+1E-10;
for k = 1:ndata
    h1 = ch(s1(k), :);  e1 = sum(-log(h1).*h1);
    h2 = ch(s2(k), :);  e2 = sum(-log(h2).*h2);
    e12 = (area1(k)*e1 + area2(k)*e2)/(area1(k)+area2(k));
    h3 = (area1(k)*h1 + area2(k)*h2)/(area1(k)+area2(k));
    e3 = sum(-log(h3).*h3);
    tx(k, f+4) = e3-e12;
end

% gradient
% ch = X.region.gradHist+1E-10;
% for k = 1:ndata
%     h1 = ch(s1(k), :);  e1 = sum(-log(h1).*h1);
%     h2 = ch(s2(k), :);  e2 = sum(-log(h2).*h2);
%     e12 = (area1(k)*e1 + area2(k)*e2)/(area1(k)+area2(k));
%     h3 = (area1(k)*h1 + area2(k)*h2)/(area1(k)+area2(k));
%     e3 = sum(-log(h3).*h3);
%     tx(k, f+5) = e3-e12;
% end

% position
tx(:, f+6) = (X.region.x(s1, 2)-X.region.x(s2, 2)) / imw;
tx(:, f+7) = (X.region.y(s1, 2)-X.region.y(s2, 2)) / imh;
% x overlap
x1 = X.region.x(s1, [1 3]);
x2 = X.region.x(s2, [1 3]);
tx(:, f+8) = (min([x1(:, 2) x2(:, 2)], [], 2)-max([x1(:, 1) x2(:, 1)], [], 2)) ./ ...
    (max([x1(:, 2) x2(:, 2)], [], 2)-min([x1(:, 1) x2(:, 1)], [], 2));
% y overlap
y1 = X.region.y(s1, [1 3]);
y2 = X.region.y(s2, [1 3]);
tx(:, f+9) = (min([y1(:, 2) y2(:, 2)], [], 2)-max([y1(:, 1) y2(:, 1)], [], 2)) ./ ...
    (max([y1(:, 2) y2(:, 2)], [], 2)-min([y1(:, 1) y2(:, 1)], [], 2));

f = f + 9;


%% 3D Geometry features

% geometric context features
gc = X.region.geomContext;

tx(:, f+(1:ng)) = gc(s1, :);
tx(:, f+ng+(1:ng)) = gc(s2, :);
tx(:, f+2*ng+(1:ng)) = tx(:, f+(1:ng))-tx(:, f+ng+(1:ng));
tx(:, f+3*ng+1) = sum(abs(tx(:, f+2*ng+(1:ng))), 2)/2;

[maxval, maxlab] = max([gc(:, 1) sum(gc(:, 2:4), 2) gc(:, 5)], [], 2);
tx(:, f+3*ng+(2:3)) = [maxlab(s1) maxlab(s2)];
categoryFeatures(end+(1:2)) = f+3*ng+(2:3);

f = f + 18;

% geometric T-junctions
% 1) j = lower junction
% 2) s3 is third region that j touches
% 3) if (s3_c below j) and (s2_c below s1_c) and (s2_c below j)
%       and (s1=vert, s2=vert, s3=gnd) then s2 is probably in front of s1
% 4) if (s3_c below j) and (s1_c below s2_c) and (s1_c below j)
%       and (s1=vert, s2=vert, s3=gnd) then s1 is probably in front of s2
spLR = bndinfo.edges.spLR;
spLR = [spLR ; spLR(:, [2 1])];
for k = 1:ndata
    [jy, jid] = max([jpos1(k, 2) jpos2(k, 2)]); % get lower of junctions        
    tk = ind(k) + ne*(jid==2);
    if ~isempty(eadj{tk})
        nexte = eadj{tk}(1);      
        ni = (spLR(nexte, :)~=s1(k)) & (spLR(nexte, :)~=s2(k));
        s3 = spLR(nexte, ni);   
        %if isempty(s3), wasempty = wasempty+1; end       
        if numel(s3)>1, s3 = s3(1); end
        if X.region.y(s3, 2) > jy  % s3 below junction                       

           if X.region.y(s2(k), 2)>X.region.y(s1(k), 2) && ...
                   X.region.y(s2(k), 3) > jy+5
               % s2 below s1 and below junction
               tx(k, f+1) = gc(s3, 1)*sum(gc(s1(k), 2:4))*sum(gc(s2(k), 2:4));
           elseif X.region.y(s1(k), 2)>X.region.y(s2(k), 2) && ...
                   X.region.y(s1(k), 3) > jy+5
                % s1 below s2 and below junction
               tx(k, f+1) = -gc(s3, 1)*sum(gc(s1(k), 2:4))*sum(gc(s2(k), 2:4)); 
           end            
        end
    end
end
tx(:, f+2) = abs(tx(:, f+1));




        
    
