function combinedFeatures = getCombinedFeatures(bndinfo, im)
%% 
imsize = bndinfo.imsize;
eindices = bndinfo.edges.indices;
ejunction = bndinfo.edges.junctions;
jPos = bndinfo.junctions.position;

%% Tjfeatures
% Tjfeatures = 
for k = 1:size(jPos)
    adjacentEdgeIndexs = find(ejunction == k);
    
    if numel(adjacentEdgeIndexs) ~= 3
        % not T-junction
        continue;
    end
    adjacentEdgeInfo = cell(3,1);
    for k1 = 1:numel(adjacentEdgeIndexs)
        [i, j] = ind2sub(size(ejunction), adjacentEdgeIndexs(k1));
        [edgeX, edgeY] = ind2sub(imsize(1:2), double(eindices{i}));
        edgeXY = [edgeX, edgeY];
        if j == 2
            edgeXY = flip(edgeXY);
        end

        dv = clacDirectionVector(edgeXY);
        ecf = getEdgeConvexityFeature(edgeXY);
        
        adjacentEdgeInfo{k1,1}.XY  = edgeXY;
        adjacentEdgeInfo{k1,1}.directionVector =  dv;
        adjacentEdgeInfo{k1,1}.edgeConvexityFeature = ecf;
        
%         vX = edgeX;
%         vY = vX * (dv(2)/ dv(1));
%         vY = vY + (edgeY(1) - vY(1));
%         plot(edgeXY(:, 1),edgeXY(:, 2), vX, vY);
%         plot(edgeX, edgeX * (dv(2)/ dv(1)),'-');
    end
end


%% Calculate the direction vector of a set of points 
function directionVector = clacDirectionVector(pos)
X = pos(:,1);
Y = pos(:,2);
p = polyfit(X, Y, 1);
vector = [1,p(1)];
if mean(X) < X(1)
    vector = -vector;
end
directionVector = vector / norm(vector);
    
%% get the convexity feature of edge
function ECFeature = getEdgeConvexityFeature(edgeXY)

lb = edgeXY(end, :) - edgeXY(1, :);
lb_atan2d = atan2d(lb(2), lb(1));

thetas = zeros(1, size(edgeXY, 1));
for i = 2:size(edgeXY, 1)
    li = edgeXY(i, :) - edgeXY(1, :);
    li_atan2d = atan2d(li(2), li(1));
    theta_i2b = li_atan2d - lb_atan2d;
    if theta_i2b > 180
        theta_i2b = theta_i2b - 360;
    end
    if theta_i2b < -180
        theta_i2b = theta_i2b + 360;
    end
    thetas(i) =  theta_i2b;
end

edges = [-180  -170  -160  -150  -140  -130  -120  -110  -100   -90   -80   -70   -60   -50   -40   -30   -20   -10     0    10    20    30    40    50    60    70    80    90   100   110   120   130 140   150   160   170   180];

ECFeature = histcounts(thetas, edges);
