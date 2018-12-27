function combinedFeatures = getCombinedFeatures(bndinfo, im)
%% 
imsize = bndinfo.imsize;
eindices = bndinfo.edges.indices;
ejunction = bndinfo.edges.junctions;
jPos = bndinfo.junctions.position;

%% Tjfeatures
Tjinfo =  {};
for k = 1:size(jPos, 1)
    adjacentEdgeIndexs = find(ejunction == k);
    
    if numel(adjacentEdgeIndexs) ~= 3
        % not T-junction
        continue;
    end
    adjacentEdgeInfo = [];
    adjacentEdgeInfo.XY = cell(3,1);
    adjacentEdgeInfo.directionVector = cell(3,1);
    adjacentEdgeInfo.angles = cell(3,1);
    adjacentEdgeInfo.edgeConvexityFeature = cell(3,1);
    for k1 = 1:numel(adjacentEdgeIndexs)
        [i, j] = ind2sub(size(ejunction), adjacentEdgeIndexs(k1));
        [edgeX, edgeY] = ind2sub(imsize(1:2), double(eindices{i}));
        edgeXY = [edgeX, edgeY];
        if j == 2
            edgeXY = flip(edgeXY);
        end

        dv = clacDirectionVector(edgeXY);
        ecf = getEdgeConvexityFeature(edgeXY);
        
        adjacentEdgeInfo.XY{k1, 1}  = edgeXY;
        adjacentEdgeInfo.directionVector{k1, 1} =  dv;
        adjacentEdgeInfo.edgeConvexityFeature{k1, 1} = ecf;
        
%         vX = edgeX;
%         vY = vX * (dv(2)/ dv(1));
%         vY = vY + (edgeY(1) - vY(1));
%         plot(edgeXY(:, 1),edgeXY(:, 2), vX, vY);
%         plot(edgeX, edgeX * (dv(2)/ dv(1)),'-');
    end
    
    dvs = adjacentEdgeInfo.directionVector;
    angle1 = clacVectorsAngle(dvs{1}, dvs{2});
    angle2 = clacVectorsAngle(dvs{1}, dvs{3});
    angle3 = clacVectorsAngle(dvs{2}, dvs{3});
    
    adjacentEdgeInfo.angles{1} =  [angle1, angle2];
    adjacentEdgeInfo.angles{2} =  [angle1, angle3];
    adjacentEdgeInfo.angles{3} =  [angle2, angle3];
    Tjinfo{end + 1, 1} = adjacentEdgeInfo;
end

X = getFeatures(bndinfo, im, bndinfo.pbim, bndinfo.gconf);

combinedFeatures.TJInfo = Tjinfo;
combinedFeatures.TJnum = size(Tjinfo, 1);

combinedFeatures.edgeInfo = X.edge;



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
    
%% calclate angle between two vectors
function angle = clacVectorsAngle(v1, v2)

angle = acos(dot(v1, v2) / norm(v1) * norm(v2));

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
ECFeature =  ECFeature / sum(ECFeature);
