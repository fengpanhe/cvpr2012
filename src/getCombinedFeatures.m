function combinedFeatures = getCombinedFeatures(bndinfo, im)

imsize = bndinfo.imsize;
eindices = bndinfo.edges.indices;
ejunction = bndinfo.edges.junctions;
jPos = bndinfo.junctions.position;

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
        edgeXY = ind2sub(imsize(1:2), double(eindices{i}));
        if j == 2
            edgeXY = flip(edgeXY);
        end
        adjacentEdgeInfo{k1,1}.XY  = edgeXY;
        adjacentEdgeInfo{k1,1}.directionVector =  clacDirectionVector(edgeXY);
    end
end


function directionVector = clacDirectionVector(pos)
    X = pos(:,1);
    Y = pos(:,2);
    p = polyfit(X, Y, 1);
    vector = [1,p(1)];
    if mean(X) < X(1)
        vector = -vector;
    end
    directionVector = vector / norm(vector);
