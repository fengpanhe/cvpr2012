function bndinfo = labelBoundaries(bndinfo)

for f = 1:numel(bndinfo)
     
    bndinfo(f).edges.spLR;
    
    if isfield(bndinfo(f), 'labels')
        labels = bndinfo(f).labels;
        types = bndinfo(f).type;         
        splabs = zeros(size(spLR));        
        splabs(spLR>0) = labels(spLR(spLR>0));
        validind = (sum(spLR==0, 2)==0); 
        ind = sub2ind(size(bndinfo(f).occludes), 
        nequal = bndinfo(f).occludes(spLR
        nequal = splabs(:, 1)~=splabs(:, 2) & splabs(:, 1)>0 & splabs(:, 2)>0;                
        neind = find(nequal);
    else
        nequal = false(nadj, 1);
        neind = [];
    end
    

    elabels = zeros(nadj, 3); % denotes type1, type2, relationship
    bndinfo(f).eRelations = {'leftOf', 'upperLeftOf', 'above', 'upperRightOf'};
 
    for k = 1:nadj
        
        if nequal(k)
            t = types(splabs(k, 1:2));
        else
            t = [0 0]; 
        end
                        
        rel = 0; % relationship
        
        if numel(edges{k})>2  % get edge orientation            
            orient = getEdgeOrientation(ejunctions{k}, imsize, edges{k});
            if abs(orient)> pi/4+pi/8 % leftOf
                rel = 1;  nx = -1;  ny = 0;
            elseif orient > pi/8 % upperLeftOf
                rel = 2;  nx = -1;  ny = -1;
            elseif orient > -pi/8 % above
                rel = 3;  nx = 0; ny = -1;
            else % upper-right of 
                rel = 4;  nx = 1; ny = -1;
            end

            [ey, ex] = ind2sub(imsize(1:2), edges{k});  
            pixr1 = wseg(ey+ny + (ex+nx-1)*imsize(1));
            r1 = mode(double(pixr1(pixr1>0)));
            if r1 == adjlist(k, 2)
                t = t([2 1]); % switch order of types
            elseif ~any(r1==adjlist(k, :)) % not sure of orientation
                rel = 0;
            end
        end
                
        if rel==0 % determine rel based on centers
            [spy1, spx1] = find(wseg==adjlist(k, 1));
            [spy2, spx2] = find(wseg==adjlist(k, 2));
            orientc = mod(atan2(-(median(spy2)-median(spy1)), ...
                median(spx2)-median(spx1)), 2*pi);
            if orientc > 15*pi/8
                minind = 1;
            else
                [mindist, minind] = min(abs((0:7)*pi/4 - orientc));
            end
            ind2rel = [1 4 3 2 1 4 3 2]; %(left up-left up, up-right)
            rel = ind2rel(minind);
            if any(minind==[2 3 4 5])
                t = t([2 1]);
            end
        end
            
        elabels(k, :) = [t rel];

    end
        
    bndinfo(f).elabels = elabels;
    
    if 0 
        eim = zeros(imsize(1:2));
        for k = 1:size(elabels,1)
            if elabels(k, 1)>0
                eim(edges{k}) = elabels(k,1) + 6*(elabels(k,2)-1);
            end
        end
        labim = zeros(size(eim));
        labim(wseg>0) = labels(wseg(wseg>0));
        figure(3), imagesc(labim), axis image, colormap jet
        figure(4), imagesc(eim), axis image, colormap jet
    end
    
end
        
