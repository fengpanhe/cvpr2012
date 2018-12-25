function [tx, categoryFeatures] = getBoundaryClassifierFeatures(X, type)

type = lower(type(1));

switch type
    
    % boundary/noboundary classifier
    case 'c'
        ng = 5; % five geometric classes
        ndata = size(X, 1);

        tx = zeros([ndata 26], 'single');
        tx(:, 1:21) = X(:, 1:21);
        tx(:, 22) = sum(abs(X(:, 1:ng)-X(:, 2*ng+1:2*ng+ng)), 2);
        tx(:, 23) = mean([max(X(:, ng+1:ng+ng), [], 2)  ...
            max(X(:, 3*ng+1:3*ng+ng), [], 2)], 2);

        [maxval1, maxlab1] = max([X(:, 1) sum(X(:, 2:4), 2)  X(:, 5)], [], 2);
        [maxval2, maxlab2] = max([X(:, ng*2+1) ...
            sum(X(:, ng*2+(2:4)), 2)  X(:, ng*2+5)], [], 2);
        tx(:, 24:26) = [maxlab1 maxlab2 X(:, 22)]; 

        categoryFeatures = [24 25 26];
        
    % classifies type of boundary
    case 's'
        
        nc = 6; % six classes of region
        ng = 5; % five geometric classes

        c2g = [1 5 2 3 4 4]; % map region classes to geometric labels

        ndata = size(X, 1);        
        
        tx = zeros([ndata*nc*nc 13], 'single');
        
        load('./data/contourPriorProb.mat');
        
        m = 0;
        for y1 = 1:nc
            g1 = c2g(y1);
            for y2 = 1:nc
                g2 = c2g(y2);
                tx(m+1:m+ndata, 1:6) = X(:, [g1 ng+g1 2*ng+g2 3*ng+g2 4*ng+1 4*ng+2]);
                tx(m+1:m+ndata, 7:8) = repmat([y1 y2], [ndata 1]);
                
                labelOrientIndex = y1 + (y2-1)*nc + nc*nc*max((X(:, 22)-1),0);
                
                tx(m+1:m+ndata, 9) = contourPriorProb(labelOrientIndex);
                
                m = m + ndata;
            end
        end

        tx(:, 10:13) = [mean(tx(:,[1 3]),2) tx(:,1).*tx(:,3) ...
            min(tx(:,[2 4]), [], 2) max(tx(:, [2 4]), [], 2)];          
        
        categoryFeatures = [7 8];
        
    otherwise        
        error(['invalid type entered: ' type]);          
end