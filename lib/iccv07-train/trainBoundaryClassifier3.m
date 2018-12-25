function classifier = trainBoundaryClassifier3(X, Y, type)

w = [];
if lower(type(1)) == 'c' % contour or not
    [tx, ty, categoryFeatures] = reformatData(X, Y);                
    %w = (ty==0) / sum(ty==0) + (ty>0) / sum(ty>0);
    %w = w/sum(w);
elseif lower(type(1)) == 's' % type of subcontour
    [tx, ty, categoryFeatures] = reformatDataSub(X, Y);   
else
    error(['Valid types: c = contour, s = contour subclass.  Given: ' type]);
end

% classifier = train_boosted_dt_2c(tx, categoryFeatures, ty*2-1, 20, 20);
nnodes = 16;
ntrees = 20;
classifier = train_boosted_dt_mc(tx, categoryFeatures, ty, ntrees, nnodes, 0, w);

classifier.prior = hist(ty, unique(ty));
classifier.prior = classifier.prior / sum(classifier.prior);


%% Reformat the data for classifying between boundary and no boundary
function [tx, ty cf] = reformatData(X, Y)


if 0  % original method for classifying among boundary types
    X = cat(1, X{:});
    Y = cat(1, Y{:});

    ind = Y(:, 1)==-1;
    X(ind, :) = [];
    Y(ind, :) = [];

    [tx, cf] = getBoundaryClassifierFeatures2(X, 'c');
    ty = Y(:,1);
end
    
if 1 % new method: each edge could be 0, foreground left, foreground right
    for f = 1:numel(X)    
        ind = (Y{f}(1:end/2)==-1) | (Y{f}(end/2+1:end)==-1);
        X{f}(ind, :) = [];
        ind = [ind ; ind];
        Y{f}(ind) = [];  
        [tx{f}, cf] = getBoundaryClassifierFeatures2(X{f}, 'c');
        ty{f} = (Y{f}(1:end/2)>0) + 2*(Y{f}(end/2+1:end)>0);
    end
    tx = cat(1, tx{:});
    ty = cat(1, ty{:});
end

hy = hist(ty, [0:2]);

disp(['Num data: ' num2str(hy)])

if 0
maxpoints = 250000;
ndata = numel(ty);
if ndata > maxpoints    
    step = floor(ndata / maxpoints);
    disp(['Reducing from ' num2str(ndata) ' to ' num2str(floor((ndata-1)/step)+1) ' points.']);    
    ty = ty(1:step:end);
    tx = tx(1:step:end, :);
end
end




%% Reformat the data for classifying the type of boundary
function [tx, ty, cf] = reformatDataSub(X, Y)


X = cat(1, X{:});
Y = cat(1, Y{:});

ind = Y(:, 1)<=0; % ignore non-contour labels
X(ind, :) = [];
Y(ind, :) = [];

ndata = size(X, 1);

[tx, cf] = getBoundaryClassifierFeatures(X, 's');

ty = false([size(tx, 1) 1]);
nc = sqrt(size(tx, 1)/ndata);

m = 0;
for y1 = 1:nc
    for y2 = 1:nc
        ty(m+1:m+ndata, :) = (Y(:, 1)==y1) & (Y(:, 2)==y2);
        m = m + ndata;
    end
end
