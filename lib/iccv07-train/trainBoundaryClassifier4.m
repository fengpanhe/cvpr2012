function [classifier, xall, yall] = ...
    trainBoundaryClassifier4(X, Y, bndinfo, featureset)

w = [];
[tx, ty, categoryFeatures, xall, yall] = reformatData(X, Y, bndinfo);                

disp('data weighting');
w = tx(:, 10); % weight by area of smaller region
w = w ./ sum(w);
disp(num2str(sum(-w.*log(w))/log(numel(w))))

if exist('featureset', 'var') && ~isempty(featureset)
    [tx, categoryFeatures] = getPartialFeatures(tx, categoryFeatures, featureset);
end

% classifier = train_boosted_dt_2c(tx, categoryFeatures, ty*2-1, 20, 20);
nnodes = 16;
ntrees = 20;
classifier = train_boosted_dt_mc(tx, categoryFeatures, ty, ntrees, nnodes, 0, w);

classifier.prior = hist(ty, unique(ty));
classifier.prior = classifier.prior / sum(classifier.prior);


%% Reformat the data for classifying between boundary and no boundary
function [tx, ty, cf, xall, yall] = reformatData(X, Y, bndinfo)

   
% new method: each edge could be off, foreground left, or foreground right
for f = 1:numel(X)    
    ind = (Y{f}(1:end/2)==-1) | (Y{f}(end/2+1:end)==-1);
    xind = find(~ind);
    ind = [ind ; ind];
    Y{f}(ind) = [];  
    [tx{f}, cf] = getBoundaryClassifierFeatures6(bndinfo(f), X(f), xind); % was Features4
    ty{f} = (Y{f}(1:end/2)>0) + 2*(Y{f}(end/2+1:end)>0);

    xall{f} = zeros(bndinfo(f).ne, size(tx{f},2));
    yall{f} = -ones(bndinfo(f).ne, 1);
    xall{f}(xind, :) = tx{f};
    yall{f}(xind, :) = ty{f};

end

tx = cat(1, tx{:});
ty = cat(1, ty{:});


hy = hist(ty, [0:2]);

disp(['Num data: ' num2str(hy)])

if 1
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
