function classifier = trainBoundaryClassifier(X, Y, type)

if lower(type(1)) == 'c' % contour or not
    [tx, ty, categoryFeatures] = reformatData(X, Y);    
elseif lower(type(1)) == 's' % type of subcontour
    [tx, ty, categoryFeatures] = reformatDataSub(X, Y);
else
    error(['Valid types: c = contour, s = contour subclass.  Given: ' type]);
end

classifier = train_boosted_dt_2c(tx, categoryFeatures, ty*2-1, 20, 20);


%% Reformat the data for classifying between boundary and no boundary
function [tx, ty cf] = reformatData(X, Y)

X = cat(1, X{:});
Y = cat(1, Y{:});

ind = Y(:, 1)==-1;
X(ind, :) = [];
Y(ind, :) = [];

[tx, cf] = getBoundaryClassifierFeatures(X, 'c');

ty = Y(:,1)>0;




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
