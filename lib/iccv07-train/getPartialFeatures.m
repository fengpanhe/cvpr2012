function [tx2, catfeatures2] = getPartialFeatures(tx, catfeatures, featureset)

catfeatures2 = [];
for k = 1:numel(catfeatures)
    if ismember(catfeatures(k), featureset)
        catfeatures2(end+1) = sum(featureset<catfeatures(k))+1;
    end
end
tx2 = tx(:, featureset);
