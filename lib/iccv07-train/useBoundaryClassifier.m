function [conf, tx] = useBoundaryClassifier(bndinfo, X, classifier, featureset)
% [conf, tx] = useBoundaryClassifier(bndinfo, X, classifier, featureset)

tx = double(getBoundaryClassifierFeatures6(bndinfo, X, (1:bndinfo.ne)));

if exist('featureset', 'var') && ~isempty(featureset)
    tx = getPartialFeatures(tx, [], featureset);
end

conf = test_boosted_dt_mc(classifier, tx);
conf = 1 ./ (1+exp(-conf));
conf = conf ./ repmat(sum(conf, 2), [1 size(conf, 2)]); 

% reformat conf for subclassification
% type = lower(type(1));
% 
% if type == 's'
%     ndata = size(X, 1);
%     nc = numel(conf) / ndata;    
%     conf = reshape(conf, [ndata nc]);    
%     mean(sum(conf, 2))
%     conf = conf ./ repmat(sum(conf, 2), [1 nc]);             
%     
% end