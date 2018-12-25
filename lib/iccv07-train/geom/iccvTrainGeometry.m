function gclassifiers = iccvTrainGeometry(fn, segmaps, ncv, gdatadir)

nim = numel(fn);

for tf = 1:numel(fn)

    disp(['Geometry train features: ' num2str(tf)])
        
    gdata = load([gdatadir strtok(fn{tf}, '.') '_gdata.mat']);
    
    if isempty(segmaps{tf})
        segmaps{tf} = (1:gdata.imsegs.nseg);
    end
    
    labdata{tf} = mcmcGetSegmentFeatures(gdata.imsegs, gdata.spdata, gdata.imdata, ...
        segmaps{tf}, (1:max(segmaps{tf})));
    [mclab{tf}, mcprc{tf}, allprc{tf}, trainw{tf}] = ...
        segmentation2labels(gdata.imsegs, segmaps{tf});
    unilabel{tf} = mclab{tf} .*(mcprc{tf}>0.95);         
end

save './data/tmp/geomtmp.mat';
for k = 1:ncv    
    disp(['Iteration: ' num2str(k)]);
    testind = (floor((k-1)*nim/ncv)+1):(floor(k*nim/ncv));
    trainind = setdiff((1:nim), testind);
    [vclassifier(k), hclassifier(k)] = ...
        mcmcTrainSegmentClassifier2(labdata(trainind),...
        unilabel(trainind), trainw(trainind), 25000);       
    save './data/tmp/geomtmp.mat';
end

gclassifiers.vclassifier = vclassifier;
gclassifiers.hclassifier = hclassifier;