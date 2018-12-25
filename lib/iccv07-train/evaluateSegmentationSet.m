function results = evaluateSegmentationSet(gtdir, resdir, fn, ext)

if ~exist('ext', 'var')
    ext = '_seg';
end

for f = 1:numel(fn)

    disp(num2str(f))
    
    bn = strtok(fn{f}, '.');
    
    gt = load([gtdir bn '_gt.mat']); 
    res = load([resdir bn ext '.mat']);
       
    try
        wseg = res.bndinfo.wseg;
    catch
        wseg = res.bndinfo_sm.wseg;
    end
    if isfield(res.bndinfo, 'result') && isfield(res.bndinfo.result, 'geomProb')
        objlab = zeros(res.bndinfo.nseg, 1);
        [maxval, maxlab] = max(res.bndinfo.result.geomProb, [], 2);
        ind = find(maxlab>=2 & maxlab<=4);
        objlab(ind) = (3:2+numel(ind));
        objlab(maxlab==1) = 1;
        objlab(maxlab==5) = 2;
        wseg = objlab(wseg);
        disp(num2str([numel(unique(objlab)) res.bndinfo.nseg]))
    end
    
    tmp = evaluateSegmentationError(gt.bndinfo, wseg); 
    tmp.imname = fn{f};
    results(f) = tmp;
    
end