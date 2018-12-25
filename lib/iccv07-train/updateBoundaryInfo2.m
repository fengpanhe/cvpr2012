function [bndinfo2, err] = updateBoundaryInfo2(bndinfo, result)

if isstruct(result) && isfield(result, 'regions')
    regions = result.regions;

    
    empty = false(numel(regions), 1);
    splab = zeros(bndinfo.nseg, 1);
    for k = 1:numel(regions)
        splab(regions{k}) = k;
        empty(k) = isempty(regions{k});        
    end
    regions(empty) = [];
else
    splab = result;
    regions = cell(max(splab), 1);
    empty = false(numel(regions), 1);
    for k = 1:numel(regions)
        regions{k} = find(splab==k);
        empty(k) = isempty(regions{k});  
    end
    regions(empty) = [];
end
    
wseg = splab(bndinfo.wseg);

% XXX only reading image until seg2fragments gets fixed
im = imread(['/run/media/he/data/Experiment/iccv07Final/GeometricContext/test_dir/images/' bndinfo.imname]);
im = im2double(im);

[edges, juncts, neighbors, wseg] = seg2fragments(wseg, im, 1);
bndinfo2 = processBoundaryInfo3(wseg, edges, neighbors);

bndinfo2.imname = bndinfo.imname;

stats = regionprops(bndinfo.wseg, 'Area');
area = cat(1, stats(:).Area);
bndinfo2.spArea = area;

if isfield(bndinfo, 'type')
    bndinfo2.type = bndinfo.type;
    bndinfo2.names = bndinfo.names;
    [bndinfo2.labels, err] = transferLabels(bndinfo.labels, regions, area);
    bndinfo2 = processGtBoundaryLabels(bndinfo2);

else
    err = nan;
end

bndinfo2 = orderfields(bndinfo2);