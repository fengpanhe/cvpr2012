function result = evaluateSegmentationError(bndinfo, wseg2)

wseg1 = bndinfo.wseg;

labim = zeros(size(wseg1));
labim(wseg1>0) = bndinfo.labels(wseg1(wseg1>0));

stats = regionprops(wseg2, 'PixelIdxList', 'Area');
spind = {stats(:).PixelIdxList};
area = [stats.Area];

lab2 = zeros(max(wseg2(:)), 1);

for k = 1:numel(lab2)
    pixlab = labim(spind{k});
    pixlab = pixlab(pixlab>0);
    if numel(pixlab)>numel(spind{k})*0.01
        lab2(k) = mode(pixlab);        
    end
end

nobjOriginal = numel(setdiff(unique(bndinfo.labels), 0));
nobjNew = numel(setdiff(unique(lab2), 0));

result.efficiency = log(nobjNew / numel(unique(wseg2(:)))) / log(2);
%result.efficiency = log(nobjNew / sum(area>0.005)) / log(2);

labim2 = zeros(size(wseg2));
labim2(wseg2>0) = lab2(wseg2(wseg2>0));

ind = (labim2 > 0) & (labim > 0);

result.conservation = mean(labim(ind)==labim2(ind));
result.nregions = numel(lab2);
result.nobjOriginal = nobjOriginal;
result.nobjNew = nobjNew;
    
