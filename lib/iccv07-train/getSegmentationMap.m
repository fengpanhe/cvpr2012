function [smap, regions] = getSegmentationMap(nsp, result)
%
% [smap, regions] = getSegmentationMap(nsp, result)
%
% Returns a segmentation map based on one or more "result" structures from
% mergeStage

regions = result(end).regions;
for s = numel(result)-1:-1:1
    prevregions = result(s).regions;    
    regions2 = cell(size(regions));
    for k = 1:numel(regions)
        for k2 = regions{k}
            regions2{k} = [regions2{k} prevregions{k2}];
        end
    end
    regions = regions2;
end

smap = zeros(nsp, 1);
for k = 1:numel(regions)
    smap(regions{k}) = k;
end