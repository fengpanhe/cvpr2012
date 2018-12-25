function err = getSegmentationError(seg1, seg2)
% err = getSegmentationError(seg1, seg2, direction)
%
% Returns the oversegmentation error of seg2 with respect to seg1.  
%% Not implemented!!
hists = zeros(max(seg2(:)), max(seg1(:)), 'uint32');

% assign each segment in seg2 to the most common segment in seg1
for k = 1:numel(seg2)
    



labim = bndinfo.labels(wseg1);


spind = {stats(:).PixelIdxList};

lab2 = zeros(max(wseg2(:)), 1);

for k = 1:numel(lab2)
    pixlab = labim(spind{k});
    pixlab = pixlab(pixlab>0);
    if numel(pixlab)>numel(spind{k})*0.01
        lab2(k) = mode(pixlab);        
    end
end

bndinfo2 = bndinfo;
bndinfo2.wseg = uint16(wseg2);
bndinfo2.labels = lab2;

labim2 = zeros(size(wseg2));
labim2(wseg2>0) = bndinfo2.labels(wseg2(wseg2>0));

ind = (labim2 > 0) & (labim > 0);

disp(num2str([max(wseg1(:)) max(wseg2(:))]))
disp(['Error: ' num2str(mean(labim(ind)~=labim2(ind)))]);
    