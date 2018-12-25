function bndinfo2 = processAllBoundaryInfo(bndinfo, doLabels)

for f = 1:numel(bndinfo)
    disp(num2str(f))
    bndinfo2(f) = processBoundaryInfo2(bndinfo(f));
end

if doLabels
    bndinfo2 = processGtBoundaryLabels(bndinfo2);
end
    