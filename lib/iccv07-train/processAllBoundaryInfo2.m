function bndinfo2 = processAllBoundaryInfo2(bndinfo, doLabels)

for f = 1:numel(bndinfo)
    disp(num2str(f))
    wseg = bndinfo(f).wseg;
    [fragments, neighbors, seg, polyfragments, poly_params] = seg2fragments(double(wseg), im, 25);
    bndinfo2(f) = processBoundaryInfo3(bndinfo(f));
end

if doLabels
    bndinfo2 = processGtBoundaryLabels(bndinfo2);
end
    