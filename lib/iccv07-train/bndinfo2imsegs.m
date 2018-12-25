function imsegs2 = bndinfo2imsegs(bndinfo, imsegs1)

tmpseg.segimage = bndinfo.wseg;
tmpseg.imname = bndinfo.imname;

stats = regionprops(bndinfo.wseg, 'Area');
tmpseg.npixels = vertcat(stats.Area);
tmpseg.nseg = numel(tmpseg.npixels);
tmpseg.imsize = bndinfo.imsize(1:2);

imsegs2 = APPtransferLabels(imsegs1, tmpseg);

imsegs2 = orderfields(imsegs2);
