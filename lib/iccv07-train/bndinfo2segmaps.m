function segmaps = bndinfo2segmaps(bndinfo, segdir)

for f = 1:numel(bndinfo)
    
    fn = [segdir '/' strtok(bndinfo(f).imname, '.') '_seg.mat'];
    tmp = load(fn);
    
    bndinfo_orig = tmp.bndinfo;
    
    bndinfo(f).labels = (1:bndinfo(f).nseg);
    bndinfotmp = transferSuperpixelLabels(bndinfo(f), bndinfo_orig.wseg);
    
    segmaps{f} = bndinfotmp.labels;

end