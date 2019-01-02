imname = 'house05';

im = imread(['resources/GeometricContextDataset/images/' imname '.jpg']);
[bndinfo, pbim, gconf, bndinfo_all] = im2boundariesTopLevel(double(im) / 255);
saldir = './result/';
savePath = [saldir imname '_seg.mat'];
% save(savePath,'bndinfo', 'pbim', 'gconf', 'bndinfo_all');
iccvWriteContourDepthResults();
