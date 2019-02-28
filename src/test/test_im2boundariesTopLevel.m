imname = 'sd3';
image_file = strcat('resources/SyntheticDataset/images/', imname, '.jpg');

im = imread(image_file);
[bndinfo, pbim, gconf, bndinfo_all] = im2boundariesTopLevel(double(im) / 255);
bndinfo_file = fullfile('result', strcat(imname, '_ex.mat'));
save(bndinfo_file, 'bndinfo', 'pbim', 'gconf', 'bndinfo_all');

dest_dir_path = fullfile('result');
writeContourDepthResults(bndinfo_file, image_file, dest_dir_path, 1);
