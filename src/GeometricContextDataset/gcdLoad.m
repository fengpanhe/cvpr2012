gcd_dir = fullfile('resources', 'GeometricContextDataset');

data_dir = fullfile('result', 'gcd', 'data');

images_dir = fullfile(gcd_dir, 'images');
pb_dir = fullfile(gcd_dir, 'pbim');
gc_dir = fullfile(gcd_dir, 'gconf');
gt_dir = fullfile(gcd_dir, 'gtsave');

segdir = fullfile(data_dir, 'segs');
gdatadir = fullfile(data_dir, 'geomdata2');
smallsegdir = fullfile(data_dir, 'smallsegs');

if ~exist(images_dir, 'dir')
    error('Directory \"%s\" does not exist.\n', images_dir);
end

load(fullfile(images_dir, 'allimsegs2.mat'));
load(fullfile(images_dir, 'rand_indices.mat'));
fn = {imsegs(:).imname};

preparePbDir(fn, images_dir, pb_dir);
prepareSegsDir(fn, images_dir, pb_dir, segdir);
prepareGdataDir(fn, imsegs, gdatadir);

tmp = load(fullfile(gcd_dir, 'data', 'classifiers_08_22_2005.mat'));
APPtestDirectoryNew(tmp.segment_density, tmp.vert_classifier, tmp.horz_classifier, images_dir, fn, gc_dir);
clear tmp;

exist_file_names = dir(fullfile(gt_dir, '*_gt.mat'));
exist_file_names = {exist_file_names.name};
exist_file_names = strtok(exist_file_names, '_');
exist_file_names_num = numel(exist_file_names);
clusterind = zeros(1, exist_file_names_num);
for i = 1 : exist_file_names_num
    name = [exist_file_names{i} '.jpg'];
    clusterind(i) = find(fn == string(name));
end
% clusterind = cluster_images;
clear cluster_images;
clear imsegs;

if 0
    load('./data/bndClassifiers2.mat');
    load('./data/trainData2.mat');
    load('./data/contClassifier2.mat');
end
