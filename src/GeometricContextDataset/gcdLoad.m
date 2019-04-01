gcd_dir = fullfile('resources', 'GeometricContextDataset');
classifiers_dir = fullfile('resources', 'ClassifierData');

imdir = fullfile(gcd_dir, 'images');
pbdir = fullfile(gcd_dir, 'pbim');
gcdir = fullfile(gcd_dir, 'gconf');
gtdir = fullfile(gcd_dir, 'gtsave');
segdir = fullfile(gcd_dir, 'segs');
gdatadir = fullfile(gcd_dir, 'geomdata2');
smallsegdir = fullfile(gcd_dir, 'smallsegs');

if ~exist(imdir, 'dir')
    error('Directory \"%s\" does not exist.\n', imdir);
end

load(fullfile(imdir, 'allimsegs2.mat'));
load(fullfile(imdir, 'rand_indices.mat'));
fn = {imsegs(:).imname};

generatePbFile(fn, imdir, pbdir);

tmp = load(fullfile(gcd_dir, 'data', 'classifiers_08_22_2005.mat'));
APPtestDirectoryNew(tmp.segment_density, tmp.vert_classifier, tmp.horz_classifier, imdir, fn, gcdir);
clear tmp;

exist_file_names = dir(fullfile(gtdir, '*_gt.mat'));
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
