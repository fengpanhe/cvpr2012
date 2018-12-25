% iccvLoad

basedir = '/home/dhoiem/src/cmu';
imdir = fullfile('../GeometricContext/images/all_images/');
pbdir = fullfile(basedir, 'iccv07/data/pb_nonmax2/');
gcdir = fullfile(basedir, 'iccv07/data/gconf/');
gtdir = fullfile(basedir, 'iccv07/iccvGroundTruth/gtsave/');
segdir = fullfile(basedir, 'iccv07/iccvGroundTruth/segs/');
gdatadir = fullfile(basedir, 'iccv07/data/geomdata2/');
smallsegdir = fullfile(basedir, 'iccv07/data4/smallsegs/');

load('../GeometricContext/data/allimsegs2.mat');
load('../GeometricContext/data/rand_indices.mat');

fn = {imsegs(:).imname};
clear imsegs;

clusterind = cluster_images;
clear cluster_images;

if 0
load('./data/bndClassifiers2.mat');
load('./data/trainData2.mat');
load('./data/contClassifier2.mat');
end