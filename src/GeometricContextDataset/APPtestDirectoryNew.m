% copy by lib/GeometricContext/src/APPtestDirectory.m
function APPtestDirectoryNew(segDensity, vClassifier, hClassifier, ...
        imdir, imfn, outdir)
    % [labels, conf_map] = APPtestDirectory(segDensity, vClassifier,
    %                                   hClassifier, imdir, imsegs, varargin)
    %
    % Gets the geometry for each image with superpixels given by imsegs.
    %
    % Input:
    %   segDensity: structure giving probability of 2 sp having same label
    %   vClassifier: segment classifier for ground/vert/sky
    %   hClassifier: segment classifier for subclassses of vert
    %   imdir: the image directory
    %   imfn: the image filenames
    %   varargin{1} (optional): the output directory for displaying results
    % Output:
    %   labels: structure containing labeling results for each image
    %   conf_map: the likelihoods for each sp for each class for each image
    %
    % Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
    % Current Version: 1.0  09/30/2005

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    file_names = strtok(imfn, '.');

    exist_file_names = dir(fullfile(outdir, '*.c.mat'));
    exist_file_names = {exist_file_names.name};
    exist_file_names = strtok(exist_file_names, '.');
    not_exist_file_names = setdiff(file_names, exist_file_names);

    not_exist_file_num = numel(not_exist_file_names);
    % spmd
    if not_exist_file_num > 0

        parfor f = 1:not_exist_file_num
            % if labindex ~= f % numlabs
            %     continue;
            % end
            fprintf('\n\n image %d / %d .\n', f, not_exist_file_num);
            bn = not_exist_file_names{f};
            fn = [bn '.jpg'];
            c_mat_file = fullfile(outdir, [bn '.c.mat']);

            if ~exist(c_mat_file, 'file')

                image = im2double(imread(fullfile(imdir, fn)));

                tmp = im2superpixels(image);
                tmp.imname = fn;
                imsegs = tmp;

                if size(image, 3) == 3

                    disp(['processing image ' bn]);

                    [glabels, gconf_map, gmaps, gpmaps] = ...
                        APPtestImage(image, imsegs, vClassifier, hClassifier, segDensity);
                    % for generating context
                    [cimages, cnames] = APPclassifierOutput2confidenceImages(imsegs, gconf_map);

                    limage = APPgetLabeledImage(image, imsegs, glabels.vert_labels, glabels.vert_conf, ...
                        glabels.horz_labels, glabels.horz_conf);
                    imwrite(limage, [outdir, '/', bn, '.l.jpg']);
                    imwrite(image, [outdir, '/', fn]);

                    % make a vrml file
                    APPwriteVrmlModel(imdir, imsegs, glabels, outdir);

                    c_mat_file_save(c_mat_file, glabels, gconf_map, cimages, cnames, gmaps, gpmaps, imsegs)
                end

                drawnow;

            end

        end

    end

end

function c_mat_file_save(c_mat_file, glabels, gconf_map, cimages, cnames, gmaps, gpmaps, imsegs)
    save(c_mat_file, 'glabels', 'gconf_map', 'cimages', 'cnames', 'gmaps', 'gpmaps', 'imsegs');
end
