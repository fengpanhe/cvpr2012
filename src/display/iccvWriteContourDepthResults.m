%{
% Filename: iccvWriteContourDepthResults.m
% Project: display
% Created Date: Tuesday January 22nd 2019
% Author: Feng Panhe
% -----
% Last Modified:
% Modified By:
% -----
% Copyright (c) 2019 Feng Panhe
%}

function iccvWriteContourDepthResults(file_name_list, ex_bndinfo_dir_path, gt_bndinfo_dir_path, image_dir_path, dest_dir_path)
    %myFun - Draw a depth image on the original image
    %
    % Syntax:  iccvWriteContourDepthResults(fileNameList)
    %
    % 是 11年iccv的 iccvWriteContourDepthResults 的修改。
    % 加入了面的深度显示，颜色越深，面越深。
    %
    % Inputs:
    %   file_name_list: The names of some pictures to be drawn
    %   ex_bndinfo_dir_path: The directory path of bndinfo information of experiment result
    %   gt_bndinfo_dir_path:
    %   image_dir_path: The directory path of the original image
    %   dest_dir_path: Result directory path

    ex_draw = false;

    if exist('ex_bndinfo_dir_path', 'var') &&~isempty(ex_bndinfo_dir_path)
        ex_draw = true;
    end

    if ~exist('gt_bndinfo_dir_path', 'var') || isempty(gt_bndinfo_dir_path)
        gt_bndinfo_dir_path = fullfile('resources', 'GeometricContextDataset', 'gtsave');
    end

    if ~exist('image_dir_path', 'var') || isempty(image_dir_path)
        image_dir_path = fullfile('resources', 'GeometricContextDataset', 'images');
    end

    if ~exist('dest_dir_path', 'var') || isempty(dest_dir_path)
        dest_dir_path = fullfile('result');
    end

    for i = 1:numel(file_name_list)

        file_name = strtok(file_name_list{i}, '.');
        im = im2double(imread(fullfile(image_dir_path, [file_name '.jpg'])));
        % gt
        gt_bndinfo_file = fullfile(gt_bndinfo_dir_path, [file_name '_gt.mat']);

        if exist(gt_bndinfo_file, 'file')
            load(gt_bndinfo_file, 'bndinfo');
            out_name = fullfile(dest_dir_path, [file_name '_gt.jpg']);

            [bndinfo2, ~] = updateBoundaryInfo2(bndinfo, bndinfo.labels);
            bndinfo2.names = bndinfo.names;
            bndinfo2.type = bndinfo.type;
            bndinfo2 = orderfields(bndinfo2);

            bndinfo2 = processGtBoundaryLabels(bndinfo2);

            lab = bndinfo2.edges.boundaryType;
            lab = (lab(1:end / 2) > 0) + 2 * (lab(end / 2 + 1:end) > 0);
            printOcclusionResult(rgb2gray(im), bndinfo2, lab, out_name, 1);
        else
            error('File \"%s\" does not exist.', gt_bndinfo_file);
        end

        if ex_draw
            ex_bndinfo_file = fullfile(ex_bndinfo_dir_path, [file_name '_ex.mat']);

            if exist(ex_bndinfo_file, 'file')
                load(ex_bndinfo_file, 'bndinfo');
                lab = bndinfo.edges.boundaryType;
                lab = lab(1:end / 2) + 2 * lab(end / 2 + 1:end);
                out_name = fullfile(dest_dir_path, [file_name '_ex.jpg']);
                printOcclusionResult(im, bndinfo, lab, out_name, 1);
            else
                error('File \"%s\" does not exist.', ex_bndinfo_file);
            end

        end

    end

end
