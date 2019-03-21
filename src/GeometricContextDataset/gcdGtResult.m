function gcdGtResult()
    dataset_path = fullfile('resources', 'GeometricContextDataset');
    fileList = dir(fullfile(dataset_path, 'gtsave', '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');

    set(0, 'DefaultFigureVisible', 'off');

    for i = 1:numel(file_name_list)
        im_name = file_name_list{i};
        gt_file = fullfile(dataset_path, 'gtsave', strcat(im_name, '_gt.mat'));
        image_file = fullfile(dataset_path, 'images', strcat(im_name, '.jpg'));
        out_path = fullfile('result', 'gcd');
        writeContourDepthResults(gt_file, image_file, out_path, 0);
        close all;
    end

    set(0, 'DefaultFigureVisible', 'on');
end
