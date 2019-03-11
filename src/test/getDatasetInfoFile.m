function [gt_file, image_file, pbim_file] = getDatasetInfoFile(dataset_path, image_name)
    %getDatasetInfoFile - Description
    %
    % Syntax: [gt_file, image_file, pbim_file] = getDatasetInfoFile(dataset_path, image_name)
    %
    % Long description

    pbim_dir_path = fullfile(dataset_path, 'pbim');

    if ~exist(pbim_dir_path, 'dir')
        mkdir(pbim_dir_path);
    end

    gt_file = fullfile(dataset_path, 'gtsave', strcat(image_name, '_gt.mat'));
    image_file = fullfile(dataset_path, 'images', strcat(image_name, '.jpg'));
    pbim_file = fullfile(dataset_path, 'pbim', strcat(image_name, '_pbim.mat'));

    if ~exist(gt_file, 'file')
        error('File \"%s\" does not exist.', gt_file);
    end

    if ~exist(image_file, 'file')
        error('File \"%s\" does not exist.', image_file);
    end

    if ~exist(pbim_file, 'file')
        image_f = imread(image_file);
        pbim = pbCGTG_nonmax(double(image_f) / 255);
        save(pbim_file, 'pbim');
    end

end
