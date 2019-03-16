function [ssvm_precision, mrf_precision, each_res] = sdSSVMTrainTest()
    %sdSSVMTrainTest - Description
    %
    % Syntax: sdSSVMTrainTest()
    %
    % Long description

    sdGtPath = fullfile('resources', 'SyntheticDataset', 'gtsave');
    fileList = dir(fullfile(sdGtPath, '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');
    file_num = numel(file_name_list);
    train_file_num = round(file_num * 2/3);
    train_file_list = file_name_list(1:train_file_num);
    test_file_list = file_name_list(train_file_num + 1:end);

    X = [];
    Y = [];
    infos = [];

    for i = 1:numel(train_file_list)
        [X_, Y_, infos_] = sd2ssvmXY(train_file_list{i});

        if size(X, 1) > 0
            X_(:, 1) = X_(:, 1) + X(end, 1);
        end

        X = [X; X_];
        Y = [Y; Y_];
        infos_(:, end + 1) = i;
        infos = [infos; infos_];
    end

    infos(:, [2, 1]) = infos(:, [1, 2]);

    % X = X(:, 1:40);

    model_file = SSVMTrain(X, Y);
    dataset_path = fullfile('resources', 'SyntheticDataset');
    [ssvm_precision, mrf_precision, each_res] = multImageTJTest(dataset_path, test_file_list, model_file);
end

