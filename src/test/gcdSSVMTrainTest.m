function [ssvm_precision, mrf_precision, each_res] = gcdSSVMTrainTest()
    %gcdSSVMTrainTest - Description
    %
    % Syntax: gcdSSVMTrainTest()
    %
    % Long description
    gcdGtPath = fullfile('resources', 'GeometricContextDataset', 'gtsave');
    fileList = dir(fullfile(gcdGtPath, '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');
    file_num = numel(file_name_list);
    train_file_num = round(file_num * 2/3);
    train_file_list = file_name_list(1:train_file_num);
    test_file_list = file_name_list(train_file_num + 1:end);

    model_file = gcdSSVMTrain(train_file_list);

    dataset_path = fullfile('resources', 'GeometricContextDataset');
    [ssvm_precision, mrf_precision, each_res] = multImageTJTest(dataset_path, test_file_list, model_file);
end
