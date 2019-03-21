function [ssvm_precision, mrf_precision, each_res] = gcdSSVMTest(model_file)
    dataset_path = fullfile('resources', 'GeometricContextDataset');
    fileList = dir(fullfile(dataset_path, 'gtsave', '*_gt.mat'));
    fileList = {fileList.name};

    file_name_list = strtok(fileList, '_');

    [ssvm_precision, mrf_precision, each_res] = multImageTJTest(dataset_path, file_name_list, model_file);
end