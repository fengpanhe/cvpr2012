function modelFilePath = SSVMTrain(X, Y, ssvm_model_file_path, ssvmlean_param)
    %trainSSVM - Description
    %
    % Syntax: modelFilePath = trainSSVM(X, Y)
    %
    % Long description

    if ~exist('ssvmlean_param', 'var') || isempty(ssvmlean_param)
        ssvmlean_param = '-c 20.0';
    end

    if ~exist('ssvm_model_file_path', 'var') || isempty(ssvm_model_file_path)
        ssvm_model_file_path = fullfile('resources', 'SSVMmodel', 'ssvm.model');
    end

    trainFile = 'result/tmp/ssvmtrainfile.txt';

    formatStr = '%d qid:%d';

    for index = 1:(size(X, 2) - 1)
        formatStr = sprintf('%s %d:%s', formatStr, index, '%f');
    end

    formatStr = [formatStr '\n'];

    trainf = fopen(trainFile, 'w');
    fprintf(trainf, formatStr, transpose([Y, X]));
    fclose(trainf);

    ssvmlearn = strcat('lib/svm_rank/build/svm_rank_learn', ssvmlean_param);
    modelFilePath = ssvm_model_file_path;
    cmd = [ssvmlearn ' ' trainFile ' ' modelFilePath];
    fprintf('Start training...\n');
    disp(cmd);
    system(cmd);
    fprintf('End of training\n\n\n')
end
