function modelFilePath = SSVMTrain(X, Y)
    %trainSSVM - Description
    %
    % Syntax: modelFilePath = trainSSVM(X, Y)
    %
    % Long description
    trainFile = 'result/tmp/ssvmtrainfile.txt';

    formatStr = '%d qid:%d';

    for index = 1:(size(X, 2) - 1)
        formatStr = sprintf('%s %d:%s', formatStr, index, '%f');
    end

    formatStr = [formatStr '\n'];

    trainf = fopen(trainFile, 'w');
    fprintf(trainf, formatStr, transpose([Y, X]));
    fclose(trainf);

    ssvmlearn = 'lib/svm_rank/build/svm_rank_learn -c 1';
    modelFilePath = 'resources/SSVMmodel/ssvm.model';
    cmd = [ssvmlearn ' ' trainFile ' ' modelFilePath];
    fprintf('Start training...\n');
    disp(cmd);
    system(cmd);
    fprintf('End of training\n\n\n')
end
