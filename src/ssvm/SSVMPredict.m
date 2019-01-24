function predictScore = SSVMPredict(X, Y, modelFile)
    %predictSSVM - Description
    %
    % Syntax: predictionFilePath = predictSSVM(input)
    %
    % Long description

    testFile = 'result/tmp/sssvmtestfile.txt';
    ssvmpredict = 'lib/svm_rank/build/svm_rank_classify';
    predictFile = 'result/tmp/predict.txt';

    if ~exist(ssvmpredict, 'file')
        error('File \"%s\" does not exist.', ssvmpredict);
    end

    if ~exist(modelFile, 'file')
        error('File \"%s\" does not exist.', modelFile);
    end

    formatStr = '%d qid:%d';

    for index = 1:(size(X, 2) - 1)
        formatStr = sprintf('%s %d:%s', formatStr, index, '%f');
    end

    formatStr = [formatStr '\n'];

    testf = fopen(testFile, 'w');
    fprintf(testf, formatStr, transpose([Y, X]));
    fclose(testf);

    cmd = sprintf("%s %s %s %s", ssvmpredict, testFile, modelFile, predictFile);
    disp(cmd);
    system(cmd);

    predictf = fopen(predictFile, 'r');
    predictScore = fscanf(predictf, '%f');
    fclose(predictf);
end
