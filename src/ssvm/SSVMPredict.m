function predictScore = SSVMPredict(X, Y, modelFile)
    %predictSSVM - Description
    %
    % Syntax: predictionFilePath = predictSSVM(input)
    %
    % Long description

    testFile = [tempname, '.txt'];
    ssvmpredict = 'lib/svm_rank/build/svm_rank_classify';
    predictFile = [tempname, '.txt'];

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

    cmd = sprintf('%s %s %s %s', ssvmpredict, testFile, modelFile, predictFile);
    fprintf('Start predicting...\n');
    disp(cmd);
    system(cmd);
    fprintf('End of predicting\n\n\n');

    predictf = fopen(predictFile, 'r');
    predictScore = fscanf(predictf, '%f');
    fclose(predictf);
    delete(testFile);
    delete(predictFile);
end
