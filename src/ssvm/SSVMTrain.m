function modelFilePath = SSVMTrain(X, Y)
%trainSSVM - Description
%
% Syntax: modelFilePath = trainSSVM(X, Y)
%
% Long description
    trainFile = 'result/tmp/ssvmtrainfile.txt';
    
    formatStr = '%d qid:%d';
    for index = 1:(size(X, 2) - 1)
        formatStr = formatStr + ' ' + string(index) + ':%f'; 
    end
    formatStr = formatStr + '\n';

    trainf = fopen(trainFile，'w');
    for index = 1:numel(X)
        fprintf(trainf, formatStr, [Y, X]);
    end
    fclose(trainf);
    
    ssvmlearn = 'lib/svm_rank/build/svm_rank_learn svm_learn -z p -c 1';
    modelFilePath = 'resources/SSVMmodel/ssvm.model';
    cmd = [ssvmlearn ' ' trainFile ' ' modelFilePath];
    disp(cmd);
    system(cmd);
end