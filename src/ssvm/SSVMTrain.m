function modelFilePath = SSVMTrain(X, Y)
%trainSSVM - Description
%
% Syntax: modelFilePath = trainSSVM(X, Y)
%
% Long description
    trainFile = 'result/tmp/ssvmtrainfile.dat';
    
    formatStr = '%f qid:%f';
    for index = 1:(size(X, 2) - 1)
        formatStr = sprintf('%s %d:%s', formatStr, index, '%f'); 
    end
    formatStr = [formatStr '\n'];

    trainf = fopen(trainFile, 'w');
    YX = double([Y, X]);
    for index = 1:size(YX, 1)
        fprintf(trainf, formatStr, YX(1, :));
    end
    fclose(trainf);
    
    ssvmlearn = 'lib/svm_rank/build/svm_rank_learn svm_learn -z p -c 1';
    modelFilePath = 'resources/SSVMmodel/ssvm.model';
    cmd = [ssvmlearn ' ' trainFile ' ' modelFilePath];
    disp(cmd);
    system(cmd);
end