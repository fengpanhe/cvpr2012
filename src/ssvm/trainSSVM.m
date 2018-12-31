function modelFilePath = trainSSVM(trainFile)
%trainSSVM - Description
%
% Syntax: modelFilePath = trainSSVM(X, Y)
%
% Long description
    ssvmlearn = 'lib/svm_rank/build/svm_rank_learn svm_learn -z p -c 1';
    modelFilePath = 'resources/SSVMmodel/ssvm.model';
    cmd = [ssvmlearn ' ' trainFile ' ' modelFilePath];
    disp(cmd);
    system(cmd);
end