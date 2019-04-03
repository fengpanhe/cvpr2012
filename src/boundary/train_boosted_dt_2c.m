function classifier = train_boosted_dt_2c(features, cat_features, ...
        labels, num_iterations, nodespertree, stopval, w)
    classifier = fitcensemble(features, labels, 'Method', 'AdaBoostM1', 'Learners', 'tree', 'HyperparameterOptimizationOptions', struct('UseParallel', true));
end
