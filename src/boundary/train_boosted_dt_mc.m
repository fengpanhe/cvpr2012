function classifier = train_boosted_dt_mc(features, cat_features, labels, ...
        num_iterations, num_nodes, stopval, init_weights, varargin)
    %
    classifier = fitcensemble(features, labels, 'Method', 'AdaBoostM2', 'Learners', 'tree', 'Weights', init_weights, 'HyperparameterOptimizationOptions', struct('UseParallel', true));
end
