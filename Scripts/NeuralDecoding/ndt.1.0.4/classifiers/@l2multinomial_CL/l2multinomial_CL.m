classdef l2multinomial_CL
% NDT-compatible L2-regularized multinomial logistic classifier.
%
% NDT interface:
%   cl = cl.train(XTr, YTr)
%   [predicted_labels, decision_values] = cl.test(XTe)
%
% XTr/XTe are features x observations. decision_values are
% observations x classes. Requires Optimization Toolbox (fminunc).
%
% For strict nested CV, do not also pass zscore_normalize_FP: this class
% standardizes within each inner fold and refits scaling on all outer
% training data.

    properties
        lambda_grid = logspace(-4, 2, 7)
        num_inner_cv_splits = 4
        inner_cv_seed = 1
        inner_cv_metric = 'balanced_accuracy'
        standardize_features = false

        max_iterations = 500
        optimality_tolerance = 1e-6
        step_tolerance = 1e-9
        optimizer_display = 'iter'

        selected_lambda = []
        class_labels = []
        weights = []
        intercepts = []
        feature_means = []
        feature_stds = []
        inner_cv_scores = []
    end

    methods
        function cl = l2multinomial_CL(varargin)
            if mod(numel(varargin), 2) ~= 0
                error('Constructor inputs must be name-value pairs.');
            end
            for iArg = 1:2:numel(varargin)
                name = char(varargin{iArg});
                if ~isprop(cl, name)
                    error('Unknown property: %s', name);
                end
                cl.(name) = varargin{iArg + 1};
            end
        end

        function cl = train(cl, XTr, YTr)
            if ~isnumeric(XTr) || ndims(XTr) ~= 2
                error('XTr must be a numeric features x observations matrix.');
            end
            YTr = YTr(:);
            if size(XTr, 2) ~= numel(YTr)
                error('Columns of XTr must equal the number of labels.');
            end
            if any(~isfinite(XTr(:)))
                error('XTr contains NaN or Inf.');
            end

            cl.class_labels = unique(YTr);
            nClasses = numel(cl.class_labels);

            [~, yIndex] = ismember(YTr, cl.class_labels);

            lambdaGrid = cl.lambda_grid(:);
            classCounts = accumarray(yIndex, 1, [nClasses, 1]);

            if numel(lambdaGrid) == 1
                cl.selected_lambda = lambdaGrid(1);
                cl.inner_cv_scores = NaN;
            else
                nFolds = min(cl.num_inner_cv_splits, min(classCounts));
                foldId = cl.make_stratified_fold_ids( ...
                    yIndex, nFolds, cl.inner_cv_seed);
                scores = nan(numel(lambdaGrid), nFolds);

                for iLambda = 1:numel(lambdaGrid)
                    for iFold = 1:nFolds
                        isVal = foldId == iFold;

                        [XInnerTrain, XInnerVal] = cl.standardize_train_and_test( ...
                                XTr(:, ~isVal), XTr(:, isVal));

                        [W, b] = cl.fit_softmax_model( ...
                            XInnerTrain, yIndex(~isVal), nClasses,lambdaGrid(iLambda));

                        [pred, ~] = cl.predict_indices(XInnerVal, W, b);
                        scores(iLambda, iFold) = cl.calculate_metric( ...
                            yIndex(isVal), pred, nClasses);
                    end
                end

                meanScore = mean(scores, 2, 'omitnan');
                [~,bestIndex] = randmax(find(abs(meanScore - max(meanScore)) <= 1e-12));

                cl.selected_lambda = lambdaGrid(bestIndex);
                cl.inner_cv_scores = scores;
            end

            if cl.standardize_features
                [cl.feature_means, cl.feature_stds] = cl.fit_standardization(XTr);
                XFit = cl.apply_standardization(XTr, cl.feature_means, cl.feature_stds);
            else
                cl.feature_means = zeros(size(XTr, 1), 1);
                cl.feature_stds = ones(size(XTr, 1), 1);
                XFit = XTr;
            end

            [cl.weights, cl.intercepts] = cl.fit_softmax_model( ...
                XFit, yIndex, nClasses, cl.selected_lambda);
        end

        function [predicted_labels, decision_values] = test(cl, XTe)
            if size(XTe, 1) ~= size(cl.weights, 1)
                error('XTe has a different number of features than XTr.');
            end
            if any(~isfinite(XTe(:)))
                error('XTe contains NaN or Inf.');
            end

            if cl.standardize_features
                XTe = cl.apply_standardization(XTe, cl.feature_means, cl.feature_stds);
            end

            [predictedIndex, decision_values] = cl.predict_indices(XTe, cl.weights, cl.intercepts);
            predicted_labels = cl.class_labels(predictedIndex);
            predicted_labels = predicted_labels(:);
        end
    end

    methods (Access = private)
        function [XTrainZ, XTestZ] = standardize_train_and_test(cl, XTrain, XTest)
            if cl.standardize_features
                [mu, sigma] = cl.fit_standardization(XTrain);
                XTrainZ = cl.apply_standardization(XTrain, mu, sigma);
                XTestZ = cl.apply_standardization(XTest, mu, sigma);
            else
                XTrainZ = XTrain;
                XTestZ = XTest;
            end
        end

        function [mu, sigma] = fit_standardization(~, X)
            mu = mean(X, 2);
            sigma = std(X, 0, 2);
            mu(~isfinite(mu)) = 0;
            sigma(~isfinite(sigma) | sigma < 1e-12) = 1;
        end

        function Xz = apply_standardization(~, X, mu, sigma)
            Xz = bsxfun(@rdivide, bsxfun(@minus, X, mu), sigma);
        end

        function [W, b] = fit_softmax_model(cl, X, y, nClasses, lambda)
            [nFeatures, ~] = size(X);
            nModeledClasses = nClasses - 1;
            nWeightParameters = nFeatures * nModeledClasses;
            theta0 = zeros(nWeightParameters + nModeledClasses, 1);

            objective = @(theta) cl.softmax_objective(theta, X, y, nClasses, lambda);

            options = optimoptions('fminunc','Algorithm', 'quasi-newton', ...
                'SpecifyObjectiveGradient', true,'Display', cl.optimizer_display, ...
                'MaxIterations', cl.max_iterations,'OptimalityTolerance', cl.optimality_tolerance, ...
                'StepTolerance', cl.step_tolerance);

            theta = fminunc(objective, theta0, options);

            W = reshape(theta(1:nWeightParameters),nFeatures, nModeledClasses);
            b = reshape(theta(nWeightParameters + 1:end),1, nModeledClasses);
        end

        function [loss, gradient] = softmax_objective(~, theta, X, y, nClasses, lambda)
            [nFeatures, nObs] = size(X);
            nModeledClasses = nClasses - 1;
            nWeightParameters = nFeatures * nModeledClasses;

            W = reshape(theta(1:nWeightParameters),nFeatures, nModeledClasses);
            b = reshape(theta(nWeightParameters + 1:end),1, nModeledClasses);

            modeledScores = X.' * W + repmat(b, nObs, 1);
            allScores = [modeledScores, zeros(nObs, 1)];

            rowMax = max(allScores, [], 2);
            expScores = exp(bsxfun(@minus, allScores, rowMax));
            probabilities = bsxfun(@rdivide,expScores, sum(expScores, 2));

            idx = sub2ind([nObs, nClasses], (1:nObs).', y);
            pTrue = probabilities(idx);

            nll = -mean(log(max(pTrue, realmin('double'))));
            loss = nll + 0.5 * lambda * sum(W(:).^2);

            target = zeros(nObs, nModeledClasses);
            modeledRows = find(y <= nModeledClasses);
            if ~isempty(modeledRows)
                targetIdx = sub2ind(size(target), modeledRows, y(modeledRows));
                target(targetIdx) = 1;
            end

            residual = probabilities(:, 1:nModeledClasses) - target;
            gradientW = (X * residual) ./ nObs + lambda .* W;
            gradientB = sum(residual, 1) ./ nObs;
            gradient = [gradientW(:); gradientB(:)];
        end

        function [predictedIndex, probabilities] = predict_indices(~, X, W, b)
            nObs = size(X, 2);
            modeledScores = X.' * W + repmat(b, nObs, 1);
            allScores = [modeledScores, zeros(nObs, 1)];

            rowMax = max(allScores, [], 2);
            expScores = exp(bsxfun(@minus, allScores, rowMax));
            probabilities = bsxfun(@rdivide,expScores, sum(expScores, 2));

            [~, predictedIndex] = max(probabilities, [], 2);
        end

        function score = calculate_metric(cl, yTrue, yPred, nClasses)
            switch lower(cl.inner_cv_metric)
                case 'accuracy'
                    score = mean(yTrue == yPred);
                case 'balanced_accuracy'
                    recall = nan(nClasses, 1);
                    for iClass = 1:nClasses
                        mask = yTrue == iClass;
                        if any(mask)
                            recall(iClass) = mean(yPred(mask) == iClass);
                        end
                    end
                    score = mean(recall, 'omitnan');
                otherwise
                    error('Unknown inner_cv_metric: %s',cl.inner_cv_metric);
            end
        end

        function foldId = make_stratified_fold_ids(~, y, nFolds, seed)
            stream = RandStream('mt19937ar', 'Seed', seed);
            foldId = zeros(numel(y), 1);
            classes = unique(y);

            for iClass = 1:numel(classes)
                indices = find(y == classes(iClass));
                order = randperm(stream, numel(indices));
                indices = indices(order);
                folds = mod(0:(numel(indices) - 1), nFolds) + 1;
                foldId(indices) = folds(:);
            end
        end
    end
end