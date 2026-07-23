classdef gaussian_naive_bayes_CL
% gaussian_naive_bayes_CL implements an independent Gaussian
% naive Bayes classifier.
%
% For every class and feature, the classifier estimates:
%   1. The feature mean
%   2. The feature variance
%
% The classifier assumes that features are conditionally independent
% within each class. XTr and XTe have dimensions:
%       number of features x number of examples
%
% YTr has one label per training example.
    properties
        means = [];
        variances = [];
        labels = [];
    end

    methods
        function cl = gaussian_naive_bayes_CL
        end

        function cl = train(cl, XTr, YTr)
            if numel(YTr) ~= size(XTr, 2)
                error('There must be exactly one label for each training example.');
            end
            YTr = YTr(:);
            cl.labels = unique(YTr);
            num_features = size(XTr, 1);
            num_classes = length(cl.labels);
            means = zeros(num_features, num_classes);
            variances = zeros(num_features, num_classes);
            % Calculate a feature-specific minimum variance from training data
            variance_floor = max(1e-12, 1e-6 .* var(XTr, 1, 2));
            for iLabel = 1:num_classes
                means(:, iLabel) = mean(XTr(:,YTr == cl.labels(iLabel)),2);
                % Maximum-likelihood variance estimate.
                variances(:, iLabel) = var(XTr(:,YTr == cl.labels(iLabel)),1,2);
                variances(:, iLabel) = max(variances(:,iLabel), variance_floor);
            end
            cl.means = means;
            cl.variances = variances;
        end

        function [predicted_labels, decision_values] = test(cl, XTe)
            if size(XTe, 1) ~= size(cl.means, 1)
                error('The test data must have the same number of features as the training data.');
            end
            num_test_points = size(XTe, 2);
            num_classes = size(cl.means, 2);
            % Dimensions after replication: features x classes x test points
            curr_means = repmat(cl.means, [1, 1, num_test_points]);
            curr_variances = repmat(cl.variances, [1, 1, num_test_points]);
            XTe_for_all_classes = permute(repmat(XTe, [1, 1, num_classes]), ...
                [1, 3, 2]);

            % Independent gauss log-likelihood for every feature,class,test
            log_feature_likelihoods = -0.5 .* ( log(2 .* pi .* curr_variances) + ...
                ((XTe_for_all_classes - curr_means).^2) ./ curr_variances);
            % Sum across features
            log_likelihoods = reshape(sum(log_feature_likelihoods, 1), ...
                num_classes,num_test_points);

            [~, inds] = randmax(log_likelihoods);
            predicted_labels = cl.labels(inds);
            decision_values = log_likelihoods;
        end
    end
end