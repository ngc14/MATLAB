classdef euclidean_CL

    % euclidean_CL is a classifier object (CL)
    %  that learns a mean population vector (template) for each class from the
    %  training data.  When the classifier is tested, the correlation coefficient
    %  is calculated between a test point and each of the class templates, and the
    %  class with the largest correlation coefficient value is selected as the label
    %  (and the decision values are the corrcoef values).  
    %
    % Like all CL objects, there are two main methods, which are:
    %  1.  cl = train(cl, XTr, YTr)
    %         This method takes the training data (XTr, YTr) and learns a mean vector
    %           (i.e., a template) for each class.
    %  2.  [predicted_labels decision_values] = test(cl, XTe)
    %         This method takes the test data and calculates the correlation coefficient
    %           value between each test point and each learned class template.  The
    %           predicted label for a test point is the class that had the largest corrcoef
    %           value with the test point, and the decision values are the actually corrcoef values.
    %
    %   XTr and XTe are in the form [num_features x num_examples]
    %   YTr is in the form [num_examples x 1]
    %
    %==========================================================================
    %==========================================================================
    properties
        templates = [];   % average of the training vectors for each class
        labels = [];  % all the unique labels for each class (one for each template)
    end

    methods
        % constructor
        function cl = euclidean_CL
        end

        function cl = train(cl, XTr, YTr)
            if size(YTr, 2) ~= size(XTr, 2)  &&  size(YTr, 1) ~= size(XTr, 2)
                error('Number of columns in YTr, and XTr must be the same (i.e., there must be one and exactly one label for each data point)')
            end

            unique_labels = unique(YTr);
            for i = 1:length(unique_labels)
                template{i} = reshape(mean(XTr(:, (YTr == unique_labels(i)),:), 2,'omitnan'),size(XTr,1),1,[]);
            end
            cl.templates = cell2mat(template);
            cl.labels = unique_labels;

        end

        function [predicted_labels, decision_values] = test(cl, XTe)
            XT = repmat(XTe,[ones(1,length(size(XTe))),length(cl.labels)]);
            tmps = permute(repmat(cl.templates,[ones(1,length(size(XTe))),size(XTe,2)]),[1 length(size(XTe))+1,3:length(size(XTe)),2]);
            template_corrcoeffs = sqrt(sum(XT-tmps,[1:length(size(XT))-2]).^2);
            [~, ind] = min(squeeze(template_corrcoeffs)');
            predicted_labels = cl.labels(ind);
            decision_values = squeeze(template_corrcoeffs)';

        end

    end  % end public methods

end   % end classdef








