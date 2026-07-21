classdef cosine_CL
    % cosine_CL is a classifier object (CL)
    % Fits an L2-regularized multinomial logistic regression model.
    %
    % INPUTS
    %   XTr          : trials x features matrix
    %   YTr          : integer labels from 1 to numClasses
    %
    % OUTPUT
    %   model.templates
    %   model.labels
    %
    % Requires the Optimization Toolbox for fminunc
    %==========================================================================
    %==========================================================================
    properties
        templates = [];   % average of the training vectors for each class
        labels = [];  % all the unique labels for each class (one for each template)
    end

    methods
        % constructor
        function cl = cosine_CL
        end

        function cl = train(cl,XTr, YTr)
            YTr = YTr(:);
            cl.labels = unique(YTr);
            temps = cell(1,length(cl.labels));
            for c =1:length(cl.labels)
                temps{c} = permute(mean(XTr(:, YTr == cl.labels(c), :), 2, 'omitnan'),[1 3 2]);
            end
            cl.templates = cell2mat(reshape(temps,1,1,[]));
        end


        function [predicted_labels, decision_values] = test(cl, X)
            decision_values = zeros(size(X,2),length(cl.labels));
            templateMatrix = reshape(cl.templates, [], length(cl.labels));
            for t = 1:size(X,2)
                testVector = reshape(X(:,t,:),[],1);
                validFeatureMask = isfinite(testVector) & all(isfinite(templateMatrix), 2);
                testVector = testVector(validFeatureMask);
                for c = 1:length(cl.labels)
                    decision_values(t,c) = sum((sum(testVector.*templateMatrix(validFeatureMask,c)) ./...
                        sqrt(sum(testVector.^2) .* sum(templateMatrix(validFeatureMask,c).^2)))/size(X,1));
                end
            end
            [~, ind] = randmax(decision_values');  
            predicted_labels = cl.labels(ind);
        end

    end  % end public methods
end   % end classdef