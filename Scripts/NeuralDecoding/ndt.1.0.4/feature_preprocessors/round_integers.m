classdef round_integers
     properties 
        avgData = [];
     end

    methods
        function fp = round_integers
        end
        function current_FP_info_to_save = get_current_info_to_save(fp)
            current_FP_info_to_save = []; 
        end
        function [fp,X_rounded] = set_properties_with_training_data(fp, XTr,tilda)
            fp.avgData = XTr;
            X_rounded = fp.preprocess_test_data(XTr);
        end
        function X_rounded = preprocess_test_data(fp, X_data)
            X_rounded = round(X_data);
        end
    end
end