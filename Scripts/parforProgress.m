function queue = parforProgress(totalIterations)
    % Args:
    %     totalIterations (int): Total number of iterations.
    % Returns:
    %     queue (parallel.pool.DataQueue): DataQueue to receive updates.
    % Initialize DataQueue and Progress Bar
    queue = parallel.pool.DataQueue;
    progressBar = waitbar(0, 'Processing...', 'Name', ['Iterating ', num2str(totalIterations), ' instances....']);
    
    % Reset persistent variable count
    persistent count
    count = 0;
    st = tic;
    function updateProgress(~)
        count = count + 1;
        shareComplete = count / totalIterations;
        tIter = toc(st);
        waitbar(shareComplete, progressBar,[num2str(count),...
            ' of ', num2str(totalIterations), ' (', num2str(fix(tIter)), ...
            's elapsed; ~', num2str(fix((totalIterations-count)*(tIter/count))), 's remaining)']);
        if count == totalIterations
            close(progressBar);
            count = [];
            disp(['Processed ', num2str(totalIterations), ' instances in ', num2str(tIter), ' seconds.']);
        end
    end
    % Add listener to the DataQueue
    afterEach(queue, @updateProgress);
end