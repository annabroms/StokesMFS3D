function estimate = bubble_svd_memory(rows,columns)
%BUBBLE_SVD_MEMORY Budget dense assembly, SVD factors, copies, and workspace.
% This conservative planning estimate is not an operating-system memory cap.
thin = min(rows,columns);
estimate.rows = rows;
estimate.columns = columns;
estimate.inputGB = 8*rows*columns/1e9;
estimate.knownGB = 8*(2*rows*columns+rows*thin+columns*thin+thin^2)/1e9;
estimate.estimatedGB = max(estimate.knownGB, ...
    8*(8*rows*columns+8*thin^2)/1e9);
end
