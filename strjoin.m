function out = strjoin(C,delimiter)
%% Implementation of strjoin to prevent legacy errors.  Joins a cell array of strings with a delimiter.
    out = cellfun(@(x)[x delimiter],C,'UniformOutput',false);
    out = cat(2,out{:});
    out = out(1:end-numel(delimiter));
end