function function_calls = find_functions_calls(filename)
    % Read the content of the file
    file_content = fileread(filename);

    % Define a regular expression pattern for function definitions
    pattern = 'function\s+(\w+)\s*\(.*?\)';

    % Find all matches of the pattern in the file content
    matches = regexp(file_content, pattern, 'tokens');

    % Extract function names from the matches
    function_calls = unique(cellfun(@(x) x{1}, matches, 'UniformOutput', false));
end