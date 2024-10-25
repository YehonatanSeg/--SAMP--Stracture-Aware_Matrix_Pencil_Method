function [at, bt] = closestNumbers(a, b, array)
    % Compute absolute differences between a and each element in the array
    diff_a = abs(array - a);
    % Compute absolute differences between b and each element in the array
    diff_b = abs(array - b);
    
    % Combine the absolute differences for a and b
    combined_diff = diff_a + diff_b;
    
    % Find the index of the minimum combined difference
    [~, min_index] = min(combined_diff);
    
    % Assign the closest numbers to c and d
    at = array(min_index);
    
    % Remove the element at min_index to find the next closest number
    array(min_index) = [];
    
    % Compute the new absolute differences without the closest number found
    diff_a = abs(array - a);
    diff_b = abs(array - b);
    
    % Combine the new absolute differences
    combined_diff = diff_a + diff_b;
    
    % Find the index of the minimum combined difference in the updated array
    [~, min_index] = min(combined_diff);
    
    % Assign the second closest number to d
    bt = array(min_index);
end