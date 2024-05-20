function result= get_images_bpp_vvc_out(charVector)
    % Split the char vector into lines
    lines = strsplit(charVector, '\n');  

    image_no = [];  % Initialize the image_no vector
    num_of_bits = [];  % Initialize the num_of_bits vector

    for i = 1:length(lines)
        line = lines{i};
        if startsWith(line, 'POC')  % Check if the line starts with 'POC'
            tokens = strsplit(line);  % Split the line into words
            image_no = [image_no, str2double(tokens{2})];  % Add the number after 'POC' to image_no
            bitsIndex = find(strcmp(tokens, 'bits'));  % Find the index of 'bits'
            num_of_bits = [num_of_bits, str2double(tokens{bitsIndex-1})];  % Add the number before 'bits' to num_of_bits
        end
    end

    image_no = image_no + 1;
    % Create a table with two columns
    result = table(image_no.', num_of_bits.', 'VariableNames', {'Image_No', 'Num_of_Bits'});

    % Sort the table in ascending order based on the Image_No column
    result = sortrows(result, 'Image_No');
end

