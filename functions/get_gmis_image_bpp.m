function img_bpp = get_gmis_image_bpp(out1, imgNo)
% Assuming 'text' is your multiline text
lines = splitlines(out1);
img_bpp = NaN;
for i = 1:length(lines)


    line = lines{i};
    if ~contains(line, "ImageNo")
        continue;
    end
    words = split(line);
    if strcmp(words{2}, "ImageNo")
        if str2double(words{4}) == double(imgNo)
            img_bpp = str2double(words{end});
            break;
        end
    end
end

end