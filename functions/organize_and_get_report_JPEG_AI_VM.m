function report = organize_and_get_report_JPEG_AI_VM(input)

report = struct();
no_of_images = input.no_of_images;
ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;
q_v_str = sprintf('%03d', q_v);
no_rows = input.height;
no_cols = input.width;
ai_compressed_folder = fullfile(input.data_folder, "bit");
ai_decompressed_folder = fullfile(input.data_folder, "rec");

raw_files = dir(fullfile(input.input_raw, "*.tiff"));

for imgNo = 1: no_of_images

    input_raw = fullfile(input.input_raw, raw_files(imgNo).name);
    [~, image_name, ~] = fileparts(input_raw);
    compressed_fname = sprintf("VM_%s_%s.bits", image_name, q_v_str);
    source_compressed_path = fullfile(ai_compressed_folder, compressed_fname);
   
    cmp_file_exists = isfile(source_compressed_path);
    if cmp_file_exists
        copyfile(source_compressed_path, input.compressed_folder);
        jpg_image_file_info=dir(source_compressed_path);
        compressed_img_size = jpg_image_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
    else
        compressed_img_size = NaN;
        img_bpp = NaN;
    end

    img_enc_time  = NaN; 

    decompressed_fname = sprintf("VM_%s_%s.png", image_name, q_v_str);
    source_decomp_path = fullfile(ai_decompressed_folder, decompressed_fname);
    dec_file_exists = isfile(source_decomp_path);
    if dec_file_exists
        decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", imgNo));

        ffmpeg_exe = fullfile(input.codec_folder, "ffmpeg.exe");
        command=sprintf('%s -i %s %s',ffmpeg_exe, source_decomp_path, decoded_image_ppm);
        system(command);
    end

    img_dec_time = NaN; %toc(tDstartimg);

    % get bpp and encoding time
    ImagesResult_T1(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end

% replace missing bpp with average
mean_bpp = mean(ImagesResult_T1.bpp, 'omitnan');
nanIndices = isnan(ImagesResult_T1.bpp);
ImagesResult_T1.bpp(nanIndices) = mean_bpp;

% replace misssing filesize with average
mean_size = mean(ImagesResult_T1.CompressedImageSize, 'omitnan');
nanIndices = isnan(ImagesResult_T1.CompressedImageSize);
ImagesResult_T1.CompressedImageSize(nanIndices) = mean_size;

report.set_enc_time = sum(ImagesResult_T1.("enc_time"));
compressed_set_size = sum(ImagesResult_T1.("CompressedImageSize"));
report.compressed_file_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / (input.width*input.height*no_of_images);
report.set_dec_time = sum(ImagesResult_T1.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult_T1;

end


