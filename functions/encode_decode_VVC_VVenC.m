function report = encode_decode_VVC_VVenC(input)
report = struct();
set_info = input.set_info;
set_name = input.set_name;
no_of_images = set_info.no_of_images;
q_v = input.q_value;
no_of_threads = 8;

encoder_exe = fullfile(input.codec_folder, "VVC_VVenC", "vvencFFapp.exe");
decoder_exe = fullfile(input.codec_folder, "VVC_VVenC/vvdecapp.exe");
ffmpeg_exe = fullfile(input.codec_folder, "ffmpeg.exe");
rnd_cfg_file = fullfile(input.codec_folder, "VVC_VVenC/cfg/randomaccess_slow.cfg");
seq_cfg_file = fullfile(input.codec_folder, "VVC_VVenC/cfg/sequence.cfg");

ImagesResult = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};

if strcmp(input.config, "intra")
    for imgNo = 1 : no_of_images
        no_rows = set_info.images_details.ImageHeight(imgNo);
        no_cols = set_info.images_details.ImageWidth(imgNo);
        image_name = set_info.images_details.NameWoExt(imgNo);
        yuv_file_name = sprintf("%s_%dx%d_1Hz_prec10_P420.yuv", image_name, no_cols, no_rows);
        input_yuv_image = fullfile(input.input_raw, yuv_file_name); % file naming format is followed in conversion.
        if ~exist(input_yuv_image, 'file')
            disp(input_yuv_image);
            error("Input File does not exists");
        end
        bin_image_file = fullfile(input.compressed_folder, sprintf("%s.bin", image_name));
        enc_cmd = sprintf("%s --preset=slow -c %s -c %s -q %d --InputFile %s -s %dx%d --InputBitDepth 10 --FramesToBeEncoded %d --BitstreamFile %s --InputChromaFormat 420 --framerate 30 --IntraPeriod=1  -t %d", encoder_exe, rnd_cfg_file, seq_cfg_file, q_v, input_yuv_image, no_cols, no_rows, no_of_images, bin_image_file, no_of_threads);
        tCstartimg= tic;
        [status, outE] = system(enc_cmd);
        tCcomplete = toc(tCstartimg);
        img_enc_time = tCcomplete;
        if(status)
            fprintf("VVC VVenC intra => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_image);
            disp(outE);
            error('compression failed');
        end
        decoded_image_yuv = fullfile(input.decompressed_folder, sprintf("Rec_%s.yuv", image_name));
        dec_cmd = sprintf("%s -b %s -o %s ", decoder_exe,  bin_image_file, decoded_image_yuv);
        tDstartimg= tic;
        [status, outE] = system(dec_cmd);
        img_dec_time = toc(tDstartimg);
        if(status)
            disp(outE);
            fprintf("VVC vvdecapp => SetName: %s, QP: %d, => Failed to decompress [%s].\n", set_name,  q_v, bin_image_file);
            error('compression failed');
        end
        compressed_file_info=dir(bin_image_file);
        compressed_img_size = compressed_file_info.bytes;
        img_bpp = (compressed_img_size*8) / (double(no_cols)*double(no_rows));
        ppm_output=fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name));
        conv_cmd = sprintf('%s -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', ffmpeg_exe, no_cols, no_rows, decoded_image_yuv, ppm_output);
        [status, out] = system(conv_cmd);
        if(status)
            disp(out);
            fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_image_yuv, ppm_output);
            error('conversion failed');
        end
        delete(decoded_image_yuv);
        ImagesResult(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
    end
    report.set_enc_time = sum(ImagesResult.("enc_time"));
    compressed_set_size = sum(ImagesResult.("CompressedImageSize"));
    report.compressed_set_size = compressed_set_size;
    report.set_bpp = (compressed_set_size*8) / double(set_info.num_pixels);
    report.set_dec_time = sum(ImagesResult.("dec_time"));
    report.set_dec_folder = input.decompressed_folder;
    report.images_result_table = ImagesResult;

else
    if ~set_info.same_resolution
        fprintf("The codec {VVC-VVenC-Inter} does not support sets {%s} with images having different resolution.", set_name);
        error("Compression Failed");
    end
    no_rows = set_info.images_details.ImageHeight(1);
    no_cols = set_info.images_details.ImageWidth(1);
    input_yuv_video = fullfile(input.input_raw, sprintf("%s_%dx%dx%d_1Hz_prec10_P420.yuv", set_name, no_cols,no_rows, no_of_images)); % file naming format is followed in conversion.
    bin_video_file = fullfile(input.compressed_folder, sprintf("%s.266", set_name));     
    enc_cmd = sprintf("%s --preset=slow -c %s -c %s -q %d --InputFile %s -s %dx%d --InputBitDepth 10 --FramesToBeEncoded %d --BitstreamFile %s --InputChromaFormat 420 --MinSearchWindow 1 --framerate 30 -t %d", encoder_exe, rnd_cfg_file, seq_cfg_file, q_v, input_yuv_video, no_cols, no_rows, no_of_images, bin_video_file, no_of_threads);
    tCstartimg= tic;
    [status, outE] = system(enc_cmd);
    report.set_enc_time = toc(tCstartimg);
    if(status)
        disp(outE);
        fprintf("VVC VVC_VVenC => SetName: %s, QP: %d, => Failed to compress [%s].\n", set_name,  q_v, input_yuv_video);
        error('compression failed');
    end

    compressed_file_info=dir(bin_video_file);
    report.compressed_file = bin_video_file;
    report.compressed_set_size = compressed_file_info.bytes;
    report.set_bpp = (compressed_file_info.bytes*8) / (double(no_cols)*double(no_rows)*double(no_of_images));   
    num_of_bits_Table = get_images_bpp_vvc_out(outE);
    % time_taken = extract_enc_time(outE);
    % report.set_enc_time = time_taken;
    % decoding
    decoded_video_yuv = fullfile(input.decompressed_folder, sprintf("Rec_%s.yuv", set_name));
    dec_cmd = sprintf("%s -b %s -o %s ", decoder_exe,  bin_video_file, decoded_video_yuv);
    tDstartimg= tic;
    [status, out] = system(dec_cmd);
    report.set_dec_time = toc(tDstartimg);
    if(status)
        disp(out);
        fprintf("VVC vvdecapp => SetName: %s, QP: %d, => Failed to decompress [%s].\n", set_name,  q_v, bin_video_file);
        error('compression failed');
    end

    % tConStartT = tic;
    ppm_output=fullfile(input.decompressed_folder,"dec_image_%03d.ppm");
    conv_cmd = sprintf('%s -f rawvideo -vcodec rawvideo -s %dx%d -pix_fmt yuv420p10le -i %s -pix_fmt rgb24 -vf scale=in_range=full:in_color_matrix=bt709:out_range=full:out_color_matrix=bt709 -color_primaries bt709 -color_trc bt709 -colorspace bt709 -y %s', ffmpeg_exe, no_cols, no_rows, decoded_video_yuv, ppm_output);
    [status, out] = system(conv_cmd);
    if(status)
        disp(out);
        fprintf("Failed to convert [%s] with [ffmpeg] to [%s].\n", decoded_video_yuv, decoded_rgb_files);
        error('conversion failed');
    end
    delete(decoded_video_yuv);
    report.set_dec_folder = input.decompressed_folder;
    for imgNo = 1: no_of_images
        image_name = set_info.images_details.NameWoExt(imgNo);
        row = num_of_bits_Table.Image_No == imgNo;
        num_of_bits = num_of_bits_Table.Num_of_Bits(row);
        num_of_pixels = no_rows*no_cols;
        compressed_img_size = num_of_bits/8; % to convert to bytes from bits
        img_bpp = num_of_bits/num_of_pixels;
        img_enc_time = report.set_enc_time/no_of_images;
        img_dec_time = report.set_dec_time/no_of_images;
        movefile(fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm",imgNo)), fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name)));
        ImagesResult(imgNo,:) ={imgNo, image_name,  uint16(no_cols), uint16(no_rows),  compressed_img_size, img_bpp, img_enc_time, img_dec_time};
    end
end
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult;
end
function time = extract_enc_time(charVector)
% Split the char vector into lines
lines = strsplit(charVector, '\n');


for i = 1:length(lines)
    line = strtrim(lines{i});  % Remove leading and trailing white spaces
    if startsWith(line, 'vvencFFapp [info]: Total Time:')  % Check if the line starts with 'Total Time:'
        tokens = strsplit(line);  % Split the line into words
        time = str2double(tokens{5});  % Add the number after 'Total Time:' to time
    end
end
end

