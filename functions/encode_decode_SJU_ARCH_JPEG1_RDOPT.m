function report = encode_decode_SJU_ARCH_JPEG1_RDOPT(input)
report = struct();
no_of_images = input.no_of_images;
no_cols = input.width;
no_rows = input.height;
codec = input.codec;
encoder_exe = fullfile(input.codec_folder, "JPEG1", "cjpeg.exe");
decoder_exe = fullfile(input.codec_folder, "JPEG1", "djpeg.exe");
rdopt_exe = fullfile(input.codec_folder, "JPEG1", "RD_Opt.exe");
config = input.config;
intra_frame = input.intra_refresh_rate;

brotli_exe = fullfile(input.codec_folder, "JPEG_XL\brotli.exe");
raw_files = dir(fullfile(input.input_raw, "*.ppm"));

ImagesResult_T1 = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult_T1.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;

tmp_folder = fullfile(input.set_q_folder, "tmp");
if(~exist(tmp_folder, 'dir')), mkdir(tmp_folder); end

%% Mother Image
mother_image_file = fullfile(input.input_raw, raw_files(1).name);
[~, image_name, ~] = fileparts(mother_image_file);
%encoding

img_qfile_name = fullfile(tmp_folder, sprintf("qtable_image_%03d.txt", 1));
if(exist(img_qfile_name, "file"))
    delete(img_qfile_name);
end
rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, mother_image_file);
[status, out1] = system(rdopt_cmd);
if(status)
    disp(out1);
    fprintf("Failed to generate [%s] with [%s] codec.\n", img_qfile_name, rdopt_exe);
    error('q-tables generation failed');
end
[command, mother_image_compressed_file] = get_encoding_command(encoder_exe, mother_image_file, input.compressed_folder, image_name, config, img_qfile_name);
tCstartimg= tic;
[status, out] = system(command);
img_enc_time  = toc(tCstartimg);
if(status)
    disp(out);
    fprintf("Failed to compress [%s] with [%s] codec.\n", mother_image_file, codec);
    error('compression failed');
end
% decoding
decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", 1));
dec_cmd = get_decoding_command(decoder_exe, mother_image_compressed_file, decoded_image_ppm);
tDstartimg = tic;
[status, out] = system(dec_cmd);
img_dec_time = toc(tDstartimg);
if(status)
    disp(out)
    fprintf("Failed to decompress [%s] with [%s] codec.\n", mother_image_compressed_file, codec);
    error('decompression failed');
end

% get bpp and encoding time
compressed_file_info=dir(mother_image_compressed_file);
compressed_img_size = compressed_file_info.bytes;
img_bpp = (compressed_img_size*8) / (no_cols*no_rows);
ImagesResult_T1(1,:) ={1, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};

Orig_IMG = imread(mother_image_file);
decoded_mother_img = imread(decoded_image_ppm);
% psnrv = psnr(uint8(decoded_mother_img), Orig_IMG);
% fprintf("%s  => BPP   %f  => PSNR = %f\n", image_name, img_bpp, psnrv);
previous_image = int16(decoded_mother_img);

%% Child Images.
for imgNo = 2: no_of_images

    input_raw = fullfile(input.input_raw, raw_files(imgNo).name);
    Orig_IMG = imread(input_raw);
    [no_rows, no_cols, ~] = size(Orig_IMG);
    [~, image_name, ~] = fileparts(input_raw);
    br1_bytes = 0;
    br2_bytes = 0;
    if mod(imgNo, intra_frame) == 0  % intra refresh rate
        input_file = input_raw;
    else
        next_image = int16(Orig_IMG);
        [diff_img, idx1, idx2] = sju_residue_dpcm_pixels_forward(previous_image, next_image);
        idx1_tmp_file = fullfile(tmp_folder, sprintf("%s_idx1.mat", image_name));
        idx2_tmp_file = fullfile(tmp_folder, sprintf("%s_idx2.mat", image_name));
        save(idx1_tmp_file, "idx1");
        save(idx2_tmp_file, "idx2");
        residue_tmp_file = fullfile(tmp_folder, sprintf("residue_%03d.ppm", imgNo));
        imwrite(uint8(diff_img), residue_tmp_file);
        idx1_brotli_file = fullfile(input.compressed_folder, sprintf("%s_idx1.br", image_name));
        idx2_brotli_file = fullfile(input.compressed_folder, sprintf("%s_idx2.br", image_name));
        brotli_cmd = sprintf("%s -f -o %s %s", brotli_exe, idx1_brotli_file, idx1_tmp_file);
        [status1, out1] = system(brotli_cmd);
        if(status1)
            disp(out1);
            error('compression failed');
        end
        brotli_cmd = sprintf("%s -f -o %s %s", brotli_exe, idx2_brotli_file, idx2_tmp_file);
        [status1, out1] = system(brotli_cmd);
        if(status1)
            disp(out1);
            error('compression failed');
        end

        br1_file_info = dir(idx1_brotli_file);
        br1_bytes = br1_file_info.bytes;
        br2_file_info = dir(idx2_brotli_file);
        br2_bytes = br2_file_info.bytes;
        input_file = residue_tmp_file;
    end

    %encoding 

    img_qfile_name = fullfile(tmp_folder, sprintf("qtable_image_%03d.txt", imgNo));
    if(exist(img_qfile_name, "file"))
        delete(img_qfile_name);
    end
    rdopt_cmd = sprintf("%s -numtables 2 -bpp %f -qfile %s -im %s", rdopt_exe, q_v, img_qfile_name, input_file);
    [status, out1] = system(rdopt_cmd);
    if(status)
        disp(out1);
        fprintf("Failed to generate [%s] with [%s] codec.\n", img_qfile_name, rdopt_exe);
        error('q-tables generation failed');
    end
    [command, compressed_file] = get_encoding_command(encoder_exe, input_file, input.compressed_folder, image_name, config, img_qfile_name);
    tCstartimg= tic;
    [status, out] = system(command);
    img_enc_time  = toc(tCstartimg);
    if(status)
        disp(out);
        fprintf("Failed to compress [%s] with [%s] codec.\n", compressed_file, codec);
        error('compression failed');
    end

    % decoding

    if mod(imgNo, intra_frame) == 0 
        decoded_residue_image_ppm = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", imgNo));
    else
        decoded_residue_image_ppm = fullfile(tmp_folder, sprintf("dec_image_%03d.ppm", imgNo));
    end

    dec_cmd = get_decoding_command(decoder_exe, compressed_file, decoded_residue_image_ppm);
    tDstartimg = tic;
    [status, out] = system(dec_cmd);
    if(status)
        disp(out);
        fprintf("Failed to decompress [%s] with [%s] codec.\n", compressed_file, codec);
        error('decompression failed');
    end
    img_dec_time = toc(tDstartimg);


    recon_residue = imread(decoded_residue_image_ppm);

    if mod(imgNo, intra_frame) == 0 
        rec_next_image = recon_residue;
    else
        rec_idx1_tmp_file = fullfile(tmp_folder, sprintf("rec_%s_idx1.mat", image_name));
        rec_idx2_tmp_file = fullfile(tmp_folder, sprintf("rec_%s_idx2.mat", image_name));
        brotli_cmd = sprintf("%s -d -f -o %s %s", brotli_exe, rec_idx1_tmp_file, idx1_brotli_file);
        [status1, out1] = system(brotli_cmd);
        if(status1)
            disp(out1);
            error('compression failed');
        end
        brotli_cmd = sprintf("%s -d -f -o %s %s", brotli_exe, rec_idx2_tmp_file, idx2_brotli_file);
        [status1, out1] = system(brotli_cmd);
        if(status1)
            disp(out1);
            error('compression failed');
        end
        rec_idx1 = load(rec_idx1_tmp_file);
        rec_idx1 = rec_idx1.idx1;
        rec_idx2 = load(rec_idx2_tmp_file);
        rec_idx2 = rec_idx2.idx2;
        rec_next_image = sju_residue_dpcm_pixels_inverse(previous_image, int16(recon_residue), rec_idx1, rec_idx2);
        decoded_image_file = fullfile(input.decompressed_folder, sprintf("dec_image_%03d.ppm", imgNo));
        imwrite(uint8(rec_next_image), decoded_image_file);
    end

    % psnrv = psnr(uint8(rec_next_image), Orig_IMG);
    previous_image = int16(rec_next_image);

    % get bpp and encoding time
    compressed_file_info=dir(compressed_file);
    compressed_img_size = compressed_file_info.bytes + br1_bytes + br2_bytes;
    img_bpp = (compressed_img_size*8) / (no_cols*no_rows);

    % fprintf("%s  => BPP   %f  => PSNR = %f\n", image_name, img_bpp, psnrv);
    ImagesResult_T1(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end


report.set_enc_time = sum(ImagesResult_T1.("enc_time"));
compressed_set_size = sum(ImagesResult_T1.("CompressedImageSize"));
report.compressed_file_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / (input.width*input.height*no_of_images);
report.set_dec_time = sum(ImagesResult_T1.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult_T1;

end

function [command, output_file] = get_encoding_command(encoder_exe, input_raw, compressed_folder, image_name, config, q_table_file)
    if(~exist('config', 'var'))
        config = "arith";
    end
    output_file = fullfile(compressed_folder, sprintf("%s.jpg", image_name));
    if strcmp(config  , "arith")
        command =sprintf('%s -qtables %s -arithmetic -outfile %s %s', encoder_exe,  q_table_file, output_file, input_raw);
    else
        command =sprintf('%s -qtables %s -outfile %s %s', encoder_exe,  q_table_file, output_file, input_raw);
    end
end

function dec_cmd = get_decoding_command(decoder_exe, compressed_file, decoded_file)
    dec_cmd=sprintf('%s -outfile %s %s',decoder_exe, decoded_file, compressed_file);
end

