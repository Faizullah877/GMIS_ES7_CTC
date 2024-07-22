function report = encode_decode_SJU_ARCHs(input)
report = struct();
set_info = input.set_info;
no_of_images = set_info.no_of_images;
set_name = set_info.set_name;
if ~set_info.same_resolution
    fprintf("The codec {gmis.exe} does not support sets {%s} with images having different resolution.", set_name);
    error("Compression Failed");
end

codec = input.codec;
config="";
if strcmp(codec, "JPEG2000")
    encoder_exe = fullfile(input.codec_folder, "JPEG_2000", "kdu_compress.exe");
    decoder_exe = fullfile(input.codec_folder, "JPEG_2000", "kdu_expand.exe");
elseif strcmp(codec, "JPEG_XL")
    encoder_exe = fullfile(input.codec_folder, "JPEG_XL", "cjxl.exe");
    decoder_exe = fullfile(input.codec_folder, "JPEG_XL", "djxl.exe");
elseif strcmp(codec, "JPEG1")
    encoder_exe = fullfile(input.codec_folder, "JPEG1", "cjpeg.exe");
    decoder_exe = fullfile(input.codec_folder, "JPEG1", "djpeg.exe");
    config = input.config;
end

intra_frame = 10;

no_rows = set_info.images_details.ImageHeight(1);
no_cols = set_info.images_details.ImageWidth(1);

brotli_exe = fullfile(input.codec_folder, "JPEG_XL\brotli.exe");

ImagesResult = table('Size', [no_of_images 8], 'VariableTypes', {'uint16', 'string','uint16', 'uint16', 'double', 'single', 'double','double'});
ImagesResult.Properties.VariableNames = {'ImageNo','ImageName', 'ImageWidth', 'ImageHeight', 'CompressedImageSize','bpp', 'enc_time', 'dec_time'};
q_v = input.q_value;

tmp_folder = fullfile(input.set_q_folder, "tmp");
if(~exist(tmp_folder, 'dir')), mkdir(tmp_folder); end

%% Mother Image

image_name_wo_ext = set_info.images_details.NameWoExt(1);
mother_image_file = fullfile(input.input_raw, sprintf("%s.ppm", image_name_wo_ext));
[~, image_name, ~] = fileparts(mother_image_file);
%encoding
[command, mother_image_compressed_file] = get_encoding_command(codec, encoder_exe, mother_image_file, input.compressed_folder, image_name, q_v, config);
tCstartimg= tic;
[status, out] = system(command);
img_enc_time  = toc(tCstartimg);
if(status)
    disp(out);
    fprintf("Failed to compress [%s] with [%s] codec.\n", mother_image_file, codec);
    error('compression failed');
end
% decoding
decoded_image_ppm = fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name_wo_ext));
dec_cmd = get_decoding_command(codec, decoder_exe, mother_image_compressed_file, decoded_image_ppm);
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
img_bpp = (compressed_img_size*8) / (double(no_cols)*double(no_rows));
ImagesResult(1,:) ={1, image_name_wo_ext, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};

% Orig_IMG = imread(mother_image_file);
decoded_mother_img = imread(decoded_image_ppm);
% psnrv = psnr(uint8(decoded_mother_img), Orig_IMG);
% fprintf("%s  => BPP   %f  => PSNR = %f\n", image_name, img_bpp, psnrv);
previous_image = int16(decoded_mother_img);

%% Child Images.
for imgNo = 2: no_of_images
    image_name_wo_ext = set_info.images_details.NameWoExt(imgNo);
    no_rows = set_info.images_details.ImageHeight(imgNo);
    no_cols = set_info.images_details.ImageWidth(imgNo);
    input_raw = fullfile(input.input_raw, sprintf("%s.ppm", image_name_wo_ext));
    Orig_IMG = imread(input_raw);
    image_name = image_name_wo_ext;
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
    [command, compressed_file] = get_encoding_command(codec, encoder_exe, input_file, input.compressed_folder, image_name, q_v, config);
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
        decoded_residue_image_ppm = fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name_wo_ext));
    else
        decoded_residue_image_ppm = fullfile(tmp_folder, sprintf("dec_image_%03d.ppm", imgNo));
    end

    dec_cmd = get_decoding_command(codec, decoder_exe, compressed_file, decoded_residue_image_ppm);
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
        decoded_image_file = fullfile(input.decompressed_folder, sprintf("%s.ppm", image_name_wo_ext));
        imwrite(uint8(rec_next_image), decoded_image_file);
    end

    % psnrv = psnr(uint8(rec_next_image), Orig_IMG);
    previous_image = int16(rec_next_image);

    % get bpp and encoding time
    compressed_file_info=dir(compressed_file);
    compressed_img_size = compressed_file_info.bytes + br1_bytes + br2_bytes;
    img_bpp = (compressed_img_size*8) / (double(no_cols)*double(no_rows));

    % fprintf("%s  => BPP   %f  => PSNR = %f\n", image_name, img_bpp, psnrv);
    ImagesResult(imgNo,:) ={imgNo, image_name, no_cols, no_rows, compressed_img_size, img_bpp, img_enc_time, img_dec_time};
end


report.set_enc_time = sum(ImagesResult.("enc_time"));
compressed_set_size = sum(ImagesResult.("CompressedImageSize"));
report.compressed_set_size = compressed_set_size;
report.set_bpp = (compressed_set_size*8) / double(set_info.num_pixels);
report.set_dec_time = sum(ImagesResult.("dec_time"));
report.set_dec_folder = input.decompressed_folder;
report.images_result_table = ImagesResult;

end

function [command, output_file] = get_encoding_command(codec, encoder_exe, input_raw, compressed_folder, image_name, q_v, config)
if strcmp(codec, "JPEG1")
    if(~exist('config', 'var'))
        config = "arith";
    end
    output_file = fullfile(compressed_folder, sprintf("%s.jpg", image_name));
    if strcmp(config  , "arith")
        command =sprintf('%s -quality %d -arithmetic -outfile %s %s', encoder_exe,  q_v, output_file, input_raw);
    else
        command =sprintf('%s -quality %d -outfile %s %s', encoder_exe,  q_v, output_file, input_raw);
    end
elseif strcmp(codec, "JPEG2000")
    output_file = fullfile(compressed_folder, sprintf("%s.j2c", image_name));
    command = sprintf("%s -i %s -o %s Qfactor=%d", encoder_exe, input_raw, output_file, q_v);
elseif strcmp(codec, "JPEG_XL")
    output_file = fullfile(compressed_folder, sprintf("%s.jxl", image_name));
    command = sprintf('%s %s %s -q %d -e 6 --num_threads 8', encoder_exe, input_raw, output_file, q_v);
end
end

function dec_cmd = get_decoding_command(codec, decoder_exe, compressed_file, decoded_file)
if strcmp(codec, "JPEG1")
    dec_cmd=sprintf('%s -outfile %s %s',decoder_exe, decoded_file, compressed_file);
elseif strcmp(codec, "JPEG2000")
    dec_cmd = sprintf("%s -i %s -o %s", decoder_exe, compressed_file, decoded_file);
elseif strcmp(codec, "JPEG_XL")
    dec_cmd=sprintf('%s %s %s --num_threads 8',decoder_exe, compressed_file, decoded_file);
end
end

