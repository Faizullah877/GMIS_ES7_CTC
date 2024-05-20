% **********************************************************
% Code to test finding models, creating quadtree and determing residue 
% Functions for testing decoder side steps are provided.
% Decoder side: 
% Add residue to prediction to produce output.
%
%
%
%
% set_name: name of the image set (e.g. 'AC_Boy_basketball')
% ModelID: model ID number (3: rigid, 6: affine, 8: projective)
% MAX_NUM_MODELS: maximum number of models for a single image (1, 2 or 3)
% ref_frm_number: reference frame number (e.g. 1)
% tgt_frm_number: target frame number (e.g. 2)
% QTree_Depth: depth of the quadtree used when encoding (4)
% SMALLEST_BLK_SIZE: smallest block size used when encoding (16)
% **********************************************************
function dec_image = inverse_residue_regis_prediction_v3_based(set_name, residue_folder, ref_recon_image, target_image, dec_resid_name, ModelID, MAX_NUM_MODELS, QTree_Depth, SMALLEST_BLK_SIZE, imgNo)
%    addpath('C:\Users\reji_\MatlabProjects\secondDWT\wt');
    addpath('pfm');

    out_folder = fullfile(residue_folder, "data");
    UP_SHIFT = 1; %1
    pix_max_input = 255;
    pix_max = pix_max_input;
    pix_min = 0;
    % flag 0: to signal rgb color mode (default)
    % flag 1: to signal xyb color mode
    % flag 2: to signal jpegYCbCr color mode
    inFlags.color_conv_flag = 2;
    % flag 0: to allow intra mode
    % flag 1: to replace intra mode with pred of mv=0 (dpcm)
    inFlags.replace_intra_with_dpcm = 0;

    COLOR_CONV_FLAG = inFlags.color_conv_flag;

    ref = ref_recon_image;
    tgt = target_image;

    % model file name
    mdl_dir = fullfile(residue_folder, "model");
    if(~exist(mdl_dir, 'dir')), mkdir(mdl_dir); end

    file_model_name =fullfile(mdl_dir, sprintf("%s_img%03d_transfrom_list_model_id_%d_num_models_%d.mat", set_name,imgNo, ModelID, MAX_NUM_MODELS));
% 
    if isfile(file_model_name)
        inList = load(file_model_name, "tformList");
        tformList = inList.tformList;
        ListSize = size(tformList, 1);
    else
        fprintf('error: unable to read model file containing warping parameters\n');
    end

    % write QuadTree decision and model index to bitsream
%     qtBitStream_name = [out_folder 'qtBitStream.bin'];
    qtBitStream_name = fullfile(out_folder,sprintf("qtBitStream_img%03d.bin",imgNo));

    % residue is saved as uint16
    % out_res_name = [out_folder 'residue.ppm'];

    % Decoder side function calls ****************************************
    % Need to decode residue at various  rates (lower or equal to ref).
    % Currently just copying the residue output (for testing).
    % An encoder and decoder should be called for the residue as noted
    % above and the result saved as 'dec_residue.ppm'.
%     dec_resid_name = [out_folder 'dec_residue.ppm'];

    % for testing copy the residue file to dec_residue.ppm
    % copyfile(out_res_name, dec_resid_name);
    
    % Decode quadtree and read models and perform warping to determine
    % warped prediction based on models and quadtree.
    % Warping is applied to ref (reference image) input 
    [qt_bits, pred, I_Blk_Map] = dec_pred_warp(pix_min, pix_max, ref, size(tgt), tformList, qtBitStream_name, QTree_Depth, SMALLEST_BLK_SIZE, inFlags);

    % read decoded residue
    % Note that residue is saved as uint16 which is then converted to
    % double here

    if COLOR_CONV_FLAG == 1
        pfm_folder = 'data\';
        out_dec_name = [pfm_folder 'dec_residue.pfm'];
        dec_residue = read_pfm(out_dec_name, 0);
    else
        dec_residue = double(imread(dec_resid_name));
    end

    % add warped ref (i.e. prediction) to residue
    [dec_tgt] = dec_add_pred_res(pix_min, pix_max, pred, dec_residue, I_Blk_Map, inFlags, UP_SHIFT);

    if COLOR_CONV_FLAG == 1
        pfm_folder = 'data\';
        out_dec_name = [pfm_folder 'dec.pfm'];
        write_pfm(dec_tgt, out_dec_name, 'l', 0);
        
        dec_tgt = conv_xyb_to_srgb(double(dec_tgt)) * 2^UP_SHIFT;
    elseif COLOR_CONV_FLAG == 2
        dec_tgt = uint8(round(dec_tgt/2^UP_SHIFT));
        tgt_ycbcr = conv_rgb2jpegycbcr(tgt);

        mse = immse(tgt_ycbcr(:,:,1),dec_tgt(:,:, 1));

        dec_tgt = conv_jpegycbcr2rgb(dec_tgt);
        dec_tgt = double(dec_tgt) * 2^UP_SHIFT;
    else

        % print MSE and bits needed for quadtree and motion models
        tgt_ycbcr = rgb2ycbcr(double(tgt));
        dec_tgt_ycbcr = rgb2ycbcr(dec_tgt/2^UP_SHIFT);
    
        mse = immse(tgt_ycbcr(:,:,1),dec_tgt_ycbcr(:,:, 1));
    end
    
    % Bits for quadtree determined from reading the quadtree bitstream.
    % For model assuming max of 3 models with 8 parameters with each
    % parameter requiring 32 bit representation
%     fprintf("\t\tqt bits %d, max_model_bits %d, mse %f\n", qt_bits, (32*8*3), mse);

    dec_image = uint8(round(dec_tgt/2^UP_SHIFT));
%     figure; imagesc(dec_image);
%     title('dec image');
%     imwrite(dec_image, [out_folder 'dec.png']);
end

%% perform a smooth roll off at warpped region boundaries
function [smth] = smooth_roll_off(ref_w, in_msk_w)
    out_w = double(ref_w) .* double(in_msk_w);
    sum_w = double(in_msk_w);
    
    out_w = imgaussfilt(out_w, 1);
    sum_w = imgaussfilt(sum_w, 1);

    new_mask_w = double(sum_w >= 0.9999);

    smth = zeros(size(ref_w));
    for c = 1:1:size(ref_w, 3)
        smth(:,:,c) = (new_mask_w(:,:) .* double(ref_w(:,:,c))) + ((1-new_mask_w(:,:)) .* out_w(:,:,c));
    end
end

%% decode residue at various rates and measure psnr
function [bpp_list, psnr_list] = dec_inter_rd(KDU_DIR, dec_layer_min, dec_layer_max, pix_min, pix_max, ref, tgt, UP_SHIFT, tformList, j2k_out, dec_name, qtBitStream_name, QTree_Depth, SMALLEST_BLK_SIZE)
    
    % dec quadtree and read models and perform warping to determine
    % warped prediction
    [qt_bits, pred, ~] = dec_pred_warp(pix_min, pix_max, ref, size(tgt), tformList, qtBitStream_name, QTree_Depth, SMALLEST_BLK_SIZE, UP_SHIFT);
    
    bpp_list = zeros(1, (dec_layer_max-dec_layer_min+1));
    psnr_list = zeros(1, (dec_layer_max-dec_layer_min+1));

    % decode residue
    for layerCnt = dec_layer_min:dec_layer_max
        DEC_COMMAND = [KDU_DIR 'kdu_expand.exe -i ' j2k_out ' -num_threads 0 -precise -layers ' num2str(layerCnt) ' -o ' dec_name];
        [status_dec,cmdout_dec] = system(DEC_COMMAND);
    
        line_start_index = strfind(cmdout_dec,'(excluding any file format) =');
        bpp = sscanf(cmdout_dec(line_start_index:end),'(excluding any file format) = %f');
        
        dec_residue = double(imread(dec_name));

        % add warped ref to residue
        [dec_tgt] = dec_add_pred_res(pix_min, pix_max, pred, dec_residue, UP_SHIFT);
        
        tgt_ycbcr = rgb2ycbcr(double(tgt));
        dec_tgt_ycbcr = rgb2ycbcr(dec_tgt/2^UP_SHIFT);

        mse = immse(tgt_ycbcr(:,:,1),dec_tgt_ycbcr(:,:, 1));
        psnr = 10 * log10((pix_max*pix_max)/mse);
        % bits for quad-tree, 8 parameter projective models x 3
        metadata_bpp = double(qt_bits + (8*3*32))/(size(ref, 2) * size(ref, 1));
        bpp_list(layerCnt) = bpp + metadata_bpp;
        psnr_list(layerCnt) = psnr;
    end
end

%% warp decoded reference to produce prediction at the decoder
function [qt_bits, pred, I_Blk_Map] = dec_pred_warp(pix_min, pix_max, ref, tgtSize, tformList, qtBitStream_name, QTree_Depth, SMALLEST_BLK_SIZE, inFlags)
    
    REPLACE_INTRA_DPCM = inFlags.replace_intra_with_dpcm;
    COLOR_CONV_FLAG = inFlags.color_conv_flag;
    DWT_METRIC_FLAG = 0;

    [Dec_Index_Map, qt_bits] = readQT_fromBitstream(QTree_Depth, SMALLEST_BLK_SIZE, tgtSize(2), tgtSize(1), qtBitStream_name);
    
    % create index map of block_size index decisions 
    I_Blk_Map = sum(Dec_Index_Map, 3);

    % determine pred
    pred = zeros(tgtSize);

     if COLOR_CONV_FLAG == 1
        ref = conv_srgb_to_xyb(double(ref));
     else
        ref = double(ref); 
     end

    msk = double(ones(size(ref,1), size(ref,2)));
    msk2 = upSample2(msk, 0, 1);

    ref2 = upSample2(ref, pix_min, pix_max);

    for i = 1:size(tformList, 1)
        msk_mdl(:,:,1) = double(I_Blk_Map == i+1);
        msk_mdl(:,:,2) = double(I_Blk_Map == i+1);
        msk_mdl(:,:,3) = double(I_Blk_Map == i+1);
    
        tform2 = tformList(i); %get tform;
        tform2.A(1:2, 3) = tform2.A(1:2, 3) * 2; 
        tform2.A(3, 1:2) = tform2.A(3, 1:2) * 0.5;
        
        ref2_w = imwarp(ref2,tform2,'OutputView',imref2d(tgtSize .* [2 2 1])); 
        % smoothing for each prediction valid region. can turn on/off
        msk2_w = imwarp(msk2,tform2,'OutputView',imref2d(tgtSize .* [2 2 1])); 
        ref2_w = smooth_roll_off(ref2_w, msk2_w);

        ref_w = downSample2(ref2_w, pix_max, pix_min);

        if (DWT_METRIC_FLAG == 1)
            A = ref_w;
            for c=1:size(ref_w,3)
                for lvl=0:3
                    spacing = 2^lvl;
                    offset = 0;
                    A(:,:,c) = analysis_9_7(analysis_9_7(A(:,:,c),spacing,offset)',spacing,offset)';
                end
            end
            ref_w = A;
        end

        % add warped ref to residue
        pred = pred + (ref_w .* msk_mdl);
    end

    if (DWT_METRIC_FLAG == 1)
        A = pred;
        for c=1:size(ref_w,3)
            for lvl=3:-1:0
                spacing = 2^lvl;
                offset = 0;
                A(:,:,c) = synthesis_9_7(synthesis_9_7(A(:,:,c),spacing,offset)',spacing,offset)';
            end
        end
        pred = A;
    end
    pred(pred<pix_min)=pix_min;
    pred(pred>pix_max)=pix_max;

    % for testing - instead of intra (zero prediction) try inter prediction
    % with corresponding reference image region (i.e. mv=0)
    if (DWT_METRIC_FLAG == 0) && (REPLACE_INTRA_DPCM == 1)
        msk_mdl(:,:,1) = double(I_Blk_Map == 1);
        msk_mdl(:,:,2) = double(I_Blk_Map == 1);
        msk_mdl(:,:,3) = double(I_Blk_Map == 1);

        pred = pred + (ref .* msk_mdl);
    end
end

%% add residue to prediction
function [dec_tgt] = dec_add_pred_res(pix_min, pix_max, pred, dec_residue, I_Blk_Map, inFlags, UP_SHIFT)
    REPLACE_INTRA_DPCM = inFlags.replace_intra_with_dpcm;
    COLOR_CONV_FLAG = inFlags.color_conv_flag;

    if COLOR_CONV_FLAG == 1
        % add xyb in domain
        dec_tgt = (pred + dec_residue);
    else

        if COLOR_CONV_FLAG == 2
            pred = conv_rgb2jpegycbcr(uint8(pred));
        end
        pix_max = pix_max * 2^UP_SHIFT;
        pred = double(pred) * 2^UP_SHIFT;

        offset = floor((pix_max+1)/2);

        if (REPLACE_INTRA_DPCM == 1)
            offset_map = repmat( ((I_Blk_Map > 0) * offset), 1, 1, 3 );
        else
            offset_map = repmat( ((I_Blk_Map > 1) * offset), 1, 1, 3 );
        end
    
        dec_tgt = (pred + dec_residue - offset_map - offset);
        dec_tgt(dec_tgt < pix_min) = pix_min;
        dec_tgt(dec_tgt > pix_max) = pix_max;
    end
end

%% decode intra frame at various rates and measure psnr
function [bpp_list, psnr_list] = decode_intra_rd(KDU_DIR, UP_SHIFT, dec_layer_min, dec_layer_max, pix_min, pix_max, tgt, j2k_out, dec_name)
    bpp_list = zeros(1, (dec_layer_max-dec_layer_min+1));
    psnr_list = zeros(1, (dec_layer_max-dec_layer_min+1));

    % tgt = tgt/2^UP_SHIFT;
    % pix_max = pix_max/2^UP_SHIFT;
    % pix_min = pix_min/2^UP_SHIFT;

    for layerCnt = dec_layer_min:dec_layer_max
        DEC_COMMAND = [KDU_DIR 'kdu_expand.exe -i ' j2k_out ' -num_threads 0 -precise -layers ' num2str(layerCnt) ' -o ' dec_name];
        [status_dec,cmdout_dec] = system(DEC_COMMAND);
    
        line_start_index = strfind(cmdout_dec,'(excluding any file format) =');
        bpp = sscanf(cmdout_dec(line_start_index:end),'(excluding any file format) = %f');
    
        dec_tgt = double(imread(dec_name));
    
        % add warped ref to residue
        dec_tgt(dec_tgt < pix_min) = pix_min;
        dec_tgt(dec_tgt > pix_max) = pix_max;
    
        tgt_ycbcr = rgb2ycbcr(double(tgt));
        dec_tgt_ycbcr = rgb2ycbcr(dec_tgt);

        mse = immse(tgt_ycbcr(:,:,1), dec_tgt_ycbcr(:,:,1));
        psnr = 10 * log10((pix_max*pix_max)/mse);
        bpp_list(layerCnt) = bpp;
        psnr_list(layerCnt) = psnr;
    end
end

%% upsample by 2 with sinc iterp
function up2 = upSample2(in_ref, pix_min, pix_max)
    % define filter at x2 resolution (half pixel interpolation)
    % filter defined of +/-4 at full precision
    N = 4;
    x = double([1:2*N]);

    intr = zeros(1,(N*4)+1);

    intr(1, (2*N)+2:(N*4)+1) = (sin(pi * x/2) ./ (pi * x/2)) .* ((1+cos(pi * x/9))/2);
    intr(1, 1:2*N) = intr(1, (N*4)+1:-1:(2*N)+2);
    intr(1, (2*N)+1) = 0;
    DCgain = sum(intr(:));
    intr = intr/DCgain;

    up2 = zeros(size(in_ref,1)*2, size(in_ref,2)*2, size(in_ref,3));

    for i = 1:size(in_ref, 3)
        
        ref = in_ref(:,:, i);

        refExt = zeros((size(ref,1)+(2*N))*2, (size(ref,2)+(2*N))*2);
        refExt((2*N)+1:2:end-(2*N)-1, (2*N)+1:2:end-(2*N)-1) = ref;
    
        refExt(1:1:(2*N), 1:1:end) = refExt((2*N)+1+(2*N):-1:(2*N)+1+1, 1:1:end);
        refExt(end:-1:end-(2*N), 1:1:end) = refExt(end-(2*N)-1-(2*N)-1:1:end-(2*N)-1-1, 1:1:end);
        refExt(1:1:end, 1:1:(2*N)) = refExt(1:1:end, (2*N)+1+(2*N):-1:(2*N)+1+1);
        refExt(1:1:end, end:-1:end-(2*N)) = refExt(1:1:end, end-(2*N)-1-(2*N)-1:1:end-(2*N)-1-1);
    
        % perform convolution - three interpolation phases
        hrz = conv2(refExt, intr, 'same');
    
        vrt = conv2(refExt, intr', 'same');
    
        hExt = hrz;
        hExt(1:1:(2*N), 1:1:end) = hExt((2*N)+1+(2*N):-1:(2*N)+1+1, 1:1:end);
        hExt(end:-1:end-(2*N), 1:1:end) = hExt(end-(2*N)-1-(2*N)-1:1:end-(2*N)-1-1, 1:1:end);
        hExt(1:1:end, 1:1:(2*N)) = hExt(1:1:end, (2*N)+1+(2*N):-1:(2*N)+1+1);
        hExt(1:1:end, end:-1:end-(2*N)) = hExt(1:1:end, end-(2*N)-1-(2*N)-1:1:end-(2*N)-1-1);
    
        ctr = conv2(hExt, intr', 'same');
    
        tmp = hrz + vrt + ctr + refExt;
    
        up2(:,:,i) = tmp((2*N)+1:1:end-(2*N), (2*N)+1:1:end-(2*N));
        up2(:,:,i) = min(pix_max, up2(:,:,i));
        up2(:,:,i) = max(pix_min, up2(:,:,i));
%         up2(up2(:,:,i)>pix_max) = pix_max;
%         up2(up2(:,:,i)<pix_min) = pix_min;
    end
    %figure; imagesc(out2);
end

%% sinc filter and down sample by 2
function down2 = downSample2(in_ref, pix_max, pix_min)
    % define filter at input resolution 
    % filter defined of +/-8 at full precision
    N = 8;
    x = double([1:N]);

    intr = zeros(1,(N*2)+1);

    intr(1, N+2:(N*2)+1) = (sin(pi * x/2) ./ (pi * x/2)) .* ((1+cos(pi * x/(N+1)))/2);
    intr(1, 1:N) = intr(1, (N*2)+1:-1:N+2);
    intr(1, N+1) = 1;
    DCgain = sum(intr(:));
    intr = intr/DCgain;

    down2 = zeros(size(in_ref,1)/2, size(in_ref,2)/2, size(in_ref,3));

    for i = 1:size(in_ref, 3)

        ref = in_ref(:,:, i);

        refExt = zeros((size(ref,1)+(2*N)), (size(ref,2)+(2*N)));
        refExt(N+1:1:end-N, N+1:1:end-N) = ref;
    
        refExt(1:1:N, 1:1:end) = refExt((2*N)+1:-1:N+2, 1:1:end);
        refExt(end:-1:end-N+1, 1:1:end) = refExt(end-(2*N):1:end-N-1, 1:1:end);
        refExt(1:1:end, 1:1:N) = refExt(1:1:end, (2*N)+1:-1:N+2);
        refExt(1:1:end, end:-1:end-N+1) = refExt(1:1:end, end-(2*N):1:end-N-1);
    
        % perform convolution - three interpolation phases
        hrz = conv2(refExt, intr, 'same');
    
        vrt = conv2(hrz, intr', 'same');
    
        down2(:,:,i) = vrt(N+1:2:end-N, N+1:2:end-N);
        down2(:,:,i) = min(pix_max, down2(:,:,i));
        down2(:,:,i) = max(pix_min, down2(:,:,i));
        %down2(down2(:,:,i)>pix_max) = pix_max;
        %down2(down2(:,:,i)<pix_min) = pix_min;
        %figure; imagesc(out2);
    end
end

%% Computes the probability density function of the 2d Gaussian
function p = multivariateGaussian(X, mu, Sigma2)

    p = zeros(size(X, 1), 1);
    k = length(mu);
    
    X = X - mu;

    N = 2000;
    part = floor(size(X, 1)/N);

    for i = 1:N
        
        start_indx = (part * (i-1))+1;

        end_indx = start_indx + part - 1;
        if end_indx > size(X, 1)
            end_indx = size(X, 1);
        end

        M = sum((X(start_indx:end_indx, :)/Sigma2) .* X(start_indx:end_indx, :), 2);
        
        p(start_indx:end_indx, :) = (2 * pi) ^ (- k / 2) * det(Sigma2) ^ (-0.5) * exp(-0.5 * M);
    end

end

% read QT from stored bitstream
function [Index_Map, qt_bits] = readQT_fromBitstream(QTree_Depth, SMALLEST_BLK_SIZE, imgW, imgH, streamFileName)

    byteStream = byteArrayStream(900);

    % write out bytes
    qtFile = fopen(streamFileName);
    inBits = uint8(fread(qtFile));
    fclose(qtFile);

    byteStream.byteArray(1:size(inBits)) = uint8(inBits);

    % for each level: from coarsest (QTree_Depth) to finest (1)
    Leaf_Map = zeros(imgH, imgW);
    Index_Map = zeros(imgH, imgW, QTree_Depth);

    for QT_Level = QTree_Depth:-1:1
    
        blk_size = SMALLEST_BLK_SIZE * 2^(QT_Level-1);
        
        for row_start = 1:blk_size:size(Index_Map, 1)
            for col_start = 1:blk_size:size(Index_Map, 2)
    
                row_end = row_start + blk_size - 1;
                if row_end > imgH
                    row_end = imgH;
                end
        
                col_end = col_start + blk_size - 1;
                if col_end > imgW
                    col_end = imgW;
                end
    
                % if not pruned away then need to read signal
                if Leaf_Map(row_start, col_start) == 0
                    % branch if zero else leaf
                    if QT_Level > 1
                        [byteStream, outVal] = byteStream.readBits(1);
                    else
                        outVal = 1;
                    end

                    if (outVal == 0) && (QT_Level > 1) 
                        Index_Map(row_start:row_end, col_start:col_end, QT_Level) = 0;
                    else
                        [byteStream, outVal] = byteStream.readBits(2);
    
                        % write model ID - 4 options so two bits
                        Index_Map(row_start:row_end, col_start:col_end, QT_Level) = outVal + 1;
                    end
                end
            end
        end
            
        Leaf_Map = Leaf_Map + double(Index_Map(:,:,QT_Level) > 0);
        Leaf_Map = double(Leaf_Map(:,:) > 0);
    
%         figure; imagesc(Index_Map(:,:,QT_Level));
%         title(['Quadtree Level ', num2str(QT_Level)]);
    end
    qt_bits = byteStream.byteCnt * 8;
end