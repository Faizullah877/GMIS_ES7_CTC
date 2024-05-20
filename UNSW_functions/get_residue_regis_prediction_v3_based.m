% **********************************************************
% Code to test finding models, creating quadtree and determing residue 
% Functions for testing encoder and decoder side steps are provided.
% Encoder side:
% Perform warping and prediction to produce a residue signal.
%
% In the code below can insert commands to code and decode residue
%
%
%
% set_name: name of the image set (e.g. 'AC_Boy_basketball')
% ModelID: model ID number (3: rigid, 6: affine, 8: projective)
% MAX_NUM_MODELS: maximum number of models for a single image (1, 2 or 3)
% ref_frm_number: reference frame number (e.g. 1)
% tgt_frm_number: target frame number (e.g. 2)
% **********************************************************
function residue = get_residue_regis_prediction_v3_based(previous_orig_rgb, ref , tgt, ModelID, MAX_NUM_MODELS, set_name, residue_dir, imgNo)
   
    % addpath('C:\Users\reji_\MatlabProjects\secondDWT\wt');
%     addpath('geom');
 
    out_folder = fullfile(residue_dir, 'data');
    if(~exist(out_folder, 'dir')), mkdir(out_folder); end

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


    % model file name
    mdl_dir= fullfile(residue_dir, "model");
    if(~exist(mdl_dir, 'dir')), mkdir(mdl_dir); end

%     file_model_name = [mdl_dir '\' set_name '_transform_list_' '_model_id_' num2str(ModelID) '_num_models_' num2str(MAX_NUM_MODELS) '.mat'];
    file_model_name =fullfile(mdl_dir, sprintf("%s_img%03d_transfrom_list_model_id_%d_num_models_%d.mat", set_name,imgNo, ModelID, MAX_NUM_MODELS));
    readModelFlag = isfile(file_model_name);


    % Encoder side function calls ****************************************
    % determine models used for warping
    [tformList] = enc_model(previous_orig_rgb, tgt, pix_min, pix_max, ModelID, MAX_NUM_MODELS, file_model_name);

    % determine quadtree and perform warping with models
    [pred, Index_Map, QTree_Depth, SMALLEST_BLK_SIZE] = enc_pred_warp(ref, tgt, pix_min, pix_max, tformList, inFlags);

    % write QuadTree decision and model index to bitsream
    qtBitStream_name = fullfile(out_folder, sprintf("qtBitStream_img%03d.bin",imgNo));
    writeQT_toBitstream(QTree_Depth, SMALLEST_BLK_SIZE, Index_Map, qtBitStream_name);

    % write out the transforms
    if readModelFlag == 0
        save(file_model_name, "tformList");
    end

    % subtract pred from target and write ppm residue file
    % residue is saved as uint16
    [residue] = enc_sub_pred_res(pix_min, pix_max, pred, tgt, Index_Map, inFlags, UP_SHIFT);
%     out_res_name = [out_folder 'residue.ppm'];
%     imwrite(residue, out_res_name, 'MaxValue', 65535);

    % *******************************************************************
    %
    % Can encode the residue to create decoded output 
    % residue.ppm -> Encoder -> Decoder -> dec_residue.ppm
    %
    % *******************************************************************
    
%     fprintf("output \n");
%     fprintf("model list MAT file, %s \t New_Model: %d \n", file_model_name, (1-readModelFlag));
%     fprintf("Quadtree bitsream, %s \n", qtBitStream_name);
% %     fprintf("Residue (uint16) PPM, %s \n", out_res_name);
%     fprintf("Residue left shift factor (uint16), %d \n", UP_SHIFT);
%     fprintf("Quadtree Depth %d, Smallest Block Size %d \n", QTree_Depth, SMALLEST_BLK_SIZE);
end

%% find or read motion models
function [tformList] = enc_model(in_ref, in_tgt, pix_min, pix_max, ModelID, MAX_NUM_MODELS, file_model_name)
    % get one or more transforms
    % check if transforms are stored and if available load transform
    % parameters. Otherwise determine transform list
    
    ref = rgb2gray(in_ref);
    tgt = rgb2gray(in_tgt);

    if isfile(file_model_name)
        inList = load(file_model_name, "tformList");
        tformList = inList.tformList;
        ListSize = size(tformList, 1);
    else
        pnts_ref = detectSIFTFeatures(ref,'ContrastThreshold',0.0333);%0.0133
        pnts_tgt = detectSIFTFeatures(tgt,'ContrastThreshold',0.0333);%0.0133
        
        % Display corners found in images A and B.
%         figure;
%         imshow(ref, [pix_min pix_max]);
%         hold on;
%         plot(pnts_ref);
%         title('Corners in ref');
        
%         figure;
%         imshow(tgt, [pix_min pix_max]);
%         hold on;
%         plot(pnts_tgt);
%         title('Corners in tgt');
%         
        % Extract descriptors for the corners
        [featrs_ref,pnts_ref] = extractFeatures(ref, pnts_ref);
        [featrs_tgt,pnts_tgt] = extractFeatures(tgt, pnts_tgt);
    
        [tformList, tformSpatialBoundsList, ListSize] = getTformList(featrs_ref, pnts_ref, featrs_tgt, pnts_tgt, ref, tgt, ModelID, MAX_NUM_MODELS);
    end

end
%% perform wapring of ref image
function [pred, Index_Map, QTree_Depth, SMALLEST_BLK_SIZE] = enc_pred_warp(in_ref, in_tgt, pix_min, pix_max, tformList, inFlags)
    
    % Turn on (1) DWT distortion metric or off (0) use image domain metric
    DWT_METRIC_FLAG = 0;
    % Turn on (1) search of model parameters or off (0) for no search
    NewModelFlag  = 0;
    % Replace INTRA with INTER_MV0 (DPCM mode - for testing purpose)
    COLOR_CONV_FLAG = inFlags.color_conv_flag;
    REPLACE_INTRA_DPCM = inFlags.replace_intra_with_dpcm;

    ref = double(rgb2gray(in_ref));
    tgt = double(rgb2gray(in_tgt));

    ListSize = size(tformList, 1);

    if (DWT_METRIC_FLAG == 1)
        tmp = matfile('synthGainSqrt16x16.mat');
        row16 = ceil(size(ref, 1)/16);
        col16 = ceil(size(ref, 2)/16);
        gainSqrt16 = repmat(tmp.synthGainSqrt, row16, col16);
        synthGainSqrt = gainSqrt16(1:size(tgt, 1), 1:size(tgt, 2));
    end

    % predict and find residue image in Y mode
    err = zeros(size(tgt, 1), size(tgt, 2), ListSize+1); % ref
    ref2 = upSample2(ref, pix_min, pix_max);
    
    % zero pred (i.e. intra mode)
    err(:,:, 1) = double(tgt);

    if (DWT_METRIC_FLAG == 1)
        A = err(:,:, 1);
        for lvl=0:3
            spacing = 2^lvl;
            offset = 0;
            A = analysis_9_7(analysis_9_7(A,spacing,offset)',spacing,offset)';
        end
        err(:,:, 1) = A .* synthGainSqrt;
    end
    
    % find residue error for each model
    for i = 1:ListSize
        tform = tformList(i); %get tform;
        tform2 = tform;
        tform2.A(1:2, 3) = tform2.A(1:2, 3) * 2; %2
        tform2.A(3, 1:2) = tform2.A(3, 1:2) * 0.5;
        
        ref2_w = imwarp(ref2,tform2,'OutputView',imref2d(size(tgt)*2)); % imref2d(size(ref2))
        
        ref_w = downSample2(ref2_w, pix_max, pix_min);
    
        % calculate residue error and the warped mask
        msk = uint8(ones(size(ref)));
        msk_w = imwarp(msk,tform,'OutputView',imref2d(size(tgt))); % imref2d(size(ref))
    
        % warped prediction and mask for where prediction is valid
        pred = smooth_roll_off(ref_w, msk_w(:,:,1));
    
        err(:,:, i+1) = double(tgt) - double(pred);
    
        if (DWT_METRIC_FLAG == 1)
            A = err(:,:, i+1);
            for lvl=0:3
                spacing = 2^lvl;
                offset = 0;
                A = analysis_9_7(analysis_9_7(A,spacing,offset)',spacing,offset)';
            end
        
            % multiply by sqrt of synth gain
            err(:,:, i+1) = A .* synthGainSqrt;
        end
        
    end
    
    LAMBDA = 1;
    MTN_MODEL_SIGNAL = 2;
    BRANCH_SIGNAL = 1;
    LEAF_SIGNAL = 1;
    SMALLEST_BLK_SIZE = 16;
    QTree_Depth = 4;
    Index_Map = zeros(size(tgt, 1), size(tgt, 2), QTree_Depth);
    MinCost = zeros(size(tgt, 1), size(tgt, 2), QTree_Depth);
    
    % find leaf nodes for each level and store cost
    % for each level: from finest (1) to coarsest (QTree_Depth)
    for QT_Level = 1:QTree_Depth
        % get summation of error for block_size for each model
    
        blk_size = SMALLEST_BLK_SIZE * 2^(QT_Level-1);
    
        scoreMetric = getBlockErrorSum(SMALLEST_BLK_SIZE, blk_size, err, size(tgt, 1), size(tgt, 2), ListSize, DWT_METRIC_FLAG);
        
        % detemine model index providing min sum error in blok_size 
        [M,I] = min(scoreMetric, [], 3);
    
        if QT_Level == 1
            r = 0;
            for row_start = 1:blk_size:size(tgt, 1)
                r = r + 1;
                c = 0;
                for col_start = 1:blk_size:size(tgt, 2)
                    c = c + 1;
    
                    row_end = row_start + blk_size - 1;
                    if row_end > size(tgt,1)
                        row_end = size(tgt,1);
                    end
            
                    col_end = col_start + blk_size - 1;
                    if col_end > size(tgt,2)
                        col_end = size(tgt,2);
                    end
    
                    MinCost(r,c,QT_Level) = M(r,c) + (LAMBDA * (MTN_MODEL_SIGNAL+LEAF_SIGNAL));
                    Index_Map(row_start:row_end, col_start:col_end,QT_Level) = I(r,c);
                end
            end
        else
            r = 0;
            for row_start = 1:blk_size:size(tgt, 1)
                r = r + 1;
                c = 0;
                for col_start = 1:blk_size:size(tgt, 2)
                    c = c + 1;
                    
                    row_end = row_start + blk_size - 1;
                    if row_end > size(tgt, 1)
                        row_end = size(tgt, 1);
                    end
            
                    col_end = col_start + blk_size - 1;
                    if col_end > size(tgt,2)
                        col_end = size(tgt,2);
                    end
    
                    clhd_r = (2*(r-1)) + 1;
                    clhd_c = (2*(c-1)) + 1;
                    childrenCost = sum(sum(MinCost(clhd_r:clhd_r+1, clhd_c:clhd_c+1, QT_Level-1))) + (LAMBDA * BRANCH_SIGNAL);
    
                    if childrenCost < (M(r,c) + (LAMBDA * (MTN_MODEL_SIGNAL+LEAF_SIGNAL)))
                        MinCost(r,c,QT_Level) = childrenCost;
                        Index_Map(row_start:row_end, col_start:col_end,QT_Level) = 0;
                    else
                        MinCost(r,c,QT_Level) = (M(r,c) + (LAMBDA * (MTN_MODEL_SIGNAL+LEAF_SIGNAL)));
                        Index_Map(row_start:row_end, col_start:col_end,QT_Level) = I(r,c);
                    end
                end
            end
        end
    end
    
    % modify each index map
    % for each level: from coarsest (QTree_Depth) to finest (1)
    Leaf_Map = zeros(size(tgt, 1), size(tgt, 2));
    
%     figure; imagesc(Index_Map(:,:,QTree_Depth));
%     title(['Quadtree Level ', num2str(QT_Level)]);

    for QT_Level = QTree_Depth:-1:2
    
        Leaf_Map = Leaf_Map + double(Index_Map(:,:,QT_Level) > 0);
        Leaf_Map = double(Leaf_Map(:,:) > 0);
    
        Index_Map(:,:,QT_Level-1) = (Index_Map(:,:,QT_Level-1) .* (1-Leaf_Map));
    
%         figure; imagesc(Index_Map(:,:,QT_Level-1));
%         title(['Quadtree Level ', num2str(QT_Level-1)]);
    end

    % create index map of block_size index decisions 
    I_Blk_Map = sum(Index_Map, 3);
    
    pred = double(zeros(size(tgt)));
    msk_mdl = double(zeros(size(tgt)));
    if COLOR_CONV_FLAG == 1
        ref = conv_srgb_to_xyb(double(in_ref));
        tgt = conv_srgb_to_xyb(double(in_tgt));
    else
        ref = double(in_ref);
        tgt = double(in_tgt);
    end
    ref2 = upSample2(ref, pix_min, pix_max);
    msk2 = upSample2(msk, 0, 1);
    
    % determine residue for each region based on region transform
    for i = 1:ListSize
        msk_mdl(:,:,1) = double(I_Blk_Map == i+1);
        msk_mdl(:,:,2) = double(I_Blk_Map == i+1);
        msk_mdl(:,:,3) = double(I_Blk_Map == i+1);
    
        tform = tformList(i); %get tform;
        tform2 = tform;
        tform2.A(1:2, 3) = tform2.A(1:2, 3) * 2; %2
        tform2.A(3, 1:2) = tform2.A(3, 1:2) * 0.5;
    
        if i >= 1 && NewModelFlag == 1
            % try altering transform parameters to minimize error for region/mask
            msk_w = double(imwarp(msk,tform,'OutputView',imref2d(size(tgt)))); % imref2d(size(ref))
            msk_combined = msk_mdl(:,:,1) .* msk_w;
            tform2 = optimizeTrans(tform2, ref2, tgt, msk_combined, synthGainSqrt, pix_min, pix_max);
    
            % copy new transform matrix, scaling appropriately specific terms
            tform.A(1:2, 1:2) = tform2.A(1:2, 1:2);
            tform.A(1:2, 3) = tform2.A(1:2, 3) * 0.5;
            tform.A(3, 1:2) = tform2.A(3, 1:2) * 2;
            tformList(i) = tform;
        end
    
        ref2_w = imwarp(ref2,tform2,'OutputView',imref2d(size(tgt) .* [2 2 1])); %imref2d(size(ref2))
        msk2_w = imwarp(msk2,tform2,'OutputView',imref2d(size(tgt) .* [2 2 1])); % imref2d(size(ref2))
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
    
       pred = pred + (ref_w .* msk_mdl); %((tgt - ref_w + pix_max) .* msk_mdl);
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

    % write out prediction as a pfm file
    pfm_folder = 'data';
    out_residue_name = fullfile(pfm_folder, 'pred.pfm');
    write_pfm(single(pred), out_residue_name, 'l', 0);
end

%% subtract prediction from target to generate residue
function [residue] = enc_sub_pred_res(pix_min, pix_max, pred, tgt, Index_Map, inFlags, UP_SHIFT)
    
    COLOR_CONV_FLAG = inFlags.color_conv_flag;
    REPLACE_INTRA_DPCM = inFlags.replace_intra_with_dpcm;

    if COLOR_CONV_FLAG == 1
        % write out residue as a pfm file
        tgt = conv_srgb_to_xyb(double(tgt));

        residue = single(double(tgt) - pred);
        pfm_folder = 'data\';
        out_residue_name = [pfm_folder 'residue.pfm'];
        write_pfm(residue, out_residue_name, 'l', 0);
    else
        if COLOR_CONV_FLAG == 2
            tgt = conv_rgb2jpegycbcr(tgt);
            pred = conv_rgb2jpegycbcr(uint8(pred));
        end
        pix_max = pix_max * 2^UP_SHIFT;
        tgt = double(tgt) * 2^UP_SHIFT;
        pred = double(pred) * 2^UP_SHIFT;

        % create index map of block_size index decisions 
        I_Blk_Map = sum(Index_Map, 3);
        offset = floor((pix_max+1)/2);

        if (REPLACE_INTRA_DPCM == 1)
            offset_map = repmat( ((I_Blk_Map > 0) * offset), 1, 1, 3 );
        else
            offset_map = repmat( ((I_Blk_Map > 1) * offset), 1, 1, 3 );
        end
    
        residue = uint16(floor(tgt - pred + offset_map + offset + 0.5));
    end

%     figure; imagesc(residue(:,:,1));
%     title('Residue Image');
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

%% match features and form transform matrix
function [tformList, tformSpatialBoundsList, ListSize] = getTformList(featrs_ref, pnts_ref, featrs_tgt, pnts_tgt, ref, tgt, ModelID, max_num_models)
    indexPairs = matchFeatures(featrs_ref,featrs_tgt);
    
    pnts_ref = pnts_ref(indexPairs(:,1),:);
    pnts_tgt = pnts_tgt(indexPairs(:,2),:);
    
%     figure;
%     showMatchedFeatures(ref,tgt,pnts_ref,pnts_tgt);
%     legend('ref','tgt');

    minFeaturePointsFlag = 1;
    listCnt = 1;
    % max_num_models = 2; %3;
    if ModelID == 3
        tformList(max_num_models,1) = rigidtform2d; 
    elseif ModelID == 6
        tformList(max_num_models,1) = affinetform2d; 
    elseif ModelID == 8
        tformList(max_num_models,1) = projtform2d;
    end
    tformSpatialBoundsList = zeros(max_num_models,5);

    while minFeaturePointsFlag
        if ModelID == 3
            [tform,inlierIdx, status] = estgeotform2d(pnts_ref,pnts_tgt,'rigid', MaxDistance=0.5, MaxNumTrials=3000); 
        elseif ModelID == 6
            [tform,inlierIdx, status] = estgeotform2d(pnts_ref,pnts_tgt,'affine', MaxDistance=0.5, MaxNumTrials=3000); 
        elseif ModelID == 8
            [tform,inlierIdx, status] = estgeotform2d_unsw(pnts_ref,pnts_tgt,'projective', MaxDistance=0.5, MaxNumTrials=3000);
        end
        pnts_inlier_ref = pnts_ref(inlierIdx,:);
        pnts_inlier_tgt = pnts_tgt(inlierIdx,:);
        ref_w = imwarp(ref,tform,'OutputView',imref2d(size(tgt))); % imref2d(size(ref))
        pnts_ref_w = transformPointsForward(tform,pnts_inlier_ref.Location);
        
%         figure;
%         showMatchedFeatures(tgt,ref_w,pnts_inlier_tgt,pnts_ref_w);
%         legend('ref','tgt');

        mean_x = mean(pnts_inlier_tgt.Location(:,1));
        mean_y = mean(pnts_inlier_tgt.Location(:,2));
        cov_x_y = cov(pnts_inlier_tgt.Location(:,1), pnts_inlier_tgt.Location(:,2));
       
        pnts_ref(inlierIdx,:)=[];
        pnts_tgt(inlierIdx,:)=[];

        tformList(listCnt) = tform;
        tformSpatialBoundsList(listCnt, :) = [mean_x, mean_y, cov_x_y(1,1), cov_x_y(2,2), cov_x_y(1,2)];

        if (size(pnts_ref,1) <= 10) || (listCnt == max_num_models)
            minFeaturePointsFlag = 0;
        end

        ListSize = listCnt;

        listCnt = listCnt + 1;

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

%% get summation of error for block_size for each model
function scoreMetric = getBlockErrorSum(SMALLEST_BLK_SIZE, blk_size, err, nrows, ncols, ListSize, DWT_METRIC_FLAG)
    r = 0;
    scoreMetric = zeros(nrows, ncols, ListSize+1);
    for row_start = 1:blk_size:nrows
        r = r + 1;
        c = 0;
        for col_start = 1:blk_size:ncols
            c = c + 1;
    
            row_end = row_start + blk_size - 1;
            if row_end > nrows
                row_end = nrows;
            end
    
            col_end = col_start + blk_size - 1;
            if col_end > ncols
                col_end = ncols;
            end
    
            % se = err(row_start:row_end, col_start:col_end, :) .* err(row_start:row_end, col_start:col_end, :);
            % LL_mask = ones(row_end - row_start + 1, col_end - col_start + 1);
            % LL_mask(1:SMALLEST_BLK_SIZE:end, 1:SMALLEST_BLK_SIZE:end) = 0;
            if DWT_METRIC_FLAG == 1
                se = double(abs(err(row_start:row_end, col_start:col_end, :)) > 10); %.* LL_mask
                scoreMetric(r,c, :) = sum(sum(se));
            else
                scoreMetric(r,c, :) = var(err(row_start:row_end, col_start:col_end, :), 1, [1 2]);
            end
        end
    end
end

%% change transform parameters to minimize error
function tform2 = optimizeTrans(tform2, ref2, tgt, msk_mdl, synthGainSqrt, pix_min, pix_max)

    %perform optimization for gray scale 8 bit unsigned
    scale = pix_max/255;
    ref2 = double(rgb2gray(uint8(ref2/scale)));
    tgt = double(rgb2gray(uint8(tgt/scale)));

    if isa(tform2, 'rigidtform2d')
        tx0(1) = double(acosd(tform2.A(1, 1)));
        tx0(2)  = double(tform2.A(1, 3));
        tx0(3)  = double(tform2.A(2, 3));
        tx0(4) = double(0);
        tx0(5) = double(0);
        tx0(6) = double(0);
    else
        tx0(1) = double(tform2.A(1, 1));
        tx0(2) = double(tform2.A(1, 2));
        tx0(3) = double(tform2.A(1, 3));
        tx0(4) = double(tform2.A(2, 1));
        tx0(5) = double(tform2.A(2, 2));
        tx0(6) = double(tform2.A(2, 3));
    end
%     tx0(1) = double(tform2.A(1, 3));
%     tx0(2) = double(tform2.A(2, 3));

    f = @(tx) getWarpPredSSE(tx, tform2, ref2, tgt, msk_mdl, synthGainSqrt, (pix_min/scale), (pix_max/scale));

    options = optimset('Display','iter', 'MaxIter', 50);

    [tx, fval] = fminsearch(f, tx0, options);

    if isa(tform2, 'rigidtform2d')
        tform2.RotationAngle = tx(1);
        tform2.Translation(1, 1) = tx(2);
        tform2.Translation(1, 2) = tx(3);
        tform2.R(1, 1) = cosd(tx(1));
        tform2.R(1, 2) = -sind(tx(1));
        tform2.R(2, 1) = sind(tx(1));
        tform2.R(2, 2) = cosd(tx(1));

        tform2.A(1, 1) = tform2.R(1, 1);
        tform2.A(1, 2) = tform2.R(1, 2);
        tform2.A(1, 3) = tform2.Translation(1, 1);
        tform2.A(2, 1) = tform2.R(2, 1);
        tform2.A(2, 2) = tform2.R(2, 2);
        tform2.A(2, 3) = tform2.Translation(1, 2);
    else
        tform2.A(1, 1) = tx(1);
        tform2.A(1, 2) = tx(2);
        tform2.A(1, 3) = tx(3);
        tform2.A(2, 1) = tx(4);
        tform2.A(2, 2) = tx(5);
        tform2.A(2, 3) = tx(6);
    end

%     tform2.A(1, 3) = tx(1);
%     tform2.A(2, 3) = tx(2);
end

function writeQT_toBitstream(QTree_Depth, SMALLEST_BLK_SIZE, Index_Map, streamFileName)
    
    byteStream = byteArrayStream(12000);

    % for each level: from coarsest (QTree_Depth) to finest (1)
    Leaf_Map = zeros(size(Index_Map, 1), size(Index_Map, 2));

    for QT_Level = QTree_Depth:-1:1

        blk_size = SMALLEST_BLK_SIZE * 2^(QT_Level-1);

        for row_start = 1:blk_size:size(Index_Map, 1)
            for col_start = 1:blk_size:size(Index_Map, 2)
                % if not pruned away then need to signal
                if Leaf_Map(row_start, col_start) == 0
                    % branch if zero else leaf
                    if Index_Map(row_start, col_start, QT_Level) == 0 
                        byteStream = byteStream.writeBits(0,1);
                    else
                        if QT_Level > 1
                            byteStream = byteStream.writeBits(1,1);
                        end
    
                        % write model ID - 4 options so two bits
                        modelID = Index_Map(row_start, col_start, QT_Level) - 1;
    
                        byteStream = byteStream.writeBits(modelID,2);
                    end
                end
            end
        end
                
        Leaf_Map = Leaf_Map + double(Index_Map(:,:,QT_Level) > 0);
        Leaf_Map = double(Leaf_Map(:,:) > 0);
    end

    % write out zeros to complete last byte and push left valid info
    byteStream = byteStream.endWriteBits();

    % write out bytes
    qtFile = fopen(streamFileName, 'w');
    fwrite(qtFile, byteStream.byteArray(1:byteStream.byteCnt));
    fclose(qtFile);

end
