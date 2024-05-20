function  image_metrics = evaluate_images_IQA(ref_rgb, decoded_rgb, IQA_m_list)
% This function calculates Objective quality metrices values for input
% ref_rgb and decoded_rgb
% ref_rgb = Original image in RGB colorspace
% metric = image objective quality metric
%           allowed values are
%           "PSNR_Y", 'MSSSIM_Y",
%           "PSNR_HVS_M_Y", "PSNR_YUV"
% decoded_rgb = reconstructed image in RGB colorspace
% developed by Faiz Ullah (faiz@sju.ac.kr)

no_of_metrics = numel(IQA_m_list);

image_metrics = struct();
for m = 1: no_of_metrics
    metric_name = IQA_m_list(m);
    if(metric_name == "PSNR_Y")
        Orig_yuv = rgb2yuv_bt709(ref_rgb);
        Deco_yuv = rgb2yuv_bt709(decoded_rgb);
        metric_v = psnr(Orig_yuv(:,:,1), Deco_yuv(:,:,1));
    elseif(metric_name == "MSSSIM_Y")
        Orig_yuv = rgb2yuv_bt709(ref_rgb);
        Deco_yuv = rgb2yuv_bt709(decoded_rgb);
        metric_v = mean(squeeze(multissim(Deco_yuv(:,:,1), Orig_yuv(:,:,1))));
    elseif(metric_name == "PSNR_HVS_M_Y")
        Orig_yuv = rgb2yuv_bt709(ref_rgb);
        Deco_yuv = rgb2yuv_bt709(decoded_rgb);
        [p_hvs_m, ~] = psnrhvsm(Orig_yuv(:,:,1), Deco_yuv(:,:,1), 8);
        metric_v = p_hvs_m;
    elseif(metric_name == "PSNR_YUV")
        Orig_yuv = rgb2yuv_bt709(ref_rgb);
        Deco_yuv = rgb2yuv_bt709(decoded_rgb);
%         metric_v =  psnr(Orig_yuv, Deco_yuv);
        psnr_y = psnr(Orig_yuv(:,:,1), Deco_yuv(:,:,1));
        psnr_u = psnr(Orig_yuv(:,:,2), Deco_yuv(:,:,2));
        psnr_v = psnr(Orig_yuv(:,:,3), Deco_yuv(:,:,3));
        metric_v = ((psnr_y*6)+psnr_u + psnr_v)/8;
    else
        error(['Image Quality Metric ' IQA_m_list 'is not supported.']);
    end
    image_metrics.(metric_name) = metric_v;

end

end


