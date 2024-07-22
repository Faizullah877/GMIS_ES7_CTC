function result = get_result_values(M, set_name, codec, metric, values)
% set_name = Name of the set present at M.(set_name)
% codec = Anchor name present at M.(set_name).(codec)
% metric = metric name for which results are requested. 
%           Possible values are {"BPPs", "PSNR_Y", "PSNR_YUV",
%           "PSNR_HVS_M_Y" and "MSSSIM_Y"}
% value =   "average" for getting averaged values over all set images.
%           IMAGE_NAME for getting values for specific image.

set_result = M.(set_name);
codec_result = set_result.(codec);
if strcmp(values, "average")
    if strcmp(metric, "BPPs")
        result = codec_result.SET_AVG_RESULTS.SET_BPP;
    else
        result = codec_result.SET_AVG_RESULTS.(metric);
    end
else
    image_name = values;
    if strcmp(metric, "BPPs")
        result = codec_result.IMAGES_RESULT.BPPs.(image_name);
    else
        result = codec_result.IMAGES_RESULT.METRICS.(metric).(image_name);
    end
end

