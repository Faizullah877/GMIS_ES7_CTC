% function to convert JPEG-YCbCr -> rgb
% code provided by FAIZULLAH <faiz@sju.ac.kr>
function rgb_image = conv_jpegycbcr2rgb(jpeg_ycbcr)
    %Input jpeg_ycbcr : uint8
    %output rgb_image : uint8
    T = inv([0.299,0.587,0.114;-0.168736,-0.331264,0.5;0.5,-0.418688,-0.081312]);
    ycbcr = double(jpeg_ycbcr);
    ycbcr(:,:,1) = ycbcr(:,:,1);
    ycbcr(:,:,2) = ycbcr(:,:,2)-127.5;
    ycbcr(:,:,3) = ycbcr(:,:,3)-127.5;
    rgb_image = zeros(size(ycbcr));
    
    rgb_image(:,:,1) = T(1)*ycbcr(:,:,1) + T(4)* ycbcr(:,:,2) + T(7)*ycbcr(:,:,3);
    rgb_image(:,:,2) = T(2)*ycbcr(:,:,1) + T(5)* ycbcr(:,:,2) + T(8)*ycbcr(:,:,3);
    rgb_image(:,:,3) = T(3)*ycbcr(:,:,1) + T(6)* ycbcr(:,:,2) + T(9)*ycbcr(:,:,3);
    rgb_image = uint8(rgb_image);
end