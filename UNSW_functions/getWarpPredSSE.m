function sse = getWarpPredSSE(tx, tform2, ref2, tgt, msk_mdl, synthGainSqrt, pix_min, pix_max)
    msk2 = double(ones(size(ref2(:,:,1))));

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

    ref2_w = imwarp(ref2,tform2,'OutputView',imref2d(size(tgt) * 2)); % imref2d(size(ref2))
    msk2_w = imwarp(msk2,tform2,'OutputView',imref2d(size(tgt) * 2)); % imref2d(size(ref2))
    ref2_w = smooth_roll_off(ref2_w, msk2_w);
    
    ref_w = downSample2(ref2_w, pix_max, pix_min);

    A = ref_w;
    for lvl=0:3
        spacing = 2^lvl;
        offset = 0;
        A = analysis_9_7(analysis_9_7(A,spacing,offset)',spacing,offset)';
    end
    pred = A;

    A = tgt;
    for lvl=0:3
        spacing = 2^lvl;
        offset = 0;
        A = analysis_9_7(analysis_9_7(A,spacing,offset)',spacing,offset)';
    end

    tgt = A;

    residue = double((abs(tgt - pred) .* msk_mdl) > 10);

    %err = residue .* residue;
    %residue(residue>pix_max/2) = pix_max/2;
    err = residue;

    sse = sum(err(:));

end

%% perform a smooth roll off at warpped region boundaries
function [smth] = smooth_roll_off(ref_w, in_msk_w)
    out_w = double(ref_w) .* double(in_msk_w);
    sum_w = double(in_msk_w);
    
    out_w = imgaussfilt(out_w, 1);
    sum_w = imgaussfilt(sum_w, 1);

    new_mask_w = double(sum_w == 1);

    smth = zeros(size(ref_w));
    for c = 1:1:size(ref_w, 3)
        smth(:,:,c) = (new_mask_w(:,:) .* double(ref_w(:,:,c))) + ((1-new_mask_w(:,:)) .* out_w(:,:,c));
    end
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
        down2(down2(:,:,i)>pix_max) = pix_max;
        down2(down2(:,:,i)<pix_min) = pix_min;
    
        %figure; imagesc(out2);
    end
end