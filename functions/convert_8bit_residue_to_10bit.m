function recon_residue = convert_8bit_residue_to_10bit(residue8bit, idx)
in8_C1  = residue8bit(:,:,1);
in8_C2 = residue8bit(:,:,2);
in8_C3 = residue8bit(:,:,3);
ss = size(in8_C1);
out16_C1 = int16(zeros(ss));
[h, w, ~] = size(residue8bit);
for m=1:h
    for n=1:w
        if(idx.i1_C1(m,n) == 1)
            out16_C1(m,n)= 510 - int16(in8_C1(m,n)) -128;
        elseif(idx.i2_C1(m,n) == 1)
            out16_C1(m,n) = 0 - int16(in8_C1(m,n)) - 128;
        else
            out16_C1(m,n) = int16(in8_C1(m,n)) - 128;
        end
    end
end
out16_C2 = int16(zeros(ss));
for m=1:h
    for n=1:w
        if(idx.i1_C2(m,n) == 1)
            out16_C2(m,n)= 510 - int16(in8_C2(m,n)) - 128;
        elseif(idx.i2_C2(m,n) == 1)
            out16_C2(m,n) = 0 - int16(in8_C2(m,n)) - 128;
        else
            out16_C2(m,n) = int16(in8_C2(m,n)) - 128;
        end
    end
end
out16_C3 = int16(zeros(ss));
for m=1:h
    for n=1:w
        if(idx.i1_C3(m,n) == 1)
            out16_C3(m,n)=510 - int16(in8_C3(m,n)) - 128;
        elseif(idx.i2_C3(m,n) == 1)
            out16_C3(m,n) = 0 - int16(in8_C3(m,n)) - 128;
        else
            out16_C3(m,n) = int16(in8_C3(m,n)) - 128;
        end
    end
end
%ch_size = size(y_cap_Y);
recon_C1 = out16_C1 ;
recon_C2 = out16_C2;
recon_C3 = out16_C3 ;
recon_residue = cat(3, int16(recon_C1), int16(recon_C2), int16(recon_C3));
end
