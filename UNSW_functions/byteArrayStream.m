classdef byteArrayStream
    properties
        byteArray
        bitCnt {mustBeNumeric}
        byteCnt {mustBeNumeric}
    end

    methods
        function obj = byteArrayStream(totalBytes)
             obj.byteArray = zeros(totalBytes, 1, "uint8");
             obj.bitCnt = uint8(0);
             obj.byteCnt = uint32(0);
        end
        function obj = writeBits(obj, inVal, inBits)
            for i = 1:inBits
                if mod(obj.bitCnt, 8) == 0
                    obj.byteCnt = obj.byteCnt + 1;
                    obj.bitCnt = uint8(0);
                end

                rghtShift = int32(inBits - i);
                bitVal = uint8(bitand(bitshift(inVal, -rghtShift), 0x00000001));
                obj.byteArray(obj.byteCnt) = bitor(bitshift(obj.byteArray(obj.byteCnt), 1), bitVal);
                
                obj.bitCnt = obj.bitCnt + 1;
            end
        end
        function obj = endWriteBits(obj)
            inBits = 8 - obj.bitCnt;
            
            for i = 1:inBits
                bitVal = 0;
                obj.byteArray(obj.byteCnt) = bitor(bitshift(obj.byteArray(obj.byteCnt), 1), bitVal);
                obj.bitCnt = obj.bitCnt + 1;
            end
        end
        function [obj, outVal] = readBits(obj, outBits)
            outVal = 0;
            for i = 1:outBits
                if mod(obj.bitCnt, 8) == 0
                    obj.byteCnt = obj.byteCnt + 1;
                    obj.bitCnt = uint8(0);
                end

                rghtShift = int32(8 - (obj.bitCnt+1));
                bitVal = uint8(bitand(bitshift(obj.byteArray(obj.byteCnt), -rghtShift), 0x01));
                outVal = bitor(bitshift(outVal, 1), bitVal);
                
                obj.bitCnt = obj.bitCnt + 1;
            end
        end

    end
end

