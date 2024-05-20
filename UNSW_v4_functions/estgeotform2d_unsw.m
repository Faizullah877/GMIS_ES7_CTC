function [tform, inlierIndex, status] = estgeotform2d_unsw(matchedPoints1, matchedPoints2, transformType, varargin)
%

% Copyright 2022 The MathWorks, Inc.

%#codegen

reportError = (nargout ~= 3);
is2D = true;

[tform, inlierIndex, status] = ...
    algEstimateGeometricTransform_unsw(...
    matchedPoints1, matchedPoints2, transformType, reportError,...
    'estgeotform2d', is2D, varargin{:});

end