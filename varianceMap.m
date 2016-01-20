function [varianceImg]=varianceMap(I, w, t)
% VARIANCEMAP - Compute Variance map of the image
%
% Usage:  [varianceImg] = varianceMap(inputImage, windowSize, thresh)
%
% Arguments:  I                 - A variance input image.
%             w                 - The size of the window to work with.
%             t                 - Threshold, typical value is 140.
% 
% Returns:    varianceImg       - The variance Map.
%
% With a fingerprint image at a 'standard' resolution of 500dpi suggested
% parameter values might be:
%
%    [varianceImg] = varianceMap(I, 5, 9);
%
% Author: Roland Bruggmann
% Created: January 2016
%
% Reference:
%   
% Siddhant Ahuja, 2008
% https://siddhantahuja.wordpress.com/2009/06/08/compute-variance-map-of-an-image/
%

if (mod(w,2)==0)
    error('The window size must be an odd number.');
end
r=(w-1)/2;

I=double(I);
[nr,nc] = size(I);
meanImg=zeros(nr, nc);
varianceImg=zeros(nr, nc);

% Compute a map of mean values
for i=1+r:1:nr-r
    for j=1+r:1:nc-r
        sum=0.0;
        for a=-r:1:r
            for b=-r:1:r
                sum=sum+I(i+a,j+b);
            end
        end
        meanImg(i,j)=sum/(w*w);
    end
end
% Compute a map of variance values
for i=1+r:1:nr-r
    for j=1+r:1:nc-r
        sum=0.0;
        for a=-r:1:r
            for b=-r:1:r
                sum=sum+((I(i+a,j+b)-meanImg(i,j))^2);
            end
        end         
        var=sum/((w*w)-1);
        % Apply threshold to produce a binarized variance map
        if (var > t)
            varianceImg(i,j) = 255;
        else
            varianceImg(i,j) = 0;
        end
    end
end
end
