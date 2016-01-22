function g = movingthresh(I, n, k)
%movingthresh Image segmentation using a moving average threshold.
%   g = movingthresh(f, n, k) segments image I by thresholding its
%   intensities based on the moving average of the intesities along
%   individual rows of the image. The average at pixel k is formed by
%   averaging the intensities of that pixel and its n-1 preceding
%   neighbours. To reduce shading bias, the scanning is done in a zig-zag
%   manner, treating the pixels as if they were a 1-D, continuous stream.
%   If the value of the image at a point exceeds K percent of the value of
%   the running average at that point, a 1 is output in that location in G.
%   Otherwise a 0 is output. At the end of the procedure, G is thus the
%   thresholded (segmented) image. k must be a scalar in the range [0, 1].
%
%   f image
%   n pixel neighbours
%   k percent, scalar in range [0,1]
%   
%   Roland Bruggmann, BSc student Information Technology
%   Specialisation in Computer Perception and Virtual Reality CPVR
%   <mailto:roland.bruggmann@students.bfh.ch>
% 
%   Bern University of Applied Sciences, Engineering and Information Technology
%   Biel/Bienne, January 2016
% 
%  
%   References:
%   Rafael C. Gonzales, Richard E. Woods und Steven L. Eddins. 
%   Digital Image Processing using MATLAB(R). 2. Ed. 
%   Subsec. 11.3.7 Image Thresholding Using Moving Averages, P. 575 -- 578. 
%   Prentice Hall, 2004. isbn: 978-0-982-08540-0.

%%  
% *Moving averages:*
% 
% $$m(k+1) = \frac{1}{n} \sum_{i=k+2-n}^{k+1} z_i$

    % Preliminaries.
    I = tofloat(I);
    [M,N] = size(I);
    if (n < 1) || (rem(n,1) ~= 0)
        error('n must be an integer >= 1.');
    end
    if (k < 0) || (k > 1)
        error('K must be a fraction in the range [0,1].');
    end

    % Flip every other row of f to produce the equivalent of a zig-zag scanning
    % pattern. Convert image to a vector.
    I(2:2:end, :) = fliplr(I(2:2:end, :));
    I = I'; % Still a matrix.
    I = I(:)'; % Convert to row vector for use in  function filter.

    % Compute the moving average.
    maf = ones(1,n)/n; % The 1-D moving average filter.
    ma = filter(maf,1,I); % Computation of moving average.

    % Perform thresholding.
    g = I > k * ma;

    % Go back to image format (indexed subscripts).
    g = reshape(g,N,M)';
    % Flip alternate rows back.
    g(2:2:end, :) = fliplr(g(2:2:end, :));

end

