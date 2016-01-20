function directionmap(D, s, I)
% DIRECTIONMAP - plot direction map
%
%  Usage:   directionmap(D, s, I)
%
%        D        - Direction image.
%        s        - Subsampling interval.
%        I        - Fingerprint image to underlay.
%
%  Roland Bruggmann, BSc student Information Technology
%  Specialisation in Computer Perception and Virtual Reality CPVR
%  <mailto:roland.bruggmann@students.bfh.ch>
% 
%  Bern University of Applied Sciences, Engineering and Information Technology
%  Biel/Bienne, January 2016
% 
%  
%  References:
%  Peter Kovesi, url: http://www.mathworks.com/matlabcentral/fileexchange/29280-fingerprint-matching-algorithm-using-shape-context-and-orientation-descriptors/content/sc_minutia/plotridgeorient.m
    
    [r, c] = size(D);
    l = 0.8*s;  % length of lines

    % Subsample the directions
    sDir = D(s:s:r-s, s:s:c-s);
    xoff = l/2*cos(sDir);
    yoff = l/2*sin(sDir);
    
    % Determine placement of direction vectors
    [x,y] = meshgrid(s:s:c-s, s:s:r-s);
    x = x-xoff;
    y = y-yoff;
    
    % Orientation vectors
    u = xoff*2;
    v = yoff*2;
    
    % plot
    imshow(I); alpha .5; hold on;
    quiver(x,y,u,v,0,'.','linewidth',1);
    axis off;
    hold off;
