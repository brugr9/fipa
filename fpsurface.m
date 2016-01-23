% Discrete fpsurface S for gray-levels g
[greyMax, greyMaxLoc] = max(I(:));
[greyMin, greyMinLoc] = min(I(:));
g = greyMax - greyMin;
S = zeros(sizeX, sizeY);
for x = 1:sizeX
    for y = 1:sizeY
        S(x,y) = g-1-I(x,y);
    end;
end;
mesh(S)
colormap(hot)
figure; imshow(S, []);         axis off; title('Surface');    hold off;