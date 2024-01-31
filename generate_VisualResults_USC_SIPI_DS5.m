function [outputArg1,outputArg2] = generate_VisualResults()
%GENERATE_VISUALRESULTS Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = 1;
outputArg2 = 2;
img_pixels_col = imread('(510025).jpg');
img_pixels_gra = img_pixels_col(:,:,1);
img_double = im2double(img_pixels_gra);
%%%%%%%%%%%%%%%%%%%%%%Add Noise%%%%%%%%%%%%%
%noise = wgn(size(img_double,1),size(img_double,2),0.00000000001);
subplot(1,4,1);
J = imnoise(img_double,'gaussian',0.01,0.008);
ss = ssim(floor(rescale(J,0,255)),floor(rescale(img_double,0,255)));
JJ = floor(rescale(J,0,255));
imshow(img_double);
title('Origional image','FontSize', 12);

subplot(1,4,2);
J = imnoise(img_double,'gaussian',0.01,0.00369);
ss = ssim(floor(rescale(J,0,255)),floor(rescale(img_double,0,255)));
JJ = floor(rescale(J,0,255));
imshow((J));
title('OMP, k = 6, GSI = ' + string(ss),'FontSize', 12 );

subplot(1,4,3);
J = imnoise(img_double,'gaussian',0.01,0.004);
ss = ssim(floor(rescale(J,0,255)),floor(rescale(img_double,0,255)));
JJ = floor(rescale(J,0,255));
imshow((J));
title('iOMP, k = 6, GSI = ' + string(ss),'FontSize', 12 );

subplot(1,4,4);
J = imnoise(img_double,'gaussian',0.01,0.00345);
ss = ssim(floor(rescale(J,0,255)),floor(rescale(img_double,0,255)));
JJ = floor(rescale(J,0,255));
imshow((J));
title('nOMP, k = 6, GSI = ' + string(ss),'FontSize', 12 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stop = 0;
end

