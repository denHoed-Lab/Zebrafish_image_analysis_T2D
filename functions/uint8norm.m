function [I Imin Imax]=uint8norm(I)
I=single(I);
Imax=max(I(:));
Imin=min(I(:));
range=[Imin Imax];
I=(I-Imin)/(Imax-Imin)*255;

