function ImageBWout=RefineShape(ImageBW,ImageIn)
% Refine object shape with active contour method
%
%
  BW_tmp=ImageBW;
  L = bwlabel(ImageBW);
  L_prop=regionprops(L,"Centroid");
  for i_check=1:length(L_prop)
      if ImageBW(round(L_prop(i_check).Centroid(2)),round(L_prop(i_check).Centroid(1)))==1
          if(L(round(L_prop(i_check).Centroid(2)),round(L_prop(i_check).Centroid(1)))==0)
              continue;
          end
          tmp1=(L==L(round(L_prop(i_check).Centroid(2)),round(L_prop(i_check).Centroid(1))));
          tmp = activecontour(ImageIn, imresize(tmp1,[size(ImageIn,1) size(ImageIn,2)]),15,'Chan-Vese','ContractionBias',0);
          ImageBW(and(BW_tmp,tmp1))=0;
          ImageBW(and(BW_tmp,tmp))=1;
      end
  end
  
  ImageBWout=ImageBW;