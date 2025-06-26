function I  = imreadMPTiff(fname,colorChannel)

info = imfinfo(fname);
rows = info.Height; 
cols = info.Width;
 
nrFrames = length(info);
I  = zeros(rows,cols,nrFrames/colorChannel,colorChannel,'single');
frameCount=1;
 
for k = 1:colorChannel:nrFrames
    for c=1:colorChannel
        frame = imread(fname,k+c-1);
        I(:,:,frameCount,c)=frame(:,:,1);
    end
    frameCount=frameCount+1;
end
% info = imfinfo(fname);
% rows = info.Height; 
% cols = info.Width;
%  
% nrFrames = length(info);
% I  = zeros(rows,cols,nrFrames,colorChannel);
% frameCount=1;
%  
% for k = 1:colorChannel:nrFrames
%     for c=1:colorChannel
%         frame = imread(fname,k+c-1);
%         I(:,:,frameCount,c)=frame(:,:,1);
%     end
%     frameCount=frameCount+1;
% end