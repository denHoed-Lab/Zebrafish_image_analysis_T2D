function bw2 = keepLargestn( bw )
%KEEPLARGESTN Summary of this function goes here
%   Detailed explanation goes here
bw2=zeros(size(bw));
if(sum(bw(:))==0)
    return;
end
bwl = bwlabeln(bw);
R = regionprops(bwl,'Area');
[val ind] = sort([R.Area],'descend');
maxind = ind(1);
inds=find(bwl==maxind);
bw2(inds)=1;
end

