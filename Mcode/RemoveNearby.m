function [ ansPic ] = RemoveNearby( Pic, Th )
% Modified by Zhang Haowen in 2015.08.26
%   given molecule list, returns list in which molecules with nearby
% molecules within Th are removed.
%   First two columns of Pic correspond to positions. Pic can contain other
% information in the following columns.
[up,vp]=sort(Pic(:,1));
PicSort=Pic(vp,:);
n=size(Pic,1);
Remove=zeros(n,1);
Th2=Th*Th;
for i=1:n
    for j=i:n
        if Pic(j,1)>Pic(i,1)+Th
            break;
        end
        if (Pic(i,1)-Pic(j,1))^2+(Pic(i,2)-Pic(j,2))^2<Th2
            Remove(i)=1;
            Remove(j)=1;
        end
    end
end

ansPic=Pic(find(1-Remove),:);

end

