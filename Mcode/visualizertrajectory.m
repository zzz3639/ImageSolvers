function [ trajectory ] = visualizertrajectory( mus, w, s, zoom ,J)
%VISUALIZER Summary of this function goes here
%   Detailed explanation goes here

m=size(mus,3);
n=size(w,1);
ng=s*zoom;
trajectory=zeros(ng,ng,m/J);

for i=1:m/J
    for j=1:n
        x=zoom*mus(j,1,i*J)+zoom/2+1/2;
        y=zoom*mus(j,2,i*J)+zoom/2+1/2;
        x=floor(x);
        y=floor(y);
        
        %light(x,y,trajectory,w(j),i);
        if x>=1+1&&x<=ng-1&&y>=1+1&&y<=ng-1
            trajectory(x,y,i)=trajectory(x,y,i)+w(j,1,i*J);
            trajectory(x-1,y,i)=trajectory(x-1,y,i)+w(j,1,i*J)/2;
            trajectory(x+1,y,i)=trajectory(x+1,y,i)+w(j,1,i*J)/2;
            trajectory(x,y-1,i)=trajectory(x,y-1,i)+w(j,1,i*J)/2;
            trajectory(x,y+1,i)=trajectory(x,y+1,i)+w(j,1,i*J)/2;
        end
    end
end

end


