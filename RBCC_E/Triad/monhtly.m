function []=monhtly(media)

figure
% ejex=datenum(B157(:,1),1,B157(:,2),0,B157(:,1),0) 
% plot(ejex,B157(:,4)); hold on
   
ejex=datenum(media{1}(:,1),1,media{1}(:,2))
plot(ejex,media{1}(:,3),'r')

for j=1:1:3
    

