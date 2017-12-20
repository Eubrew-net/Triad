
clear all;  
addpath(genpath('..\iberonesia\matlab'));
% addpath(genpath('.\temp'));

Cal.dir_figs=fullfile('.');

files=dir('*.mat')
for i=1:length(files)
    load(files(i).name);
    files(i).name
end

%% Mean differences
s=who('-regexp', '\<tablblind.*\d\>');
aux_blind=[];
for i=1:length(s)
    aux_blind=[aux_blind;evalin('base',s{i})];
end
aux_blind(isnan(aux_blind(:,3)),:)=[];
aux_blind(aux_blind(:,1)==157,:)=[];
aux_blind(aux_blind(:,1)==183,:)=[];
aux_blind(aux_blind(:,1)==185,:)=[];
figure; hist(aux_blind(:,end),-5:0.5:5);
% figure; h=nhist(aux_blind(:,end),-5:0.5:5);

s=who('-regexp', '\<tablfinal.*\d\>');
aux_final=[];
for i=1:length(s)
    aux_final=[aux_final;evalin('base',s{i})];
end
aux_final(isnan(aux_final(:,3)),:)=[]; aux_final(aux_final(:,1)==17,:)=[];
aux_final(aux_final(:,1)==157,:)=[];
aux_final(aux_final(:,1)==183,:)=[];
aux_final(aux_final(:,1)==185,:)=[];
% figure; hist(aux_final(:,3),-5:0.1:5);
% figure; h=nhist(aux_final(:,3),-5:0.1:5,'p');

%% Campaigns plot
% [a b]=grpstats(aux_blind(:,[1 2]),aux_blind(:,2),{'mean','numel'});
% camp=NaN*ones(size(a,1),4);
% camp(:,1)=a(:,2);
% camp(logical([1 0 0 0 1 0 0 1 0]),2)=b(logical([1 0 0 0 1 0 0 1 0])');
% camp(logical([0 0 1 0 0 0 1 0 0]),3)=b(logical([0 0 1 0 0 0 1 0 0])');
% camp(logical([0 1 0 1 0 1 0 0 1]),4)=b(logical([0 1 0 1 0 1 0 0 1])');
% 
% chk_bar

%%
% Quitamos los datos de Sodankyla
aux_blind(aux_blind(:,2)==734578 & aux_blind(:,1)==6,:)=[]; 
aux_blind(aux_blind(:,2)==734578 & aux_blind(:,1)==37,:)=[]; 

aux_final(aux_final(:,2)==734578 & aux_final(:,1)==6,:)=[]; 
aux_final(aux_final(:,2)==734578 & aux_final(:,1)==37,:)=[]; 

A={aux_blind(:,end),aux_final(:,3)};
figure; 
[tex y X]=nhist(A,'minbins',20,'maxbins',40,'separate','samebins','noerror',...
                  'xlabel','Ozone Deviation (%)','ylabel','Ocurrences'); 

ax=findobj(gcf,'Type','axes'); set(findobj(gcf,'Type','Line'),'LineWidth',1)
set(ax,'XTick',-4:1:6,'box','On','XGrid','On','YGrid','On');
axes(ax(2));text(3,12,'Initial Status','BackgroundColor','w','FontWeight','Bold');
axes(ax(1));text(3,30,'Final Days','BackgroundColor','w','FontWeight','Bold'); 
suptitle(sprintf('RBCC-E Intercomparison Campaigns\r\n2009 - 2015 (%d Ozone Cal.)',size(aux_blind,1)));

stats_blind=cat(2,round(100*length(find(abs(aux_blind(:,end))<1))/size(aux_blind,1)),...
                  round(100*length(find(abs(aux_blind(:,end))<0.5))/size(aux_blind,1)));
stats_final=cat(2,round(100*length(find(abs(aux_final(:,3))<1))/size(aux_final,1)),...
                  round(100*length(find(abs(aux_final(:,3))<0.5))/size(aux_final,1)));

fprintf('\r\nAll Campaigns Statistics\r\n');
displaytable(num2cell([stats_blind;stats_final]),...
             {'O3 Dev. < 1%','O3 Dev. < 0.5%'},15,'d',{'Initial Status','Final Days'}); 
         
%%
% printfiles_report(gcf,'.','aux_pattern',{'all'},'Height',9,'Width',13,'Format','tiff');
% close all

%% Tabla
idx=group_time(aux_blind(:,2),unique(aux_blind(:,2)));
aux_b=[]; idx_=unique(idx);
for j=1:length(idx_)
    aux_b=[aux_b,length(find(idx==idx_(j)))];  
end
blind=mat2cell(sortrows(aux_blind,2),aux_b,5);

idx=group_time(aux_final(:,2),unique(aux_final(:,2)));
aux_f=[]; idx_=unique(idx);
for j=1:length(idx_)
    aux_f=[aux_f,length(find(idx==idx_(j)))];  
end
final=mat2cell(sortrows(aux_final,2),aux_f,5);

data_blind=unique(aux_blind(:,1)); data_final=unique(aux_final(:,1));
for ii=1:length(final)
     summ_b=blind{ii}; summ_f=final{ii};  
               
     data_blind=scan_join(data_blind,summ_b(:,[1 end]));  
     data_final=scan_join(data_final,summ_f(:,[1 3]));  
end 

fprintf('\r\nAll Campaigns: Blind Days Statistics\r\n');
chk=num2cell(data_blind(:,2:end)); chk(cellfun(@(x) isnan(x),chk))={''};
displaytable(chk,cellstr(datestr(unique(aux_blind(:,2)),12))',7,'.2f',...
             cellfun(@(x) num2str(x),num2cell(unique(aux_blind(:,1))),'UniformOutput',0)); 
  
fprintf('\r\nAll Campaigns: Final Days Statistics\r\n');
chk=num2cell(data_final(:,2:end)); chk(cellfun(@(x) isnan(x),chk))={''};
displaytable(chk,cellstr(datestr(unique(aux_final(:,2)),12))',7,'.2f',...
             cellfun(@(x) num2str(x),num2cell(unique(aux_final(:,1))),'UniformOutput',0)); 
                
%% Variability
s=who('-regexp', '\<tablblind.*std\>');
aux_blind=[];
for i=1:length(s)
    aux_blind=[aux_blind;evalin('base',s{i})];
end
figure; hist(aux_blind(:,end),-5:0.5:5);

s=who('-regexp', '\<tablfinal.*std\>');
aux_final=[];
for i=1:length(s)
    aux_final=[aux_final;evalin('base',s{i})];
end
aux_final(aux_final(:,1)==17,:)=[];
figure; hist(aux_final(:,3),[-5:0.1:5])

% save('tabl_izo2011','-APPEND','tablblind_izo2011','tablblind_izo2011_std','tablfinal_izo2011','tablfinal_izo2011_std');
