function [error,desviacion]=histog_mensual(A2,A3,AM)

% Función empleada para el ajuste Histograma mensual y función gaussiana. 
% Futuro Paper de la Triada.

%Parámetros entrada
    
    % media=celda que contiene, media{}=[año, día, A2, A3, media].
    % para todas las medidas y las simultáneas.
    
    % EjeX1 y EjeX2 son los valores de corte para los histogramas

figure

% Eliminos meses donde el poco numero de datos no permite calcular la media
% mensual
A2new=[]; A3new=[]; AMnew=[];
for j=1:1:size(A2,1)
    if isnan(A2(j,3))==0; A2new=[A2new;A2(j,:)]; end
    if isnan(A3(j,3))==0; A3new=[A3new;A3(j,:)]; end
    if isnan(AM(j,3))==0; AMnew=[AMnew;AM(j,:)]; end
end

for j=1:1:3
    subplot(3,1,j)
    ejex=datenum(A2new(:,1),A2new(:,2),1);
    if j==1;plot(ejex,A2new(:,7),ejex,A2new(:,8),ejex,A2new(:,9),'linewidth',2); hold on; end
    if j==2;plot(ejex,A3new(:,7),ejex,A3new(:,8),ejex,A3new(:,9),'linewidth',2); hold on; end
    if j==3;plot(ejex,AMnew(:,7),ejex,AMnew(:,8),ejex,AMnew(:,9),'linewidth',2); end
    %startDate = datenum('01-05-2005');
    %endDate = datenum('01-01-2017');
    %xData = linspace(startDate,endDate,144);
    datetick('x','yyyy','keeplimits')
    xlim([datenum(2010,5,1),datenum(2019,6,1)])
    ylim([0.985 1.015])
end


%Calculos los std mensuales.
for j=1:1:3 % 1->A2, 2->A3, 3->AM  
error.mediaA2(1,j)=nanmean(A2new(:,6+j))
error.mediaA3(1,j)=nanmean(A3new(:,6+j))
error.mediaAM(1,j)=nanmean(AMnew(:,6+j))
 
desviacion.A2(1,j)=nanstd(A2new(:,6+j))
desviacion.A3(1,j)=nanstd(A3new(:,6+j))
desviacion.AM(1,j)=nanstd(AMnew(:,6+j))
end

% %Calculos los errores mensuales relativos.
% for j=1:1:6
% error_relativo.diario157(1,j)=mean(abs(Histograma{j}(:,6)./Histograma{j}(:,9)*100))
% error_relativo.diario183(1,j)=mean(abs(Histograma{j}(:,7)./Histograma{j}(:,9)*100))
% error_relativo.diario185(1,j)=mean(abs(Histograma{j}(:,8)./Histograma{j}(:,9)*100))
% 
% error_relativo.std157(1,j)=std(Histograma{j}(:,6)./Histograma{j}(:,9)*100)
% error_relativo.std183(1,j)=std(Histograma{j}(:,7)./Histograma{j}(:,9)*100)
% error_relativo.std185(1,j)=std(Histograma{j}(:,8)./Histograma{j}(:,9)*100)
% end
% 



% % arreglo datos para que en vez de ratio, aparezcan medias en las columnas 
% 
% subplot(3,4,4)
% 
% 
% 
% 
% A2gauss=[];A3gauss=[];AMgauss=[];
% for i=size(A2,1):-1:1
%     if sum(isnan(A2(i,3:5)))==0
%         A2gauss=[i,A2(i,7:9);A2gauss]
%     end
%     if sum(isnan(A3(i,3:5)))==0
%         A3gauss=[i,A3(i,7:9);A3gauss]
%     end
%     if sum(isnan(AM(i,3:5)))==0
%         AMgauss=[i,AM(i,7:9);AMgauss]
%     end
% end
% 
% nbins=20
% [count, centers]=hist(A2gauss(:,2:4),nbins)
% gaussEqn = 'a*exp(-((x-b)/c)^2)+1'
% startPoints = [15, 1, 0.03]
% f1 = fit(centers,count(:,1),gaussEqn,'Start',startPoints)
% f2 = fit(centers,count(:,2),gaussEqn,'Start',startPoints)
% f3 = fit(centers,count(:,3),gaussEqn,'Start',startPoints)
% a=f1.a;b=f1.b;c=f1.c;
% f11= a*exp(-((centers-b)/c).^2)+1
% a=f2.a;b=f2.b;c=f2.c;
% f22= a*exp(-((centers-b)/c).^2)+1
% a=f3.a;b=f3.b;c=f3.c;
% f33= a*exp(-((centers-b)/c).^2)+1
% subplot(3,4,4)
% plot(centers,f11,centers,f22,centers,f33); 
% 
% [count, centers]=hist(A3gauss(:,2:4),nbins)
% gaussEqn = 'a*exp(-((x-b)/c)^2)'
% startPoints = [15, 1, 0.03]
% f1 = fit(centers,count(:,1),gaussEqn,'Start',startPoints)
% f2 = fit(centers,count(:,2),gaussEqn,'Start',startPoints)
% f3 = fit(centers,count(:,3),gaussEqn,'Start',startPoints)
% a=f1.a;b=f1.b;c=f1.c;
% f11= a*exp(-((centers-b)/c).^2)+1
% a=f2.a;b=f2.b;c=f2.c;
% f22= a*exp(-((centers-b)/c).^2)+1
% a=f3.a;b=f3.b;c=f3.c;
% f33= a*exp(-((centers-b)/c).^2)+1
% subplot(3,4,8)
% plot(centers,f11,centers,f22,centers,f33); 
% 
% 
% [count, centers]=hist(AMgauss(:,2:4),nbins)
% gaussEqn = 'a*exp(-((x-b)/c)^2)+1'
% startPoints = [15, 1, 0.03]
% f1 = fit(centers,count(:,1),gaussEqn,'Start',startPoints)
% f2 = fit(centers,count(:,2),gaussEqn,'Start',startPoints)
% f3 = fit(centers,count(:,3),gaussEqn,'Start',startPoints)
% a=f1.a;b=f1.b;c=f1.c;
% f11= a*exp(-((centers-b)/c).^2)+1
% a=f2.a;b=f2.b;c=f2.c;
% f22= a*exp(-((centers-b)/c).^2)+1
% a=f3.a;b=f3.b;c=f3.c;
% f33= a*exp(-((centers-b)/c).^2)+1
% subplot(3,4,12)
% plot(centers,f11,centers,f22,centers,s); 
% 
% s=smooth(count(:,1))





% 
% 
%     plot(ejex,mensual(:,4),'*-','linewidth',2); hold on
%     plot(ejex,mensual(:,5),'*-','linewidth',2); hold on
%     
% if j==1; subplot(4,3,[1 3]); end
% if j==2; subplot(4,3,[5 7]); end
% if j==3; subplot(4,3,[9 11]); end
% plot
% 
% datos_Histograma=Histograma{j}(:,6:8)
% hist(datos_Histograma,nbins)
% [count, centers]=hist(datos_Histograma,nbins)
% xlim([-7 7])
% ylim([0 ejeX1])
% if j==3; legend('BR#157','BR#183','BR#185');title(sprintf('Figure 1'));end
% subplot(4,3,9+j)
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d'
% startPoints = [250, 0.1, 1, 0.1]
% f1 = fit(centers,count(:,1),gaussEqn,'Start',startPoints)
% f2 = fit(centers,count(:,2),gaussEqn,'Start',startPoints)
% f3 = fit(centers,count(:,3),gaussEqn,'Start',startPoints)
% 
% a=f1.a;b=f1.b;c=f1.c;d=f1.d;
% f11= a*exp(-((x-b)/c).^2)+d
% plot(x,f11); hold on
% 
% a=f2.a;b=f2.b;c=f2.c;d=f2.d;
% f22= a*exp(-((x-b)/c).^2)+d
% plot(x,f22); hold on
% 
% a=f3.a;b=f3.b;c=f3.c;d=f3.d;
% f33= a*exp(-((x-b)/c).^2)+d
% plot(x,f33); 
% 
% xlim([-7 7])
% ylim([0 ejeX2])
% end   
%  
% 
% mensual(85:144,:)=[]
% Emensual(1:60,:)=[]
% save Emensual.mat Emensual
% figure
% subplot(2,1,1)
% ejex=datenum(mensual(:,1),1,1);
% ejex=datenum(mensual(:,1),mensual(:,2),1);
% plot(ejex,mensual(:,3),'o-','linewidth',2); hold on
% plot(ejex,mensual(:,4),'*-','linewidth',2); hold on
% plot(ejex,mensual(:,5),'*-','linewidth',2); hold on
% legend('BR#157','BR#183','BR#185');
% ylabel('Total ozone (DU)')
% ylim([240 340])
% startDate = datenum('01-01-2005');
% endDate = datenum('01-01-2017');
% xData = linspace(startDate,endDate,84);
% datetick('x','yyyy','keeplimits')
% xlim([datenum(2010,1,1),datenum(2017,1,1)])
% 
% 
% %title(sprintf('O3 Monthly in the Period 2005-2015'));
% subplot(2,1,2)
% ejex(:,1)=datenum(Emensual(:,1),Emensual(:,2),1);
% plot(ejex,Emensual(:,4),'+-',ejex,Emensual(:,5),'o-',ejex,Emensual(:,6),'*-','linewidth',2);
% legend('BR#157','BR#183','BR#185');
% startDate = datenum('01-01-2010');
% endDate = datenum('01-01-2017');
% xData = linspace(startDate,endDate,84);
% datetick('x','yyyy','keeplimits')
% ylabel('Deviation (%)')
% % x=ejex;
% % datetick('x','yyyy','Keepticks','keeplimits');
% % xlim([datenum(2010,1,1),datenum(2016,12,31)])
% % legend('BR#157','BR#183','BR#185');
% xlim([datenum(2010,1,1),datenum(2017,1,1)])
% %title(sprintf('O3 Error Monthly in the Period 2005-2015'));
% 
% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',12);
% 
% 
% 
% 
% 
% a=findobj(gcf);
% 
% 
% figure
% subplot(1,2,1)
% datos_Histograma=Histograma{1}(:,6:8)
% nbins = 14;
% [count, centers]=hist(datos_Histograma,nbins)
% %title(sprintf('Daily Error respect to mean value in the Period 2005-2016'))
% title(sprintf('Difference respect to mean value in the Period 2010-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations (Days)')
% xlabel('Average daily difference (DU')
% xlim([-4 4])
% 
%  % get the handles associated with the current figure
% 
% 
% subplot(1,2,2)
% datos_Histograma=Histograma{2}(:,6:8)
% nbins = 20;
% hist(datos_Histograma,nbins)
% title(sprintf('A3  O3 Daily Error respect to mean value in the Period 2005-2016'))
% legend('BR#157','BR#183','BR#185');
% ylabel('Number of observations (Days)')
% xlim([-5 5])
end