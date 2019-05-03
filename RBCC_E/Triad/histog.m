function [Histograma,error,error_relativo]=histog(media,ejeX1,ejeX2,ejeX3,ejeX4)

% Función empleada para el ajuste Histograma y función gaussiana. 
% Futuro Paper de la Triada.

%Parámetros entrada
    
    % media=celda que contiene, media{}=[año, día, A2, A3, media].
    % para todas las medidas y las simultáneas.
    
    % EjeX1 y EjeX2 son los valores de corte para los histogramas
    
    
%Paso 1: Calculo el valor medio 

% Todas las medidas!!! Coeficiente A2, A3 y la media
aux{1}=[media{1}(:,1:3),media{2}(:,3),media{3}(:,3)] % Coeficiente A2 (ajuste Fieletov), día y año.
aux{2}=[media{1}(:,1:2),media{1}(:,4),media{2}(:,4),media{3}(:,4)]%Coeficiente A3 (ajuste Fieletov), día y año.
aux{3}=[media{1}(:,1:2),media{1}(:,5),media{2}(:,5),media{3}(:,5)]% media, día y año.

% Las medidas simultaneas!!! Coeficiente A2, A3 y la media
% Periodo  a partir del 2010.
aux{4}=[media{4}(:,1:3),media{5}(:,3),media{6}(:,3)] % Coeficiente A2 (ajuste Fieletov), día y año.
aux{5}=[media{4}(:,1:2),media{4}(:,4),media{5}(:,4),media{6}(:,4)]%Coeficiente A3 (ajuste Fieletov), día y año.
aux{6}=[media{4}(:,1:2),media{4}(:,5),media{5}(:,5),media{6}(:,5)]%media, día y año.

% Histograma=[año,día,O3_157,O3_183,O3_185,dif_157,dif_183,dif_185,media_valores_ozono]
Histograma{1}=[]; Histograma{2}=[]; Histograma{3}=[];
Histograma{4}=[]; Histograma{5}=[]; Histograma{6}=[];
descartes=0
for j=1:1:6
    b=[]
    for i=1:1:size(aux{j},1)
        if sum(isnan(aux{j}(i,:)))==0
            a=mean(aux{j}(i,3:5))
            %if abs(aux{j}(i,3)-a)<=5 && abs(aux{j}(i,4)-a)<=5 && abs(aux{j}(i,5)-a)<=5
            b=[aux{j}(i,1:2),aux{j}(i,3),aux{j}(i,4),aux{j}(i,5),aux{j}(i,3)-a,aux{j}(i,4)-a,aux{j}(i,5)-a,a;b];
         %end
        end
    end
    Histograma{j}=b   
end

% figure
% title(sprintf('Figure 1'))
% A=0 ; B=6;
% 
% x=[-5:0.1:5]'
% for j=1:1:3
% subplot(4,3,[A+j B+j])
% datos_Histograma=Histograma{j}(:,6:8)
% %nbins = 20;
% xvalues1 = -5:1:5
% hist(datos_Histograma,xvalues1)
% [count, centers]=hist(datos_Histograma,xvalues1)
% xlim([-5 5])
% ylim([0 ejeX1])
% if j==3; legend('BR#157','BR#183','BR#185');title(sprintf('Figure 1'));end
% subplot(4,3,9+j)
% %gaussEqn = 'a*exp(-((x-b)/c)^2)+d'
% %startPoints = [250, 0.1, 1, 0.1]
% %f1 = fit(centers,count(:,1),gaussEqn,'Start',startPoints)
% %f2 = fit(centers,count(:,2),gaussEqn,'Start',startPoints)
% %f3 = fit(centers,count(:,3),gaussEqn,'Start',startPoints)
% 
% %a=f1.a;b=f1.b;c=f1.c;d=f1.d;
% %f11= a*exp(-((x-b)/c).^2)+d
% %f11=smooth(f11)
% plot(centers,smooth(count(:,1)),'b'); hold on
% 
% % a=f2.a;b=f2.b;c=f2.c;d=f2.d;
% % f22= a*exp(-((x-b)/c).^2)+d
% % f22=smooth(f22)
% % plot(x,f22); hold on
% plot(centers,count(:,2)); hold on
% 
% % a=f3.a;b=f3.b;c=f3.c;d=f3.d;
% % f33= a*exp(-((x-b)/c).^2)+d
% % f3=smooth(f33)
% % plot(x,f33); 
% plot(centers,count(:,3)); hold on
% xlim([-5 5])
% ylim([0 ejeX2])
% end
% 
% % a=findobj(gcf); % get the handles associated with the current figure
% % 
% % allaxes=findall(a,'Type','axes');
% % alllines=findall(a,'Type','line');
% % alltext=findall(a,'Type','text');
% % set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
% %     'FontSize',12);
% 
% % Repetimos figura pero para las medidas simultaneas
% 
% figure
% title(sprintf('Figure 2'))
% A=0 ; B=6;
% nbins = 20;
% x=[-5:0.1:5]'
% pos=1
% for j=4:1:6
% nbins = 20;
% subplot(4,3,[A+pos B+pos])
% datos_Histograma=Histograma{j}(:,6:8)
% xvalues1 = -5:1:5
% hist(datos_Histograma,xvalues1)
% [count, centers]=hist(datos_Histograma,xvalues1)
% xlim([-7 7])
% ylim([0 ejeX3])
% if j==4; legend('BR#157','BR#183','BR#185');title(sprintf('Figure 2'));end
% subplot(4,3,9+pos)
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d'
% startPoints = [250, 0.1, 1, 0.1]
% f1 = fit(centers,count(:,1),gaussEqn,'Start',startPoints)
% f2 = fit(centers,count(:,2),gaussEqn,'Start',startPoints)
% f3 = fit(centers,count(:,3),gaussEqn,'Start',startPoints)
% 
% % a=f1.a;b=f1.b;c=f1.c;d=f1.d;
% % f11= a*exp(-((x-b)/c).^2)+d
% % f11=smooth(f11)
% % plot(x,f11); hold on
% plot(centers,count(:,1)); hold on
% 
% % a=f2.a;b=f2.b;c=f2.c;d=f2.d;
% % f22= a*exp(-((x-b)/c).^2)+d
% % f22=smooth(f22)
% % plot(x,f22); hold on
% plot(centers,count(:,2)); hold on
% 
% % a=f3.a;b=f3.b;c=f3.c;d=f3.d;
% % f33= a*exp(-((x-b)/c).^2)+d
% % f33=smooth(f33)
% % plot(x,f33); 
% plot(centers,count(:,3)); hold on
% 
% xlim([-5 5])
% ylim([0 ejeX4])
% pos=pos+1;
% end

figure
title(sprintf('Figure 3'))
A=0 ; B=1;
x=[-5:0.1:5]'
pos=[1,3,5,2,4,6]

for j=1:1:6    
subplot(3,2,pos(1,j))
datos_Histograma=Histograma{j}(:,6:8)
xvalues1 = -5:1:5
hist(datos_Histograma,xvalues1)
[count, centers]=hist(datos_Histograma,xvalues1)
xlim([-5 5])

%ylim([0 ejeX3])
if j==3; xlabel('DU');end
if j==6; xlabel('DU');end
end

%Calculos los errores diarios.
for j=1:1:6
error.diario157(1,j)=mean(abs(Histograma{j}(:,6)))
error.diario183(1,j)=mean(abs(Histograma{j}(:,7)))
error.diario185(1,j)=mean(abs(Histograma{j}(:,8)))

error.std157(1,j)=std(Histograma{j}(:,6))
error.std183(1,j)=std(Histograma{j}(:,7))
error.std185(1,j)=std(Histograma{j}(:,8))
end

%Calculos los errores diarios relativos.
for j=1:1:6
error_relativo.diario157(1,j)=mean(abs(Histograma{j}(:,6)./Histograma{j}(:,9)*100))
error_relativo.diario183(1,j)=mean(abs(Histograma{j}(:,7)./Histograma{j}(:,9)*100))
error_relativo.diario185(1,j)=mean(abs(Histograma{j}(:,8)./Histograma{j}(:,9)*100))

error_relativo.std157(1,j)=std(Histograma{j}(:,6)./Histograma{j}(:,9)*100)
error_relativo.std183(1,j)=std(Histograma{j}(:,7)./Histograma{j}(:,9)*100)
error_relativo.std185(1,j)=std(Histograma{j}(:,8)./Histograma{j}(:,9)*100)
end

% a=findobj(gcf); % get the handles associated with the current figure
% 
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',12);
end