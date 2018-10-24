% Rutina que copia y pega los archivos UV y B de cada brewer, desde el
% directorio RBCC_E al directorio del Brewer_proccesing sotfware.

% Este programa solo actualiza el año actual. Si se quiere actualizar años
% previos, la variable Year debe ser modificada con el año deseado 
% Year=2011, por ejemplo
% Sergio Leon

% Averiguamos el año.
a=datevec(now)
Year=a(1,1) 
% Year=2011
year=num2str(Year)

%Brewer#157
%Bfile
origen=strcat('C:\CODE\iberonesia\RBCC_E\',year,'\bdata157\B*.157')
destino=strcat('C:\Brewer_Processing_Software\data\br#157\data\',year)
copyfile(origen,destino)

%UV file
origen=strcat('C:\CODE\iberonesia\RBCC_E\',year,'\bdata157\UV\UV*.157')
destino=strcat('C:\Brewer_Processing_Software\data\br#157\data\',year)
copyfile(origen,destino)
%La instrucciones anteriores copian el arcchivo UVOAVG, lo borramos de la
%carpeta del Brewer_processing_sotfware
del=strcat('C:\Brewer_Processing_Software\data\br#157\data\',year,'\UVOAVG.157')
delete(del)

%Brewer#183
%Bfile
origen=strcat('C:\CODE\iberonesia\RBCC_E\',year,'\bdata183\B*.183')
destino=strcat('C:\Brewer_Processing_Software\data\br#183\data\',year)
copyfile(origen,destino)

%UV file
origen=strcat('C:\CODE\iberonesia\RBCC_E\',year,'\bdata183\UV\UV*.183')
destino=strcat('C:\Brewer_Processing_Software\data\br#183\data\',year)
copyfile(origen,destino)
%La instrucciones anteriores copian el arcchivo UVOAVG, lo borramos de la
%carpeta del Brewer_processing_sotfware
del=strcat('C:\Brewer_Processing_Software\data\br#183\data\',year,'\UVOAVG.183')
delete(del)

%Brewer#185
%Bfile
origen=strcat('C:\CODE\iberonesia\RBCC_E\',year,'\bdata185\B*.185')
destino=strcat('C:\Brewer_Processing_Software\data\br#185\data\',year)
copyfile(origen,destino)

%UV file
origen=strcat('C:\CODE\iberonesia\RBCC_E\',year,'\bdata185\UV\UV*.185')
destino=strcat('C:\Brewer_Processing_Software\data\br#185\data\',year)
copyfile(origen,destino)
%La instrucciones anteriores copian el arcchivo UVOAVG, lo borramos de la
%carpeta del Brewer_processing_sotfware
del=strcat('C:\Brewer_Processing_Software\data\br#185\data\',year,'\UVOAVG.185')
delete(del)

