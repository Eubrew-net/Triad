function cloud
% Dentro de la carpeta Triad se debe de crear una carpeta con nombre BSRN
% donde se guyardan los archivos de radiaci�n ( pedirlos a Rosa) y el
% archivo de aerosles (descargarlode aeronet)


% la opci�n plot permite o no pintar gr�ficas
cloud_screening(fullfile(pwd,'BSRN'),'170101_171231_Izana.lev15','date_range',datenum(2017,1,[1 366]),'plot',0)
end
