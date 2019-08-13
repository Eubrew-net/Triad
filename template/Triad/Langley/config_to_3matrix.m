function c=config_to_3matrix(config)
%
% convert cell config variable to a 3D matrix
% imput configuration cells config{Cal.n_inst};
% output 3d matrix with dimensions
% c=(icf fields, n_of_days, 3);
% last fiedl:
%  1 opertive configuration
%  2 alternative configuration
%  3 b file configuration file


cm=cell2mat(cxc);
c=reshape(cm,53,size(cx,1),3);
