% function string_with_number=Word2Number(string_without_number,library_set)
% 
% if nargin<2
%     library_set='ENGLISH';
% end
% 
% %constructing language library
% LIB.ENGLISH={'Zero','One','Two','Three','Four','Five','Six','Seven',...
%     'Eight','Nine'};
% LIB.SPANISH={'Cero','Uno','Dos','Tres','Cuatro','Cinco','Seis','Siete',...
%     'Ocho','Nueve'};
% LIB.JAPANESE={'Rei','Ichi','Nii','San','Shi','Go','Roku','Shichi',...
%     'Hachi','Kyuu'};
% 
% %select library
% Library=LIB.(upper(library_set));
% 
% 
% for m = 1:length(Library)
%     string_without_number=strrep(string_without_numbe,upper(Library{m}),...
%         num2str(m-1));
% end
% string_with_number=string_without_number;
% 
% end