function string_without_number=Number2Word(string_with_number,library_set)
if nargin<2
    library_set='ENGLISH';
end

%constructing language library
LIB.ENGLISH={'Zero','One','Two','Three','Four','Five','Six','Seven',...
    'Eight','Nine'};
LIB.SPANISH={'Cero','Uno','Dos','Tres','Cuatro','Cinco','Seis','Siete',...
    'Ocho','Nueve'};
LIB.JAPANESE={'Rei','Ichi','Nii','San','Shi','Go','Roku','Shichi',...
    'Hachi','Kyuu'};

%select library
Library=LIB.(upper(library_set));

num_IND=regexp(string_with_number,'\d');
num_list=unique(string_with_number(num_IND));
string_without_number=string_with_number;
for n = 1:length(num_list)
    string_without_number=strrep(string_without_number,num_list(n),...
        Library{str2num(num_list(n))+1});
end
end