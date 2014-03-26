function logical_vect=MakeLogical(L,IND)
%make a vector of logicals of length L and with 1s at IND
logical_vect=logical(zeros(1,L));
logical_vect(IND)=1;
end