function output=recur_function()
A=input('enter something','s');
if isempty(A)
    recur_function();
else
    output=A;
end
end