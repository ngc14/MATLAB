function s=twodigitsH(t)

if t<10
    s=['0', num2str(t)];
elseif t<100
    s=[num2str(t)];
else
    fprintf('Error: More than two digits!');
end
return