function s=threedigitsH(t)

if t<10
    s=['00', num2str(t)];
elseif t<100
    s=['0', num2str(t)];
elseif t<1000
    s=[num2str(t)];
else
    fprintf('Error: More than three digits!');
end
return