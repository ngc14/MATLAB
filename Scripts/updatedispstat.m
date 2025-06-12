function updatedispstat(varargin)
txt = varargin{1};
a = varargin{2};
b = varargin{3};
persistent count total;
if(a==0)
    count = 0;
    total = b;
    dispstat('','init');
    fprintf(strcat("Iterating ", num2str(total), " instances.\n"),'keepthis');
end
fprintf('%s\n',strcat(txt,' ',num2str(count), " of ", num2str(total),": ",...
    num2str(100*count/total,'%.2f'),"% done."));
count = count + 1;
end