clear all
close all

B_campo = 10.0;

directory = dir('./res01032017/');

output = sprintf('./output-B_campo_%0.4f.dat', B_campo);

f = fopen(output, 'w');

for i = 3:size(directory, 1)/2+1
  display([num2str(i) ') ' directory(i).name]);

  name = directory(i).name;

  Rmax = str2num(name(1,52:57));

  file = sprintf('./res01032017/%s', directory(i).name);
  data = load(file);

  I = (data(:,1) - B_campo == 0);
  row = find(I);

  fprintf(f, '%0.4f', Rmax);
  for j = 2:size(data,2)
    fprintf(f, '   %0.15f', data(row,j));
  end
  fprintf(f, '\n');


end

fclose(f);
