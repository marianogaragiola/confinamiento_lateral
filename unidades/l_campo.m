clear all
close all

alpha = 658.4092645439;

b_campo = linspace(1, 50, 100);

l = sqrt(2.0*alpha./b_campo);

salida = [b_campo', l'];

save('-ascii', 'l_campo.dat', 'salida');
