x = importdata('x.m');
y = importdata('y.m');
z = importdata('z.m');

f_x = fft(x);
f_y = fft(y);
f_z = fft(z);

psd_x = (f_x.*conj(f_x))/4140;
psd_y = (f_y.*conj(f_y))/4140;
psd_z = (f_z.*conj(f_z))/4140;

freq = 1/(0.000655*4138)*(0:4138);

figure;
plot(freq,psd_x);
hold on
plot(freq,psd_y);
plot(freq,psd_z);
