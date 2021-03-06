x = linspace(1.93, 2.08, 151);
p = [1 -18 144 -672 2016 -4032 5376 -4608 2304 -512];
f_horner = Horner(p, x);
h = figure();
plot(x, f_horner, 'red');
hold on;
f = (x-2).^9;
plot(x, f);
legend('Horner Method', 'Direct Evaluation', 2);
title('Comparison of Horner and Direct Evaluation method'); 
xlabel('x'); 
ylabel('Value Calculated using Horner and Direct Evaluation'); 
saveas (h,'q4','jpeg');