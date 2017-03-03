close all
IC = [0.015, -0.015,0,0]';
dofComp = 1;
tol = [1e-2, 5e-2, 1e-3, 5e-3, 1e-4, 5e-4, 1e-5, 5e-5, 1e-6, 5e-6, 1e-7, 5e-7, 1e-8, 5e-8, 1e-10];


n = length(tol);
x4 = cell(n,1);
x1 = cell(n,1);
t = cell(n,1);
error = cell(n,1);
errorf = zeros(n,1);
time = zeros(n,1);

tic
[t{n},x4{n}] = integration(M,K,zeros(2),a,IC,100*1/27.5664,0,0,0,tol(n),tol(n));
time(n) = toc;
x1{n} = x4{n}(:,dofComp);
error{n} = zeros(length(t{n}),1);
errorf(n) = 0;
for i = n-1:-1:1;
    tic 
     [t{i},x4{i}] = integration(M,K,zeros(2),a,IC,100*1/27.5664,0,0,0,tol(i),tol(i));
     time(i) = toc;
     x1{i} = x4{i}(:,dofComp);
     error{i} = interp1(t{n},x1{n},t{i},'spline') - x1{i};
     errorf(i) = abs(error{i}(end));
end

figure;
title('error for 1 dof in function time and tolerance')
xlabel('time [s]')
ylabel('log tolerance')
zlabel('error [m]')
hold on
for i=1:n
    y = log10(tol(i));
    ymat = y*ones(length(error{i}),1);
    
    plot3(t{i},ymat,error{i})
end
view(-31,36)
hold off

figure;
loglog(tol,errorf);
title('error final in function tolerance')
xlabel('tolerance')
ylabel('error [m]')

figure;
semilogx(tol,time)
title('duration integration in function tolerance')
xlabel('tolerance')
ylabel('time [s]')
  
figure;
plot3(log10(tol),log10(errorf),time);
title('duration in function tolerance and error')
xlabel('log tolerance')
ylabel('log error [m]')
zlabel('time [s]')
