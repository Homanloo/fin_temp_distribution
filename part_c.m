clearvars
clc

%info required
T_initial = 400;
T_all_max = 1350;
h = 800;
k = 20;
n = 51;
m = 7;
leng = 0.05;
alpha = 3.5e-6;

%assuming delta X = delta Y = 1mm
delta = leng/(n-1);
Bi = h * delta / k;


% SET TIME IN HERE
time_h = 2;
tf = 3600 * time_h;


%critical forier number is 0.5; for a better result
%we use lower forier number, like 0.05
Fo = 0.05;

dt = Fo * delta ^ 2 / alpha;

%number of time steps
tt = round(tf / dt + 1);


%matrix for tempretures
T = 400 .* ones(m, n, tt);

for p = 1:tt
    time_passed = ((tt-1) * dt) / 3600;
    T_inf = 400 + (1200 * ((exp(-time_passed)) * ((1 * time_passed) + exp(time_passed) - 1)));
    for i = 1:m
        for j = 1:n
            
            %base tempreture
            if j == 1
                T(i, j, p+1) = 400;
            
            %upper corner of the tip    
            elseif j == n && i == 1
                T(i, j, p+1) = 2 * Fo * (T(i+1, j, p) + T(i, j-1, p) - (2 * T(i, j, p)) + (2 * Bi * (T_inf - T(i, j, p)))) + T(i, j, p);

            %lower corner of the tip
            elseif j == n && i == m
                T(i, j, p+1) = 2 * Fo * (T(i-1, j, p) + T(i, j-1, p) - 2 * T(i, j, p) + 2 * Bi * (T_inf - T(i, j, p))) + T(i, j, p);

            %other nodes of the tip
            elseif j == n && i ~= 1 && i ~= m
                T(i, j, p+1) = Fo * ((2 * Bi * (T_inf - T(i, j, p))) + 2*T(i, j-1, p) + T(i+1, j, p) + T(i-1, j, p) - 4 * T(i, j, p)) + T(i, j, p);

            %upper nodes with convection
            elseif i == 1 && j ~= 1 && j ~= n
                T(i, j, p+1) = Fo * ((2 * Bi * (T_inf - T(i, j, p))) + 2*T(i+1, j, p) + T(i, j-1, p) + T(i, j+1, p) - 4 * T(i, j, p)) + T(i, j, p);

            %lower nodes with convection
            elseif i == m && j ~= 1 && j ~= n
                T(i, j, p+1) = Fo * ((2 * Bi * (T_inf - T(i, j, p))) + 2*T(i-1, j, p) + T(i, j-1, p) + T(i, j+1, p) - 4 * T(i, j, p)) + T(i, j, p);

            %inner nodes
            elseif i ~= 1 && i ~= m && j ~= 1 && j ~= n
                T(i, j, p+1) = Fo * (T(i+1, j, p) + T(i-1, j, p) + T(i, j+1, p) + T(i, j-1, p) - (4 * T(i, j, p))) + T(i, j, p);



            end
        end
    end
end


%check if error critia is considered or not
delta_T = T(round(m/2), round(n/2), tt+1) - T(round(m/2), round(n/2), tt);
fprintf('error critia is %d \n', delta_T)

%check tempretures above the allowable temp
for jj = 1:n
    if T(round(m/2), jj, tt+1) > T_all_max
        fprintf('there is higher tempretures than the allowable temp \n')
        break
    end
end


x = 1:m;
y = 1:n;

for i = 1:m
    for j = 1:n
        z(i, j) = T(i, j, tt+1);
    end
end

contourf(y, x, z, 'ShowText','on')
xlabel('Length (mm)')
ylabel('Thickness (mm)')
title('Trancient Conduction after 2 hour')





