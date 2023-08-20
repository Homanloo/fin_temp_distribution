clearvars
clc

%info required
T_inf = 1600;
T_base = 673;
T_all_max = 1350;
h = 800;
k = 20;
n = 51;
m = 7;
leng = 0.05;

%assuming delta X = delta Y = 1mm
delta = leng/(n-1);

%calculating biot number to use in
%Gauss-Seidel iterations
Bi = h * delta / k;

%large enough counter for all iterations

%calculations will stop exaclly when we
%reach error critia
counter = 100000;

%can change error critia in here
error_critia = 0.01;

%setting the middle node of the fin to test
% if the error critia has been accomplished
mid_m = round(m/2);
mid_n = round(n/2);

T_initial = (T_base + T_inf)/2;

T = ones(m, n, counter).*T_initial;

%"i" counts vertically
%"j" counts horizontally

for e = 1:counter
    for i = 1:m
        for j = 1:n

            %node in the joint to the base
            if j == 1
                T(i, j, e+1) = T_base;
            
            %insulated tip of the fin
            elseif j == n && i ~= 1 && i ~= m
                T(i, j, e+1) = (T(i+1, j, e) + T(i-1, j, e) + (2 * T(i, j-1, e)))/4;

            %upper nodes with convection
            elseif i == 1 && j ~= 1 && j ~= n
                T(i, j, e+1) = (T(i+1, j , e) + ((T(i, j+1, e) + T(i, j-1, e)) / 2) + Bi * T_inf) / (2 + Bi);
            
            %lower nodes with convection
            elseif i == m && j ~= 1 && j ~= n
                T(i, j, e+1) = (T(i-1, j , e) + ((T(i, j+1, e) + T(i, j-1, e)) / 2) + Bi * T_inf) / (2 + Bi);
            

            %middle nodes
            elseif i ~= 1 && i ~= m
                if j ~= 1 && j ~= n
                    T(i, j, e+1) = (T(i+1, j, e) + T(i-1, j, e) + T(i, j+1, e) + T(i, j-1, e)) / 4;

                end
            end

        end
    end

    
    if e > 10 && T(mid_m, mid_n, e) - T(mid_m, mid_n, e-1) < error_critia
        %printing number of required iterations
        fprintf('number of iterations: %d ', e)
        break
    end
end


x = 1:m;
y = 1:n;

for i = 1:m
    for j = 1:n
        z(i, j) = T(i, j, e);
    end
end

contourf(y, x, z, 'ShowText','on')
xlabel('Length (mm)')
ylabel('Thickness (mm)')
title('St-St conduction with insulated tip (part b)')