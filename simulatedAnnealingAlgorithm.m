function [bestsol,fmin] = Final(alpha)
clc
clear
close all
if nargin<1
    alpha = 0.95;% Cooling factor
end
% Initializing parameters and settings
T_init= 1; % Initial temperature
T = T_init; % Temp variable
T_min = 1e-10; % Minimum temperature
max_rej = 2500; % Max # of rejections
max_run = 500; % Max # of runs
max_accept = 250; % Max # of acceptances
s= randi(3,1,5); % Channels
i= 0;j=0; accept = 0;
G = rand (5,5);
G(logical(eye(size(G)))) = 0;
G = G./100;
E_init = fun(s,G);
E_old = E_init; count = 0;
best = s; fminimum = []; counter = [];
while ((T>T_min) || (j<=max_rej))
    i = i+1;
    %Check if max numbers of run/accept are met
    if(i>=max_run) || (accept >=max_accept)
        %reset the counters
        i=1; accept =1;
        %Cooling according to a cooling schedule
        T = cooling(alpha,T);
    end
    %new solution for channel assignment for all APs
    ns=randi(3,1,5);
    E_new = fun(ns,G);    
    % Decide to accept new solution
    %Accept if improved
    DeltaE = E_new - E_old;
    if(E_new < E_old) 
        best = ns; E_old = E_new;
        accept = accept + 1; j=0;
    end
    p= min(1,exp(-DeltaE/T)); % Probability of accepting new solution
    %Accept with a small probability p if not improved
    if ((E_new >= E_old) && (p > 0))
        best = ns; E_old = E_new;
        accept = accept + 1;
    else
        j=j+1;
    end
%Update the estimated optimal solution
fmin = E_old;
count = count + 1;
fminimum(count) = fmin; %#ok<AGROW>
counter(count) = count;  %#ok<AGROW>
%display(fmin)

end
plot (counter,fminimum)
xlabel ('Iteration');
ylabel ('Interference');
title ('Variation in interference')
display(count)
display(E_init)
bestsol = best;
bestfunctionvalue = fmin;
display(bestfunctionvalue)
end
function I=fun(s,G)
E_norm = 25; % Normalization
omega = zeros([5 5]);
if (s(1,1) == s(1,2))
    omega(1,2) = 1;
    omega(2,1) = 1;
end
if(s(1,3)==s(1,4))
    omega(3,4)=1;
    omega(4,3)=1;
end
rho = zeros([3 5]);
for y=1:5
    switch s(1,y)
        case 1
            rho(1,y)=1;
        case 2
            rho(2,y)=1;
        case 3
            rho(3,y)=1;
    end
end
noise = rand(1,5);
noise= sum(noise);
I=omega.*G + noise;
I(logical(eye(size(I)))) = 0;
I= sum(I);
I=sum(I,2)/E_norm;
end
function T=cooling(alpha,T)
T=alpha*T;
end
