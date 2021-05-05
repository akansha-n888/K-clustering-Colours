clear
clc
close all

%% Lab 4 - K clustering
I= imread('house.tiff');    %load house image
figure,imshow(I);           %display image
M = length(I); N = length(I);   % assuming how is a symmetric image
X = reshape(I, M*N, 3);
X = double(X);
figure, plot3(X(:,1), X(:,2), X(:,3),'.','Color',[1, 0, 0])    %plot showing red
%% Part A - Randomly Initialize Paramters
disp("Part A")
n_samples = size(X,1);
c = 2;
u1 = [rand()*max(X(:,1)) rand()*max(X(:,2)) rand()*max(X(:,3))];
u2 = [rand()*max(X(:,1)) rand()*max(X(:,2)) rand()*max(X(:,3))];
u1_prev = [0 0 0];
u2_prev = [0 0 0];
initial = [u1;u2];

disp('Initial value of mean 1 for part A:');
disp(u1);
disp('Initial value of mean 2 for part A:');
disp(u2);

%% PART A - clustering 
cluster_id_a = zeros(n_samples,1);
attemptA = 0;
e = [];
Error_u1 = [];
Error_u2 = [];

while((~isequal(u1,u1_prev) && ~isequal(u2,u2_prev)))
    attemptA = attemptA + 1;
    for i = 1:n_samples
        norm_u1 = norm(X(i,:)-u1);
        norm_u2 = norm(X(i,:)-u2);
        
        if(norm_u1<norm_u2)
            cluster_id_a(i) = 1;
        else
            cluster_id_a(i) = 2;
        end
    end

    u1_prev = u1;
    u2_prev = u2;

    u1 = [mean(X(find(cluster_id_a==1),1)) mean(X(find(cluster_id_a==1),2)) mean(X(find(cluster_id_a==1),3))];
    u2 = [mean(X(find(cluster_id_a==2),1)) mean(X(find(cluster_id_a==2),2)) mean(X(find(cluster_id_a==2),3))];
    
    Error_u1 = [Error_u1 ;norm_u1^2];
    Error_u2 = [Error_u2 ;norm_u2^2];
    error = (Error_u1 + Error_u2);
    error = [e; error ];
    
end

finalA = [u1;u2];
disp("Final value of updated mean for part A")
disp(finalA);

disp("Number of Iterations");
disp(attemptA);

colour_u1 = [u1(1)/255 u1(2)/255 u1(3)/255];
colour_u2 = [u2(1)/255 u2(2)/255 u2(3)/255];

%% plot of Cluster Means (two stages of cluster process)
figure
plot3(X(find(cluster_id_a==1),1), X(find(cluster_id_a==1),2), X(find(cluster_id_a==1),3),'.','Color',colour_u1)
hold on
plot3(X(find(cluster_id_a==2),1), X(find(cluster_id_a==2),2), X(find(cluster_id_a==2),3),'.','Color',colour_u2)
hold off

%% plot of error
plot(1:attemptA,Error_u2)
title("Error Plot");
xlabel("Iteration");
ylabel("Error");
legend('Total Error')
hold off
%% PART B - Clustering
disp("Part B")
n_samples = size(X,1);
c = 5;
u = [];
attempt = 0;
for i = 1:c
    u = [u;rand()*max(X(:,1)) rand()*max(X(:,2)) rand()*max(X(:,3))];

end

u_prev = zeros(5,3);

initial_b = u;

cluster_id_b = zeros(n_samples,1);

while(~isequal(u,u_prev))
    disp(attempt)
    attempt = attempt + 1;
    for i = 1:n_samples
        for j = 1:c
            norm_u(j) = norm(X(i,:)-u(j,:));
        end
        [m,id] = min(norm_u);
        cluster_id_b(i) = id;
    end
    u_prev = u;

    for i = 1:c
       u(i,:) = [mean(X(find(cluster_id_b==i),1)) mean(X(find(cluster_id_b==i),2)) mean(X(find(cluster_id_b==i),3))];
    end
end

for i = 1:c
    colour_u(i,:) = [u(i,1)/255 u(i,2)/255 u(i,3)/255];
end

figure
hold on
for i = 1:c
    plot3(X(find(cluster_id_b==i),1), X(find(cluster_id_b==i),2), X(find(cluster_id_b==i),3),'.','Color',colour_u(i,:))
end
hold off

disp("Initial value of updated mean for part B")
disp(initial_b)

disp("Final value of updated mean for part B")
disp(u)

%% Part C - Xie-Beni Index
% compare between part a and b
% part b should be better
disp("Part C");

k = 0;
min_sep_a = norm(u1-u2);
XB_a = 0;
u_a = [u1;u2];
for k = 1:n_samples
    for j = 1:2
        if ~(cluster_id_a(k) == j)
            continue
        end
        
        XB_a = XB_a + norm(X(k,:)-u_a(j,:))/min_sep_a;
    end
    

end

XB_a = XB_a/n_samples;


XB = 0;
sep = [];

for i = 1:c
    for j = 1:c
        if i == j 
            continue 
        end
        
        sep = [sep;norm(u(j,:)-u(i,:))];
    end
end

min_sep= min(sep);


for k = 1:n_samples
    for j = 1:c
        if ~(cluster_id_b(k) == j)
            continue
        end
        
        XB = XB + norm(X(k,:)-u(j,:))/min_sep;
    end
end

XB = XB/n_samples;
disp("XB is:")
disp(XB)