clear all
close all

%% Make the covariance matrices of some random data
C1 = cov(randn(100,50));
C2 = cov(randn(100,50));
C3 = cov(randn(100,50));
C4 = cov(randn(100,50));
C5 = cov(randn(100,50));
C6 = cov(randn(100,50));

%% Make input for COVSCA
Cinput = [C1 C2 C4 C4 C5 C6];

%% Input parameters
nanal = 10;  %Number of analysis

%Rank and number of prototypes. 2 matrices of  this case, all of Rank 2
%of rank 1
Q= [2 2]';
L = length(Q);


%% Run COVSCA
[loadings, scores,fp,dys, func] = covsca(Cinput,L,Q,1,1,nanal);


% Fit percentages
disp(fp)

%% Plot scores and Loadings
figure(1)
plot(scores(:,1),scores(:,2),'k.','MarkerSize',22);
xlabel('1st COVSCA component','FontSize',13);
ylabel('2nd COVSCA component','FontSize',13);
title('COVSCA scores (weights)','FontSize',16);

%% Plot scores and Loadings
figure(2)
subplot(4,1,1)
bar(loadings(:,1)');
ylabel('Loadings 1st comp','FontSize',13);
title('COVSCA loadings','FontSize',16);
subplot(4,1,2)
bar(loadings(:,2)');
ylabel('Loadings 2nd comp','FontSize',13);
subplot(4,1,3)
bar(loadings(:,3)');
ylabel('Loadings 3rd comp','FontSize',13);
subplot(4,1,4)
bar(loadings(:,4)');
xlabel('Variables','FontSize',13);
ylabel('Loadings 4th comp','FontSize',13);
