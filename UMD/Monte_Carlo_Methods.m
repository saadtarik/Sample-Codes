%ENCE689E PS#3
%
%Author: Saad Tarik

clear;clc

%% Problem 1

%Define graphical parameters
nbins = 100; %number of bins
fontsize = 10;
%(a)
%Specify parameters for exponential distribution
mu = 1.0; %Given mean = 1.0
lambda = 1/mu; %parameter of the exponential distribution

%Specify variable
x = 1:0.2:10;

%Plot the PDF
figure (1);clf
plot(x,exppdf(x,lambda),'LineWidth',2)
xlabel('x','FontSize',fontsize);ylabel('p_X(x)','FontSize',fontsize)
title('PDF of Exponential Distribution')
%Save image to file
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','PS3_Prob_1a');close

%(c),(d) and (e)
%Simulation
N = [10 100 1000 10000]; %Sample sizes
mu_x1 = []; %initialize mean vector for Problem 1(c)
var_x1 = []; %initialize variance vector Problem 1(c)
mu_x2 = []; %initialize mean vector for Problem 1(e)
var_x2 = []; %initialize variance vector for Problem 1(e)

for i = 1:length(N)
    %Generate uniformly distributed random numbers (Pr(X<=x)) and calculate corresponding x
    y = 0 + (1-0)*rand(1,N(i));
    x = -(1/lambda)*(log(1-y));
    %Calculate sample statistics
    mu_x1 = [mu_x1 mean(x)];
    var_x1 = [var_x1 var(x)];
    %Plot histogram
    figure(2)
    subplot(2,2,i)
    [count,center] = hist(x);
    rel_freq = 100 * count / sum(count);
    bar(center,rel_freq);hold on
    xlabel('x','FontSize',fontsize);ylabel('Relative Frequency [%]','FontSize',fontsize)
    title_string = ['N = ',num2str(N(i))];
    title(title_string,'FontSize',fontsize)
end
%Save image to file
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','PS3_Prob_1d');close

for i = 1:length(N)
    %Generate exponentially distributed random numbers for variable x
    x = exprnd(mu,1,N(i));
    %Calculate sample statistics
    mu_x2 = [mu_x2 mean(x)];
    var_x2 = [var_x2 var(x)];
    %Plot histogram
    figure(3)
    subplot(2,2,i)
    [count,center] = hist(x);
    rel_freq = 100 * count / sum(count);
    bar(center,rel_freq)
    xlabel('x','FontSize',fontsize);ylabel('Relative Frequency [%]','FontSize',fontsize)
    title_string = ['N = ',num2str(N(i))];
    title(title_string,'FontSize',fontsize)
end
%Save image to file
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','PS3_Prob_1e');close 

%% Problem 2

%(b) Generate random vector of x, given w = ln(x) is normally distributed

n = 1000; % ensemble size
m_w = [-0.0431 0.6628]; % mean vector of w
C_ww = [0.0862 0;0 0.0606]; % covariance matrix of w

w = mvnrnd(m_w,C_ww,n); % nomrally distributed vector of w
x = exp(w); % lognormally distributed vector of x

%Calculate sample mean and covariance of x
m_x = mean(x,1); % mean vector of x
C_xx = cov(x(:,1),x(:,2)); % covariance matrix of x

%Plot the results
figure(4);clf
plot(x(:,1),x(:,2),'+k')
xlabel('x_1','FontSize',fontsize);ylabel('x_2','FontSize',fontsize)
%Save image to file
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','PS3_Prob_2c');close

%Generate ensemble of vector y
y = [1 + 0.25*x(:,1) - 0.1*x(:,2).^2 0.05*x(:,1) + 0.1*x(:,1).^2 + 0.1*x(:,2).^2] ;

%Calculate sample mean and covariance of y
m_y = mean(y,1); % mean vector of y
C_yy = cov(y(:,1),y(:,2)); % covariance matrix of y

%Plot the results
figure(5);clf
plot(y(:,1),y(:,2),'+k')
xlabel('y_1','FontSize',fontsize);ylabel('y_2','FontSize',fontsize)
%Save image to file
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','PS3_Prob_2e');close
