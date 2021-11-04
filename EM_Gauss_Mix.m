%Load and plot "old faithful" dataset
raw_data = dlmread('old_faithful.txt', '\t', 26, 0);

x1 = raw_data(:,2); x2 = raw_data(:,3);
figure(1); clf
scatter(x1, x2); axis square; box on
xlabel('eruption duration (min)'); ylabel('time to next eruption (min)')

data = [x1 x2];

%%%% The EM algorithm %%%%

% Initialization
K = 2;
N = length(data);
D = 2;

old_means = zeros(K, D);
resp = zeros(N,K);
sum_gauss = zeros(N, 1);
SIGMAS = {[0.1 0;0 22.9],[0.1 0;0 22.9]}; %starting sigmas

%generate starting means
means = [(max(x1)-min(x1))*rand(K,1) + min(x1) (max(x2)-min(x2))*rand(K,1) + min(x2)];
PI = [0.5, 0.5];
V1 = zeros(2, N);	%storage variable
V2 = V1;

% Repeat until no change in means

while (sum(sum(abs(means - old_means))) > 1e-2)
	
% Assignment step
	%calculating responsibilities nominator part
	resp(:,1) = PI(1) * mvnpdf(data, means(1,:), SIGMAS{1});
	resp(:,2) = PI(2) * mvnpdf(data, means(2,:), SIGMAS{2});
	
	%adding sums of columns in resp vector into single column vector
	sum_gauss = (sum(resp'))';
	%div resp by total resp
	resp  = resp./sum_gauss;

% Update step
	old_means = means;
	means(1,:) = 1/sum(resp(:,1))*(sum(resp(:,1).*data));
	means(2,:) = 1/sum(resp(:,2))*(sum(resp(:,2).*data));
		
	PI(1) = 1/N*sum(resp(:,1));
	PI(2) = 1/N*sum(resp(:,2));
		
	for i = 1:N
		V1(:,i) = (resp(i, 1)*(data(i,:) - means(1,:)));
		V2(:,i) = (resp(i, 2)*(data(i,:) - means(2,:)));
	end
	%updating full covariance matrices
	SIGMAS{1} = V1*V1'./(sum(resp(:,1))); 	
	SIGMAS{2} = V2*V2'./(sum(resp(:,2))); 
end

MU1 = means(1,:);
MU2 = means(2,:);

p = @(xvals) PI(1)*mvnpdf(xvals,MU1,SIGMAS{1}) + PI(2)*mvnpdf(xvals,MU2,SIGMAS{2});

%Plot final means
figure(1); hold on
y1 = linspace(min(data(:,1)), max(data(:,1)));
y2 = linspace(min(data(:,2)), max(data(:,2)));
[Y1 Y2] = meshgrid(y1, y2);
hold on; contour(Y1, Y2, reshape(p([Y1(:) Y2(:)]), 100, 100), 5)
title('EM Algorithm'); xlabel('Eruption Time (min)'); ylabel('Waiting Time (min)'); colorbar; box on; axis square
