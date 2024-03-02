clearvars; 
close all; 
clc;

m = 50; 
sigma = 0.5;
A = 0.6;
deltat = 0.5;
thipsiw = [(deltat^2)/2 deltat 1]';
thipsiz = [(deltat^2)/2 deltat 0]';
thiphi = [1 deltat (deltat^2)/2; 0 1 deltat; 0 0 A];
psiw = blkdiag(thipsiw, thipsiw);
psiz = blkdiag(thipsiz, thipsiz);
phi = blkdiag(thiphi, thiphi); 

P = (1/20)*[16,1,1,1,1;1,16,1,1,1;1,1,16,1,1;1,1,1,16,1;1,1,1,1,16];

x_0 = [normrnd(0, 500); normrnd(0, 5); normrnd(0, 5); normrnd(0, 200); normrnd(0, 5); normrnd(0, 5)];
Z = [0 0; 3.5 0; 0 3.5; 0 -3.5; -3.5 0]';

X = repmat(x_0, 1, m);
vec2 = randi(length(Z));

for i = 1:m-1
    W = normrnd(zeros(2,1),[sigma^2 sigma^2]');
    
    X(:,i+1) = phi*X(:,i) + psiz*Z(:,vec2) + psiw*W;
    
    vec2 = randsample([1 2 3 4 5], 1, true, P(vec2,:));    
end

% Plot results
%figure(1)
%plot(X(1,:), X(4,:))
%title('Path')
%xlabel('X^1')
%ylabel('X^2')

load stations.mat
load RSSI-measurements.mat

N = 10000;
m = length(Y); 
eta = 3;
C = 1.5;
v = 90;
tau = zeros(2,m);
ppzero = zeros(N,m);

func1 = @(x,X) mvnpdf(x,v-10*eta*log10(calculateDistance(X,pos_vec,N)),diag(ones(1,6)*C^2));

vec2 = zeros(N,m);
vec2(:,1) = randi(length(Z),[N 1]);

X = [mvnrnd(zeros(1,6),diag([500;5;5;200;5;5]), N) Z(:,vec2(:,1))']; 
ppzero(:,1) = func1(Y(:,1)',X);
tau(1,1) = sum(X(:,1).*ppzero(:,1))/sum(ppzero(:,1));
tau(2,1) = sum(X(:,4).*ppzero(:,1))/sum(ppzero(:,1));

[w, et, indices] = SEQSAMP(func1, ppzero,tau, X, Y, Z, vec2, sigma, psiw, psiz, phi);

% figure(2)
% plot(et(1,:), et(2,:))
% hold on 
% scatter(pos_vec(1,:), pos_vec(2,:),'filled')
% title('Expected path w/ stations')
% xlabel('x^1')
% ylabel('x^2')



% subplot(2,2,1)
% histogram(log10(w(:,1)))
% title('n = 1')
% 
% 
% 
% subplot(2,2,2)
% histogram(log10(w(:,20)))
% title('n = 20')
% 
% 
% subplot(2,2,3)
% histogram(log10(w(:,40)))
% title('n = 40')
% 
% 
% subplot(2,2,4)
% histogram(log10(w(:,60)))
% title('n = 60')


% for i = 1:80
%     CV = (1/N)*sum(((N * w(:,i)./sum(w(:,i))) - 1).^2);
%     ESC(i)=(10000)/(1+CV)
% end
% plot(ESC)


%resetting
tau2 = zeros(2,m);
omega2 = zeros(N,m);


vec1 = zeros(N,m);
vec1(:,1) = randi(length(Z),[N 1]);

X2 = [mvnrnd(zeros(1,6),diag([500;5;5;200;5;5]), N) Z(:,vec1(:,1))']; 
omega2(:,1) = func1(Y(:,1)',X2);
tau2(1,1) = sum(X2(:,1).*omega2(:,1))/sum(omega2(:,1));
tau2(2,1) = sum(X2(:,4).*omega2(:,1))/sum(omega2(:,1));

[w2, et2, indices2] = SISR(func1,omega2,tau2, X2, Y, Z, vec1, psiw, psiz, phi);
% figure(1)
% plot(et2(1,:), et2(2,:))
% hold on 
% scatter(pos_vec(1,:), pos_vec(2,:),'filled')
% title('Expected trajectory with SISR')
% xlabel('x^1')
% ylabel('x^2')

totw = zeros(m,5);
proba = zeros(5,m);
for i = 1:5
    totw(:,i) = [indices2 == i]'*omega2/sum(omega2);
end

% totw = totw';
% [maxi, tt] = max(totw);
% figure(2)
% histogram(tt)
% title('Highest prob direction')
%xticklabels({' ','None', ' ', 'East', ' ', 'North', ' ', 'South', ' ', 'West'})

%5

load RSSI-measurements-unknown-sigma.mat
% 
% 
% sigma = 0:0.5:3;
% 
% maxlike = zeros(1,length(sigma));
% 
% count = 1;
% 
%     taunew = zeros(2,m);
%     omeganew = zeros(N,m);
%       vec2 = zeros(N,m);
%     vec2(:,1) = randi(length(Z),[N 1]);
% for i = sigma
% 
% 
% 
%     X = [mvnrnd(zeros(1,6),diag([500;5;5;200;5;5]), N) Z(:,vec2(:,1))']; 
%     omeganew(:,1) = PP(Y(:,1)',X);
%     taunew(1,1) = sum(X(:,1).*omeganew(:,1))/sum(omeganew(:,1));
%     taunew(2,1) = sum(X(:,4).*omeganew(:,1))/sum(omeganew(:,1));
%     [wnew, u,uu] = SISR(PP,omeganew,taunew, X, Y, Z, vec2, psiw, psiz, phi);    
%         maxlike(count) = sum(log(sum(wnew)));
%         count = count + 1
% end


%scatter(sigma,maxlike,'filled')
%title('Maximum likelihood of C^2')

Cnew= 2.5;


tau = zeros(2,m);
omega = zeros(N,m);

blank = zeros(N,m);
blank(:,1) = randi(length(Z),[N 1]);

X_new = [mvnrnd(zeros(1,6),diag([500;5;5;200;5;5]), N) Z(:,blank(:,1))']; 
omega(:,1) = func1(Y(:,1)',X);
tau(1,1) = sum(X(:,1).*omega(:,1))/sum(omega(:,1));
tau(2,1) = sum(X(:,4).*omega(:,1))/sum(omega(:,1));
PPP = @(x,X_new) mvnpdf(x,v-10*eta*log10(calculateDistance(X_new,pos_vec,N)),diag(ones(1,6)*Cnew^2));

[w3, et3, indices3] = SISR(PPP,omega,tau, X_new, Y, Z, blank, psiw, psiz, phi);
figure(1)
hold on
plot(et3(1,:), et3(2,:))
scatter(pos_vec(1,:), pos_vec(2,:), 'filled')

title('Expected trajectory with C = 2.5')
xlabel('x^1')
ylabel('x^2')

function [w2, et2, indices2] = SISR(PP, og,tau, X, Y, Z, index, psiw, psiz, phi)
    samp = ones(length(X),1).*[1 2 3 4 5];
    [N,m] = size(og);
    
        z = X(:,7:8);
        [~,old] = ismember(z, Z', 'rows'); 
        ny = [ones(N,15).*old samp];

    for n = 2:m     

        index(:,n) = ny(sub2ind(size(ny),(1:N)',randi(20,N,1)));
        
        W = mvnrnd(zeros(2,1),[0.5^2 0.5^2],N);
        X = [[(phi*X(:,1:6)' + psiz*X(:,7:8)' + psiw*W')]' Z(:,index(:,n))'];
        
        oldold = randsample(N, N, true, og(:,n-1));
        X = X(oldold,:);
        og(:,n) = PP(Y(:,n)',X);
       
        tau(1,n+1) = sum(X(:,1).*og(:,n))/sum(og(:,n));


        tau(2,n+1) = sum(X(:,4).*og(:,n))/sum(og(:,n));
    end  
    
    w2 = og;
    et2 = tau;  
    indices2 = index;
end

function distancefunc = calculateDistance(inputMatrix, positionVector, N)
nn = length(positionVector);
distanceMatrix = zeros(nn,N);
for i = 1:nn
distanceMatrix(i,:) = vecnorm([inputMatrix(:,1) inputMatrix(:,4)]'-positionVector(:,i));
end
distancefunc = distanceMatrix';
end

function [w, et, indices] = SEQSAMP(PP, C,tau, X, Y, Z, indipp, sigma, psiw, psiz, phi)
    samp = repmat(1:5, length(X), 3);


    [N,m] = size(C);
            z = X(:,7:8);

        [u,old] = ismember(z, Z', 'rows'); 
        new = [ones(N,15).*old samp];
    
    for i = 2:m     

        indipp(:,i) = new(sub2ind(size(new),(1:N)',randi(20,N,1)));
       
        W = mvnrnd(zeros(2,1),[sigma^2 sigma^2],N);
        X = [[(phi*X(:,1:6)' + psiz*X(:,7:8)' + psiw*W')]' Z(:,indipp(:,i))'];
        C(:,i) = PP(Y(:,i)',X).*C(:,i-1);
 
        tau(1,i) = sum(X(:,1).*C(:,i))/sum(C(:,i));
        tau(2,i) = sum(X(:,4).*C(:,i))/sum(C(:,i));
    end  
    
    et = tau; 
    w = C; 
    indices = indipp;
end