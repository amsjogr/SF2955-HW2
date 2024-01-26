clear vars
close all
clc

                            
%Germany
infected = load('germany_infected.csv');      
recovered = load('germany_removed.csv');       
Population = 83*10^6;                            
%Iran
infected = load('iran_infected.csv');         
recovered = load('iran_removed.csv');          
Population = 84*10^6; 
Sus = Population-infected-recovered;
td = length(infected);

nmbr = 4; %nmbr lambdas

lambdas = rand(1,nmbr);
t_int = round(linspace(0,td,nmbr+1));
tt = t_int(2:nmbr);
N = 12000;

sigma = 0.05;
M = 3;
a = 3;
b = 3;


bi = 3*ones(length(lambdas),1);




lambdanew = zeros(length(lambdas),N);


lambdanew(:,1) = lambdas;


tnew = zeros(length(tt),N);
tnew(:,1) = tt;


pirnew = zeros(1,N);
pirnew(1) = pirfunc(a,b,Sus,infected);
for k = 1:(12000-1)
    
    for i = 1:length(lambdas)
        pote = lambdas(i) + sigma*randn;
        if pote <= 0
           lambdanew(:,k+1) = lambdas; 
           continue
        end
        lambdastar = lambdas;
        lambdastar(i) = pote;
        
        alfa = min([0, lambdafunc(i,lambdastar,tt,Sus,infected,Population,bi) - lambdafunc(i,lambdas,tt,Sus,infected,Population,bi)]); % logarithm of alpha;
        
        if rand <= exp(alfa)
            lambdanew(:,k+1) = lambdastar;
            lambdas = lambdastar;
        else
            lambdanew(:,k+1) = lambdas;
        end
    end
    
    for i = 1:length(tt)
        pote = tt(i) + randi([-M M]);
        tstar = tt;
        tstar(i) = pote;
        if length(tstar) == 1
            rul = false;
        else
            rul = min(diff(tstar)) <= 0;
        end
        
        if rul || pote >= td || pote <= 0    
            tnew(:,k+1) = tt;
            continue
        end
        
        alfa = min([0, tfunc(i,tstar,lambdas,Sus,infected,Population) - tfunc(i,tt,lambdas,Sus,infected,Population)]);
        
        if rand <= exp(alfa)
            tnew(:,k+1) = tstar;
            tt = tstar;
        else
            tnew(:,k+1) = tt;
        end
    end
    
    pirnew(k+1) = pirfunc(a,b,Sus,infected);
    
end 




for i = 1:height(lambdanew)
    figure(1)
    ylabel('Lambdas')
    hold on
  plot(lambdanew(i,:),'DisplayName',num2str(i))
  legend
end


for i = 1:height(tnew)
    figure(2)
    ylabel('Breakpoint')
    hold on
  plot(tnew(i,:))
end






function func = lambdafunc(j,lambdas,tt,Sus,Infected,Population,bi)
phi = 0.995;
td = length(Infected); 
ai = 2*ones(length(lambdas),1);     
lambda = lambdas(j);

    if j == 1                   
        tid = 1:tt(j)-1;
    elseif j == length(lambdas)    
        tid = tt(end):td-1;
    else                        
        tid = tt(j-1):tt(j)-1;
    end
    ttt = Sus(tid,:) - Sus(tid+1,:);
    kappa = (1/phi -1)*Sus(tid).*(1-exp(-lambda*Infected(tid)/Population));
    pdt = gammaln(ttt+kappa) - gammaln(ttt+1) - gammaln(kappa) + kappa*log(1-phi) + ttt*log(phi);
    
    func = sum(ai.*log(bi) - gammaln(ai) + (ai-1)*log(lambda) - bi*lambda) + sum(pdt);
end

function func = tfunc(i,tt,lambdas,Sus,Infected,Population)
td = length(Infected); 
t_int = [1 tt td];
func = 0;
phi = 0.995;

for j=[i,i+1]  
    time = t_int(j):(t_int(j+1)-1);
    lambda = lambdas(j);

    dts = Sus(time,:) - Sus(time+1,:);
    cons = (1/phi -1)*Sus(time).*(1-exp(-lambda*Infected(time)/Population));
    P_dtI = gammaln(dts+cons) - gammaln(dts+1) - gammaln(cons) + cons*log(1-phi) + dts*log(phi);
    
    func = func + sum(P_dtI);
    
end
func = func + log(sum((min(diff(t_int))>0)));
end



function func = pirfunc(a, b, Sus, Infected)
    tinf = -diff(Infected);

    tsus = -diff(Sus);
    pir_alfa = sum(tsus + tinf) + a;
    pir_beta = sum(Infected(1:end-1) - tsus - tinf) + b;
    func = betarnd(pir_alfa, pir_beta);
end