clc;
clear all;
%%%%%%%%%%%%************Model***************%%%%%%%%%%%  
R_data = 12000000;          %(12 Mbps) replace this for 54Mbps evaluation
R_phy  = 6000000;          %(PHY rate)

%Model Parameters (in MicroSec)
SIFS = 0.00001; %10e-06; 
slot_time = 0.00002; %20e-06; 
DIFS = 0.00005; %50e-06;
delta= 0.000001; %1e-06;                %propagation delay(s)
payload = 12000;            %bits (1500B)
mac_header = 272;           %bits (34B)
phy_header = 128;           %bits (16B)
ack_p = 112;                %bits (14B)
packet = payload + mac_header; %+ phy_header;
ACK = ack_p + phy_header;

%*********Data packet and ACK time
T_data = (packet/R_data)+ (phy_header/R_phy);
T_ack = (ACK/R_data)+((phy_header/R_phy));          %includes ack and header pckt

%delay times (Basic access mechanism)
T_s = DIFS + T_data + 2*delta + SIFS + T_ack;
T_c = DIFS + T_data + delta;
T_ce = T_data + T_ack;
%*******************************************************
%Variable Number of Stations and Payload Size
N=0:8:200;
C=0.2*N;           
N(1)=1; 
C(1)=1; 
%parameters
W = 32;
m = 6;
Rtx=slot_time*(4*((2^0)*((W-1)/2)+(2^1)*((W-1)/2)+(2^2)*((W-1)/2)+(2^3)*((W-1)/2))); %Packets are sent in the 3rd tTx attempt assumption
E_s = delta + T_s+T_c + Rtx; %(crude approximation of the E[slot time])
lamda= 4;       %Data Traffic: 4 packets per second
mu = 1/E_s;
Q=lamda/mu;
%**************************************************************************
%Collision prob., Throughput and MAC delay calculated for different Num. of Users 

%*****Collision Probabilities*****#
P_fp=[];                    %Makama

%***Throughput********
S = [];                     %Done

 S_TrfMdl=[];               %Makama Simulation
 %r=[];
%******MAC Delay*********
Dmac = [];                  


Dmac_sim = [];              %Makama Simulation
%******Packet drop time/Delay*********
%Dmac_drp = [];                  

%******************************************************
%iterate for every possible number of Geophone
%for i=1:size(N, 2)
for k = 1:size(N, 2)        %Done
    if (N(k) == 1)
        p = 0.0;
        tau = 2/(W + 1);
      else
        inc = 0.000001;
        error_MAX = 0.000005;
        p= 0.0;
        tau1 = 0.0;
        tau2 = tau1 + 2*error_MAX;

        while abs(tau1 - tau2) > error_MAX
            p = p + inc;
            tau1 = 2*Q*P_cd*(1-2*p)*(1-(p).^(m+1))./(Q*((W*(1-p)*(1-(2*p).^(m+1))+ (1-2*p)*(2*P_cd-1)*(1-(p).^(m+1))))+2*(1-Q)*P_cd*(1-p)*(1-2*p));
            tau2 = 1 - (1-p).^(1./(N-1));
        end
        tau = (tau1 + tau2)./2;          %avg value of both tau
        
    end
%*******************Probabilities***************

    %Probability of collision
    P_c = 1- (1-tau).^(N-1);
    %Probability of failure due to channel error
    P_ce = 0.02;          % (2%)
    %Probability of packet failure (P_f)
    P_f = P_c+P_ce-(P_c*P_ce);
    %Packet drop probability
    P_D = P_f.^(m+1);
    %Probability of Tx (probability that there is at least one transmission in a the considered time slot)
    P_tx = 1 - (1- tau).^N(1,k);
    %Probability of successful Tx (prob. that exactly one geophone Txmits)
    P_TxG = N(1,k).*tau.*(1- tau).^(N(1,k)-1);
    %Probability that a packet transmission ongoing on the channel is successful
    P_suc=P_TxG./P_tx;
    %P_suc=(N(1,k)*tau*(1- tau)^(N(1,k)-1))./(1 - (1- tau)^N(1,k));
    %Probability of empty slot time
    P_e = 1 - P_tx;
    %Probability of collision
    %P_c = P_tr - P_s;
    %Probability counter decrements
    P_cd = (1-tau).^(N(1, k)-1);             
    
%****************************************
%Average slot duration
%#####################################
% CW_avgS over transmission all stages
 for i=1:m
    for j= 1:i
     if j<m
   CW_avgS = (((1-P_f).*(P_f)./(1-P_D)).*((2^j)*(W-1)));
     else
       CW_avgS =(((1-P_f).*(P_f)./(1-P_D)).*((2^m)*(W-1)));  
     end
    end
 end
%###############################################
%**********************Probabilities (pci,pcs and pcc)*************************

pci= (1-(1./CW_avgS)).^C;
pcs= C.*(1./CW_avgS).*(1-(1./CW_avgS)).^(C-1);
pcc= 1-pci-pcs;

%**********************Probabilities (pbi,pbs and pbc)*************************
pbi=(1-tau).^(N-1);
bs=nchoosek(199,1);             %199 = N-1          %**fix the problem wit N-1 here
pbs= bs.*tau.*(1-tau).^(N-2);
pbc=1-pbi-pbs;

%**********************Probabilities (pss and psi)*************************
pss=1/W;                   %1/CWmin
psi=1-pss;
%**********************Durations sigmaI, sigmaS and sigmaC*************************
sigmaI= slot_time;
sigmaS= (1/(1-pss))*T_s + (psi/(1-pss))*sigmaI;

for i=1:m
    for j=1:i
        if j<m
        sumC = j.*(pcc.^j);
        else
        sumC = m.*(pcc.^m);
        end
    end
end  

sigmaC= sumC*T_c + (pcs./(1-pcc))*sigmaS + (pci./(1-pcc))*sigmaI;

%**Durations for entering backoff stages from previous bckoff and Tx stage deltaB and deltaT***
 %Probability that the counter decrements
Pcd=(1-tau).^(N-1);     %Pcd=pbi
del=sigmaI.*pbi+sigmaS.*pbs+sigmaC.*pbc;

deltaB=(1./Pcd).*(del);       

deltaT=Q*((1-(1./CW_avgS))).*del;
%********Durations a geophone spends in a backoff stage (delta)********
D= (1-tau).*deltaB + tau.*deltaT;
%#Working (Showing curve behavior)
dmac_t2=0.0;
for i=1:m
    %avg_CW=167.5;
    
    for j= 1:i
        if j<m
             avg_CW = ((2^j)*(W-1))/2;
             dmac_t3=D.*(avg_CW*slot_time);
             dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^i).*(T_s+ i.*T_c+ dmac_t3)); %sum over i term of equation 40 in paper
             
        else
            avg_CW = ((2^m)*(W-1))/2;
            dmac_t3=D.*(avg_CW*slot_time);
            dmac_t2=(1./(1-P_D)).*(dmac_t2 + (P_suc.^m).*(T_s+ m.*T_c + dmac_t3)); 
        end
        
    end
    
 end


%*************************************
Dmac = dmac_t2;              %Makama 
r=rand;
S = (P_suc.*P_tx.*packet)./Dmac; 
 S_TrfMdl = (((packet).*20).*N)./(5.0+((Dmac).*20)); %Througput according to WSDA traffic model
 q_prob_12= lamda./(1./Dmac);
 %q_prob_12=lamda./Dmac;
% S_TrfMdl = ((20*(packet).*N)./((4.75)+(Dmac.*20))); 
end
%################Plots#################
%%%%%%%%%%%%%12 Mbps From Omnet simulator
 Thr_sim_12= [50526.31579 346559.7772	677232.9294	1013293.494	1348215.6	1681574.87	2015603.242	2350901.091	2686116.317	3021803.279	3355558.937	3690966.507	4025741.811	4360951.209	4693850.676	5029016.977	5364288.39	5698813.891	6033716.279	6362119.75	6695256.253	7030020.187	7364785.378	7096755.56	6528505.829	6435368.303
]; %throughput from simulation
err_thr_12= [5.48645E-12 10506.02327	4269.603179	5233.50614	5460.1618	3929.684887	4104.245662	4846.631886	5473.162887	6122.497077	6474.946105	7078.807353	7612.324778	8306.373826	8595.238641	9157.420181	9769.427893	10453.24257	11117.69836	9248.70781	7034.521138	7386.020488	7736.986717	388023.4144	76925.3259	23747.3509
]; %error from simulation    
OL=[50526.31579 404210.5263	808421.0526	1212631.579	1616842.105	2021052.632	2425263.158	2829473.684	3233684.211	3637894.737	4042105.263	4446315.789	4850526.316	5254736.842	5658947.368	6063157.895	6467368.421	6871578.947	7275789.474	7680000	8084210.526	8488421.053	8892631.579	9296842.105	9701052.632	10105263.16
];
%%%%%%%%%%%%%54 Mbps
Thr_sim_54=[50526.31579 346559.7772	677242.7216	1013304.877	1348236.496	1681598.871	2015634.902	2350934.692	2686158.531	3021859.348	3355609.058	3691014.652	4025741.811	4360951.209	4693850.676	5029073.787	5364345.606	5698872.885	6033790.375	6362185.533	6695354.891	7030122.656	7364890.369	7698417.559	8033131.281	8367083.454
]; %throughput from simulation for 54Mbps
err_thr_54=[5.48645E-12 10506.02327	4271.679503	5227.245297	5451.643617	3925.490118	4107.878068	4852.857182	5483.956931	6137.418485	6497.778885	7101.621953	7612.324778	8306.373826	8595.238641	9183.696341	9795.942744	10482.71473	11155.62902	9227.628351	7003.115637	7353.434962	7703.683135	7887.303157	8230.175436	8637.158945
]; %error from simulation for 54Mbps

%throughput from analytical model for 54Mbps
S_TrfMdl_54=[49088,389160.144075096,776145.195885316,1163062.70837431,1549956.89905444,1936820.56168648,2323623.00069602,2710306.58817647,3096779.58620775,3482907.91616028,3868505.52275812,4253321.73440562,4637023.04650988,5019165.68353143,5399153.78746408,5776175.67870939,6149106.69932434,6516360.74683253,6875662.36774493,7223695.27518230,7555559.07434883,7863932.38857529,8137801.43462799,8360590.92207193,8507603.18115206,8543026.12421342];
q_prob_54= [1.98156443161368e-60,0.00910642052856314,0.0119343702232391,0.0129393639030254,0.0134578587045112,0.0137851894190273,0.0140302447058821,0.0142498313629946,0.0144835667323613,0.0147658465505507,0.0151310310654923,0.0156165627107588,0.0162658138063144,0.0171315955459896,0.0182810522576872,0.0198027774176302,0.0218173642507136,0.0244933114307784,0.0280714238035472,0.0329029279009419,0.0395101041119045,0.0486844993709401,0.0616489071897514,0.0803293791297751,0.107820357780689,0.149194659744053];
%####################################
 %******* Sim vs Analytical **********
figure();
subplot(2,1,1);
plot(N, OL, 'kp-','LineWidth', 2,'MarkerSize',8);hold on;
plot(N, S_TrfMdl, 's--','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);hold on;
errorbar(N,Thr_sim_12,err_thr_12, 'o--', 'color', [0.8500 0.3250 0.0980],'LineWidth', 2,'MarkerSize',8);hold on;
xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput (12 Mbps)','Simulated throughput (12 Mbps)');
%title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput (12 Mbps)','Simulated throughput (12 Mbps)');
grid on;  
subplot(2,1,2);
plot(N, OL, 'kp-','LineWidth', 2,'MarkerSize',8);hold on;
plot(N, S_TrfMdl_54, 's--','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);hold on;
errorbar(N,Thr_sim_54,err_thr_54, 'o--','color', [0.8500 0.3250 0.0980],'LineWidth', 2,'MarkerSize',8);
xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput (54 Mbps)','Simulated throughput (54 Mbps)');
%title('Analytical Vs Simulation Throughput');xlabel('Number of geophones');ylabel('Throughput (bps)'); legend('Theoritical Offered Load','Analytical throughput (54 Mbps)','Simulated throughput (54 Mbps)');
grid on; 
%plot(N, S_TrfMdl, 'ks--');
%**********Collision probability Plots********
figure();
subplot(2,1,1);
plot(N, P_f(:,:,end),'s-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 1.00]);
title('Packet transmission failure probability');xlabel('Number of geophones');ylabel('Probability'); 
%######################################
grid on;  
subplot(2,1,2);
%figure();
plot(N, q_prob_12(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 0.3]);hold on;
plot(N, q_prob_54(:,:,end), 's-','color', [0.8500 0.3250 0.0980],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 0.3]);hold on;
xlabel('Number of geophones');ylabel('Probability q'); legend('12 Mbps','54 Mbps');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delay Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%N1=0:8:176;
N1=0:8:200;
N2=0:8:200;
N1(1)=1; 
C1(1)=1; 
for k1 = 1:size(N1, 2) 
    if (N1(k1) == 1)
        p1 = 0.0;
        tau_d = 2/(W + 1);
        
      else
        inc = 0.000001;
        error_MAX = 0.000005;
        p1= 0.0;
        tau1_d = 0.0;
        tau2_d = tau1_d + 2*error_MAX;

        while abs(tau1_d - tau2_d) > error_MAX
            p1 = p1 + inc;
%             tau1 = 2*Q*(1-2*p)*(1-(p).^(m+1))./(Q*((W*(1-p)*(1-(2*p).^(m+1))+ (1-2*p)*(1-(p).^(m+1))))+2*(1-Q)*(1-p)*(1-2*p));
            tau1_d = 2*Q*P_cd1*(1-2*p1)*(1-(p1).^(m+1))./(Q*((W*(1-p1)*(1-(2*p1).^(m+1))+ (1-2*p1)*(2*P_cd1-1)*(1-(p1).^(m+1))))+2*(1-Q)*P_cd1*(1-p1)*(1-2*p1));
            tau2_d = 1 - (1-p1).^(1./(N1-1));
            
        end
        tau_d = (tau1_d + tau2_d)./2;          %avg value of both tau
        
    end
%*******************Probabilities***************

    %Probability of collision
    P_c1 = 1- (1-tau_d).^(N1-1);
    %Probability of failure due to channel error
    P_ce1 = 0.0000000000000001;          % (2%)
    %Probability of packet failure (P_f)
    P_f1 = P_c1+P_ce1-(P_c1*P_ce1);
    %Packet drop probability
    P_D1 = P_f1.^(m+1);
    %Probability of Tx (probability that there is at least one transmission in a the considered time slot)
    P_tx1 = 1 - (1- tau_d).^N1(1,k1);
    %Probability of successful Tx (prob. that exactly one geophone Txmits)
    P_TxG1 = N1(1,k1).*tau_d.*(1- tau_d).^(N1(1,k1)-1);
    %Probability that a packet transmission ongoing on the channel is successful
    % P_suc1=1-P_f1;
    P_suc1=P_TxG1./P_tx1;
    %P_suc=(N(1,k)*tau*(1- tau)^(N(1,k)-1))./(1 - (1- tau)^N(1,k));
    %Probability of empty slot time
    P_e1 = 1 - P_tx1;
    %Probability of collision
    %P_c = P_tr - P_s;
    %Probability counter decrements
    P_cd1 = (1-tau_d).^(N1(1, k1)-1);             
    
  %  E_s = P_e1.*slot_time + (P_suc1.*P_tx1).*T_s + (P_tx1.*(1-P_suc1)).*T_c;
%****************************************
%CW_avgS1 = 0.0;
 for i=1:m
    for j= 1:i
     if j<m
   CW_avgS1 = (((1-P_f1).*(P_f1)./(1-P_D1)).*((2^j)*(W-1)));
     else
       CW_avgS1 =(((1-P_f1).*(P_f1)./(1-P_D1)).*((2^m)*(W-1)));  
     end
    end
 end

%CW_avgS1 = 167.5;        %calculated manually
%**********************Probabilities (pci,pcs and pcc)*************************

pci1= (1-(1./CW_avgS1)).^C1;
pcs1= C1.*(1./CW_avgS1).*(1-(1./CW_avgS1)).^(C1-1);
pcc1= 1-pci1-pcs1;

%**********************Probabilities (pbi,pbs and pbc)*************************
pbi1=(1-tau_d).^(N1-1);
bs1=nchoosek(175,1);             %199 = N-1          %**fix the problem wit N-1 here
pbs1= bs1.*tau_d.*(1-tau_d).^(N1-2);
pbc1=1-pbi1-pbs1;

%**********************Probabilities (pss and psi)*************************
pss1=1/W;                   %1/CWmin
psi1=1-pss1;
%**********************Durations sigmaI, sigmaS and sigmaC*************************
%sigmaI= 1; %*slot_time;
sigmaI= slot_time;
sigmaS1= (1/(1-pss1))*T_s + (psi1/(1-pss1))*sigmaI;

%sumC1=0.0;
for i=1:m
    for j=1:i
        if j<m
        sumC1 = j.*(pcc1.^j);
        %sumC1 = i.*pcc1.^i;
        else
         sumC1 = m.*(pcc1.^m);
        %sumC1 = m.*pcc1.^m;
        end
    end
end  

 sigmaC1= sumC1*T_c + (pcs1./(1-pcc1))*sigmaS1 + (pci1./(1-pcc1))*sigmaI;

%**Durations for entering backoff stages from previous bckoff and Tx stage deltaB and deltaT***
 %Probability that the counter decrements
 Pcd1 = (1-tau_d).^(N1-1); 
%Pcd1=(1-tau_d).^(N1-1);     %Pcd=pbi
del1=sigmaI.*pbi1+sigmaS1.*pbs1+sigmaC1.*pbc1;


deltaB1=(1./Pcd1).*(del1);       

deltaT1=Q*((1-(1./CW_avgS1))).*del1;
%********Durations a geophone spends in a backoff stage (delta)********
D1= (1-tau_d).*deltaB1 + tau_d.*deltaT1;

%************************************
dmac_t21=0.0;
% avg_CW1 = 167.5;
% dmac_t31=D1*(avg_CW1*slot_time);
for i=1:1:m
    for j= 1:i
        %i=i+1
        if j<m
             avg_CW1 = ((2^j)*(W-1))/2;
             dmac_t31=D1*(avg_CW1*slot_time);
             dmac_t21=(1./(1-P_D1)).*((dmac_t21+((P_suc1.^j).*(T_s+ j.*T_c + dmac_t31)))); 
              %sum over i term of equation 40 in paper
        else
            avg_CW1 = ((2^m)*(W-1))/2;
            dmac_t31=D1.*(avg_CW1*slot_time);
            dmac_t21=(1./(1-P_D1)).*((dmac_t21+((P_suc1.^m).*(T_s+ m.*T_c + dmac_t31)))); 
%             dmac_t31=D1.*(avg_CW1*slot_time);
%             dmac_t21=(1./(1-P_D1)).*((P_suc1.^m).*(T_s+ m.*T_c + dmac_t31));
        end    
        
    end

 end

%************************************

Dmac1 = dmac_t21;              %Makama
Dmac_sim = ((Dmac1).*20)./(N1(1,k1));
%Dmac_sim = (Dmac1.*20)./N1(1,k1);
end

%####### PLOTs #################
%12 Mbps
%Dmac_sim_12_N200=[1.01625707286532e-38,0.000797960999389173,0.00104516483022206,0.00113262307010078,0.00117736792010706,0.00120497993676577,0.00122467727136979,0.00124110461747166,0.00125740082438030,0.00127624536421375,0.00130029167614999,0.00133239840054321,0.00137581527487062,0.00143440208697078,0.00151293522795474,0.00161755709728334,0.00175644270458846,0.00194079800668918,0.00218637537524804,0.00251581224939738,0.00296230498437642,0.00357548616138383,0.00443099979174240,0.00564638760637528,0.00740793134523236,0.0100168461110371];
Dmac_sim_12=[1.70499644229892e-31,0.000992521352729925,0.00130094350517346,0.00141100591930581,0.00146762205470548,0.00150292801163491,0.00152878961908168,0.00155148440891787,0.00157553907691578,0.00160504582434040,0.00164423178368627,0.00169781493950418,0.00177136291996024,0.00187177809675687,0.00200801684193743,0.00219218042737664,0.00244119078369092,0.00277940583132661,0.00324277977182226,0.00388561810787892,0.00479177680714593,0.00609362048376930,0.00800479080836595, 0.00864638760637528,0.00940793134523236,0.0100168461110371];%0.00564638760637528,0.00740793134523236,0.0100168461110371];
MAC_D_sim= [0.00	0.001161475	0.001200727	0.001245023	0.001267362	0.001300007	0.001351038	0.00139213	0.001454977	0.001494873	0.00152785	0.001584056	0.001670921	0.001791917	0.001870966	0.002113992	0.002241747	0.002564038	0.002815672	0.003061233	0.003653249	0.004916805	0.006269095	0.278585242	0.652331463	0.835003954];
err_MAC_D= [0.00	7.07595E-05	4.45748E-05	6.15561E-05	5.91491E-05	5.27806E-05	6.35494E-05	4.88364E-05	8.72931E-05	8.35872E-05	8.97175E-05	0.000106468	0.000108018	0.000103798	9.77835E-05	0.000170373	0.000170857	0.000300435	0.000342006	0.000406274	0.000538557	0.00088955	0.00138065	0.175426983	0.042389356	0.023986584];
%N2=[184,192,200];

%******54 Mbps
Dmac_sim_54 =[4.93319144259841e-62,0.000227638449293516,0.000298332455839981,0.000323455752275035,0.000336417373189102,0.000344600194942903,0.000350726256816008,0.000356215636379443,0.000362058701803115,0.000369115288883930,0.000378244345900255,0.000390381870758273,0.000406612099237482,0.000428255224378705,0.000456989745575851,0.000495030358187006,0.000545391702727839,0.000612285983760474,0.000701732952107983,0.000822512737288617,0.000987681564548440,0.00121702693224098,0.00154111688438443,0.00200809991816095,0.00269533266643549,0.00372962836593021]; %Analytical
MAC_D_sim_54= [0.000	0.000311657	0.000315745	0.000322432	0.000325114	0.000327515	0.000331204	0.00033532	0.000341607	0.000346181	0.00034733	0.000349679	0.00035536	0.000359303	0.0003645	0.000373867	0.000377744	0.000384968	0.000389881	0.000395269	0.00040129	0.00040691	0.000411299	0.000420276	0.00042498	0.000431316];
err_MAC_D_54= [0.00	9.77124E-06	6.5672E-06	1.04476E-05	1.07374E-05	8.15356E-06	9.89266E-06	8.83824E-06	1.42677E-05	1.21098E-05	1.41912E-05	1.3998E-05	1.42944E-05	1.20029E-05	1.00803E-05	1.10767E-05	1.33059E-05	1.39478E-05	1.28095E-05	1.71911E-05	1.70533E-05	1.25725E-05	1.34697E-05	1.50259E-05	1.32219E-05	1.30039E-05];

Dmac_sim_54_ms=Dmac_sim_54.*1000;
MAC_D_sim_54_ms=MAC_D_sim_54.*1000;
err_MAC_D_54_ms=err_MAC_D_54.*1000;
%N2=[184,192,200];


%***************************************************
%******* Plot with collision prob. sub-plot *********************
figure();
subplot(2,1,1);
%plot(N1, Dmac_sim(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 1.0]);hold on;
plot(N1, Dmac_sim_12(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 0.15]);hold on;
errorbar(N1,MAC_D_sim,err_MAC_D, 'o-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 2,'MarkerSize',8);axis([0 200 0.00 0.900]);
xlabel('Number of geophones');ylabel('Delay (s)'); legend('Analytical 12Mbps','Simulation 12Mbps');

%%% Zoom plot new
% axes('position', [0.25 0.25 0.5 0.5]);
% box on
% my_index= 180<=N1 & N1<=200;  %Position to zoom on x-axis
% plot(N1(my_index),Dmac_sim_12(my_index),'s-','color', [0, 0.4470, 0.7410],'LineWidth', 1.5,'MarkerSize',10);axis([180 200 0.2 1.0]);hold on    %Analytical plot
% errorbar(N1(my_index),MAC_D_sim(my_index),err_MAC_D(my_index),'o-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 1.5,'MarkerSize',7);axis([180 200 0.2 1.0]);  %Simulation plot
% axis tight

%title('Delay: Analytical Vs Simulation')
% Zoom plot
axes('position', [0.25 0.25 0.5 0.5]);
box on
my_index= 0.0<N1 & N1<=181;  %Position to zoom on x-axis
%my_index_err= 0.0 <= N1 & N1<180;
plot(N1(my_index),Dmac_sim_12(my_index),'s-','color', [0, 0.4470, 0.7410],'LineWidth', 1.5,'MarkerSize',10);axis([1 180 0.00 0.008]);hold on    %Analytical plot
errorbar(N1(my_index),MAC_D_sim(my_index),err_MAC_D(my_index),'o-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 1.5,'MarkerSize',7);  %Simulation plot
axis tight
%title('Delay: Analytical Vs Simulation');
ax = gca;
ax.YAxis.Exponent = -3;

grid on; 
subplot(2,1,2);
%%%In milli seconds
plot(N2, Dmac_sim_54_ms(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 9.00]);hold on;%axis([0 180 0.00 0.014]);hold on;
%plot(N2, Dmac_sim_54(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 9.00]);hold on;
errorbar(N2,MAC_D_sim_54_ms,err_MAC_D_54_ms, 'o-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 2,'MarkerSize',8);
xlabel('Number of geophones');ylabel('Delay (ms)'); legend('Analytical 54Mbps','Simulation 54Mbps');
%title('Delay: Analytical Vs Simulation');
%ax = gca;
%ax.YAxis.Exponent = -3;
grid on; 
%%%In seconds
% plot(N2, Dmac_sim_54(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 0.009]);hold on;%axis([0 180 0.00 0.014]);hold on;
% %plot(N2, Dmac_sim_54(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 2,'MarkerSize',12);axis([0 200 0.00 9.00]);hold on;
% errorbar(N2,MAC_D_sim_54,err_MAC_D_54, 'o-','color', [0.8500, 0.3250, 0.0980],'LineWidth', 2,'MarkerSize',8);
% xlabel('Number of geophones');ylabel('Delay (ms)'); legend('Analytical 54Mbps','Simulation 54Mbps');
% %title('Delay: Analytical Vs Simulation');
% %ax = gca;
% %ax.YAxis.Exponent = -3;
% grid on; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(N1, P_f1(:,:,end), 's-','color', [0, 0.4470, 0.7410],'LineWidth', 1.5,'MarkerSize',6);axis([0 180 0.00 1.0]);
% title('Collision Probability');xlabel('Number of geophones');ylabel('Collision Probability'); %legend('Analytical','Simulation');
% grid on; 


%Collision prob.
% figure();
% plot(N1, P_f1(:,:,end), 'ks-');%axis([0 180 ]);