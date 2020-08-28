% *************************************************************************
%                       Anantharam & Schenk
%               Splanchnic Nerve Chromafin Cell Model
% *************************************************************************


%
%   This model simulates the electrophysiology of a single splanchnic nerve
%   - chromaffin cell synapse.  This model incorporates selected voltage gated Na, K and Ca
%   channels, buffering of Ca, and two populations of synaptic vesicles in the splanchnic buton.
%   This model also incorporates AChE and ACh spillover/diffusion in the synaptic space.
%   Lastly this model incorportates nAChR dynamics to measure expected EPSCs.
%   See accociated excel sheet for aditional details.

%   Model originally created on     22 April 2020
%   Model last modfied on           28 August 2020
%
%   Developed by        Arun Anantharam
%                       Noah Schenk
%

clear
clc

%% Frequency Selections
fs = [16.66]; %Set desired experimental stimulation frequency(s)

for freq = 1:length(fs)
    %% Experiment Type
    fs(freq)
    
    Calcium_channel_inwardTF = 0;
    %0= Calcium current comletly controled by channel condution and driving force; 1 = eliminates all outward Ca currents
    
    VGCC_TF = 1;
    %0 = No VGCC; 1 = VGCC present;
    
    Fast_bufferTF = 1;
    %0 = no fast buffer; 1 = fast buffer present;
    
    Slow_bufferTF = 1;
    %0 = no slow buffer; 1 = slow buffer present;
    
    Ca_ExtrusionTF = 1;
    %0 = no PM calcium extrusion; 1 = PM calcium pumps/exchangers present
    
    Ca_leakTF = 1;
    %0 = no Ca leak; 1 = Ca leak present
    
    newGran_TF = 1;
    %0 = no replenishment; 1 = granule replenishment
    
    syt7KO_TF = 0;
    %0 = WT; 1 = KO
    
    EGTA_TF = 0;
    %0 = no EGTA; 1 = EGTA present
    
    % **indicates Troubleshooting plot
    Plot_TF = [ 
        1 %1 ephys pannel
        1 %2 Intracellular Ca vs Time
        1 %3 Cumualitive Fusion vs Time
        1 %4 Ion channel steady state Probs**
        1 %5 Calcium Channel Probability and Vm vs Time
        1 %6 Calcium Extrusion rate vs [Ca]**
        1 %7 Calcium Driving force vs [Ca]**
        1 %8 Average [Ca] and Vm vs Time
        1 %9 Fusion Rates vs Time
        1 %10 Slow Buffer
        1 %11 Fast Buffer
        1 %12 EPSCs
        1 %13 EPSCs Extra**
        1 %14 ACh cobnc
        1 %15 nAChR Kinet
        ];
    
    %%
    
    %% Paramater Selection
    
    dt=0.001; % Time is in ms
    
    num_pulse = 2; 
    start_time = 17; %ms
    length_pulse = 0.1; %ms
    frequency_pulse = fs(freq); 
    pulse_strength = 600; %um/cm^2 
    Num_granules = 148;  
    Syt1percent = 0.53; 
    Syt7percent = 1 - Syt1percent;
    
    if( syt7KO_TF == 1)
        Num_granules = 81; 
        Syt7percent = 0;
        Syt1percent = 1;
    end
    
    Starting_Ca_SV = 50; %nM Intracellular [Ca] near plasma membrane
    Starting_Ca_Ct = 50; %nM Intracellular [Ca] in central cytoplasm
    Extracellular_Ca = 2*10^6; %nM 
    
    
    
    
    gap_pulse = 1/frequency_pulse*1000; %Length of gap between electrical stimulations
    T_end = start_time+(num_pulse)*gap_pulse+300; %Last time point of simulaton
    t= 0:dt:T_end;
    I_ext = zeros(1,length(t)); %External current applied
    for i = 1:num_pulse
        I_ext(t>(start_time+(i-1)*gap_pulse) & t< (start_time+(i-1)*gap_pulse+length_pulse)) = pulse_strength;
    end
    
    injectACh = zeros(1,length(t));
    
    
    
    v0 = -77; % mV; starting Vm
    
    C_m=1; % uF/cm^2; Membrane Capacitance is 
    
    E_K=-90; %mV
    g_K=150; %mS/cm^2
    
    E_Na=60; %mV
    g_Na=300; %mS/cm^2
    
    E_Ca_list = zeros(1, length(t))+(61.5/2*log( Extracellular_Ca)/log(Starting_Ca_SV));
    E_Ca_listCt = zeros(1, length(t))+(61.5/2*log( Extracellular_Ca)/log(Starting_Ca_Ct));
    if(VGCC_TF == 1)
        g_Ca = .032; % mS; Max Ca conduction
    else
        g_Ca = 0;
    end
    E_L=-77; %mV
    g_L=1; %mS/cm^2
    
    if( Ca_leakTF == 1)
        g_Ca_leakSV = 0.000001; %mS; max Ca leak conduction
        g_Ca_leakCt = 0.000001;
    else
        g_Ca_leakSV = 0;
        g_Ca_leakCt = 0;
    end
    
    
 
    
    AZ_SA = 9*10^-8; %cm^2 = 9 um^2 
    AZ_vol = AZ_SA*250*10^-6; %cm^3
    PreSyn_Vol = 2.45*10^-10; % cm^3
    PreSyn_SA = 1.81*10^-6; %cm^2
    Vol_rat = AZ_vol/PreSyn_Vol;
    
    PercentCaEqMap_slopes = diff(1 - exp(-0.2*(t)));
    PercentBFEqMap_slopes = diff(1-exp(-0.5/300*(t)));
    
    PercentEGTAMap_slopes = diff(1 - exp(-0.15*(t)));
    
    if (Ca_ExtrusionTF == 1)
        Ca_extrus_rate = 100*2*10^-6; %nmol/cm^2*ms
    else
        Ca_extrus_rate = 0;
    end
    
    Km_ca = 830; %nmol
    
    if (Slow_bufferTF == 1)
        startingConc_Sbuf = 250000; %nM; starting [Slow Ca buffer]
    else
        startingConc_Sbuf = 0;
    end
    if (Fast_bufferTF ==1)
        startingConc_Fbuf = 45000; %nM starting [Fast Ca buffer]
    else
        startingConc_Fbuf = 0;
    end
    K_dF = 1000;% nM  
    K_dS = 40;%nM   
    k_plusF = 1*10^-4; %/nM/mS  
    k_plusS = 5*10^-7; % /nM/mS 
    k_minF = K_dF*k_plusF;
    k_minS = K_dS*k_plusS;
    
    Syt7_b = 1/2.8; 
    Syt7_KonCa = 0.06; % per uM per ms
    Syt7_KoffCa = 0.6; % per ms
    Syt7_Krelease = 3; % per ms
    

    
    Syt1_b = 1/1.9; 
    Syt1_KonCa = 0.19;  % per uM per ms .15
    Syt1_KoffCa = 11.97; % per ms
    Syt1_Krelease = 12; %per ms 
    
    %ACh stuff
    quantalACh = 6000;
    avNum = 6.022*10^23;
    SC_Vol = AZ_SA*2.5*10^-6; %cm^3
    
    
    startingACh = 0.000*SC_Vol*10^-3*avNum; %numMoleculues
    startingACh_conc = startingACh/avNum/SC_Vol*10^12; %nM
    PercentAChMap_slopes = diff(1-exp(-0.5*(t)));
    
    

    AChE_kon1 = 2*10^-4; %/nM/ms
    AChE_koff1 = 1; %/ms; 
    AChE_kon2 = 110; %/ms; 
    AChE_koff2 = 1; %/ms; 
    AChE_k3 = 20; %/ms; 
    AChE_k4 = 10; %/ms; 
    startingAChE_conc = 750000; %nM
    

    nACHR_a = 1.2; % /ms 
    nACHR_B = 0.46; % /ms 
    nACHR_kd = 0.024; % /ms 
    nACHR_kon = 3*10^-3; %/ms/nmol  
    nACHR_Koff = 930; % /ms
    nACHR_Krec = 0.005; 
    nACHR_Kdead = 0.0000;
    nACHR_Krevive = 0.00000;
    
    num_nAChR = 2600;
    nAChR_con = 40; %pS
    inacAChR = zeros(1,length(t))+num_nAChR;
    acAChR = zeros(1,length(t));
    opAChR = zeros(1,length(t));
    desentAChR = zeros(1,length(t));
    deadAChR = zeros(1,length(t));
    
    
    ACh = zeros(1, length(t))+startingACh_conc;
    AChE = zeros(1, length(t))+startingAChE_conc;
    RateACh = zeros(1,length(t));
    MMc1 = zeros(1,length(t));
    MMc2 = zeros(1,length(t));
    AChE_A = zeros(1,length(t));
    
    diff2 = 0;
    
    %Ca Chanel Kinetics
    
    num_CaCh = 2007; 
    PQa1 = 5.89;
    PQa2 = 9.21;
    PQa3 = 5.20;
    PQa4 = 1823.18;
    PQa = 247.71;
    PQb1 = 14.99;
    PQb2 = 6.63;
    PQb3 = 132.80;
    PQb4 = 248.58;
    PQb = 8.28;
    PQk1 = 62.61;
    PQk2 = 33.92;
    PQk3 = 135.08;
    PQk4 = 20.86;
    
    EGTA_starting = 1000000; %nM
    EGTA_Kd = 200; %nM
    EGTA_Kp = 1.5*10^-6; %/nM/ms 
    EGTA_Km = EGTA_Kd*EGTA_Kp;
    
    %% Vector Setup
    v=zeros(1,length(t))+v0;
    n=zeros(1,length(t))+n_inf(v0);
    m=zeros(1,length(t))+m_inf(v0);
    h=zeros(1,length(t))+h_inf(v0);
    
    Gna = zeros(1,length(t));
    Gk = zeros(1,length(t));
    Gca = zeros(1,length(t));
    Ik = zeros(1,length(t));
    Ina = zeros(1,length(t));
    Ica = zeros(1,length(t));
    Il = zeros(1,length(t));
    
    Ca_freeIonSV = zeros(1,length(t))+ Starting_Ca_SV;
    Ca_freeIonCt = zeros(1,length(t))+ Starting_Ca_Ct;
    
    RateBfbound = zeros(1,length(t));
    RateBfnot = zeros(1,length(t));
    
    RateCaDiff = zeros(1,length(t));
    PercentCaEq = zeros(1,length(t));
    PercentCaEqMap = zeros(1,length(t));
    
    PercentCaEq(1) = 1-(abs(Ca_freeIonSV(1)-Ca_freeIonCt(1)))/(Ca_freeIonCt(1)+Ca_freeIonSV(1));
    PercentCaEqMap(1) = round(-5*log(1-PercentCaEq(1)+0.0000001)); %12.5
    RateCaDiff(1) = PercentCaEqMap_slopes(PercentCaEqMap(1)/dt+1);
    
    ExtrusionSV = zeros(1,length(t));
    ExtrusionCt = zeros(1,length(t));
    Ica_leak = zeros(1,length(t));
    Ica_leakSV = zeros(1,length(t));
    Ica_leakCt = zeros(1,length(t));
    
    Cyto_leak  = zeros(1,length(t));
    
    FastBufferChangeSV  = zeros(1,length(t));
    FastBufferChangeCt = zeros(1,length(t));
    
    SlowBufferChangeSV  = zeros(1,length(t));
    SlowBufferChangeCt  = zeros(1,length(t));
    
    FastBufferconcSV  = zeros(1,length(t))+startingConc_Fbuf - 1100;
    FastBufferconcCt  = zeros(1,length(t))+startingConc_Fbuf - 1700;
    
    SlowBufferconcSV  = zeros(1,length(t))+startingConc_Sbuf;
    SlowBufferconcCt = zeros(1,length(t))+startingConc_Sbuf;
    
    CaBound_fastSV  = zeros(1,length(t))+1100;
    CaBound_fastCt  = zeros(1,length(t))+1700;
    
    CaBound_slowSV  = zeros(1,length(t));
    CaBound_slowCt  = zeros(1,length(t));
    
    Syt7_conc = zeros(1,length(t)) + Num_granules.*Syt7percent; %Normalized
    Syt7Ca1 = zeros(1,length(t));
    Syt7Ca2 = zeros(1,length(t));
    Syt7Ca3 = zeros(1,length(t));
    Syt7Ca4 = zeros(1,length(t));
    Syt7Ca5 = zeros(1,length(t));
    Syt7Ca6 = zeros(1,length(t));
    Syt7released = zeros(1,length(t));
    
    Syt1_conc = zeros(1,length(t)) + Num_granules.*Syt1percent; %Normalized
    Syt1Ca1 = zeros(1,length(t));
    Syt1Ca2 = zeros(1,length(t));
    Syt1Ca3 = zeros(1,length(t));
    Syt1Ca4 = zeros(1,length(t));
    Syt1Ca5 = zeros(1,length(t));
    Syt1released = zeros(1,length(t));
    
    caCh0 = zeros(1,length(t)) +num_CaCh;
    caCh1 = zeros(1, length(t));
    caCh2 = zeros(1,length(t));
    caCh3 = zeros(1, length(t));
    caCh4 = zeros(1,length(t));
    caChO = zeros(1,length(t));
    
    nicI = zeros(1,length(t));
    
    if (EGTA_TF == 1)
        EGTA_freeCt = zeros(1,length(t)) + EGTA_starting;
        EGTA_freeSV = zeros(1,length(t)) + EGTA_starting;
    else
        EGTA_freeCt = zeros(1,length(t)) + 0;
        EGTA_freeSV = zeros(1,length(t)) + 0;
    end
    
    EGTA_BoundCt = zeros(1,length(t));
    EGTA_BoundSV = zeros(1,length(t));
    
    %% Simualation
    for i = 1:length(t)-1
        
        
        
        %Calcium Driving Force
        E_Ca_list(i) = 61.5/2*log( Extracellular_Ca/Ca_freeIonSV(i));
        E_Ca = E_Ca_list(i);
        E_Ca_listCt(i) = 61.5/2*log( Extracellular_Ca/Ca_freeIonCt(i));
        
        
        % VGCC State changes
        Ch0_1(i+1) = PQa1.*exp(v(i)./PQk1).*caCh0(i).*dt;
        if(Ch0_1(i+1) > caCh0(i))
            Ch0_1(i+1) = caCh0(i);
        end
        if(Ch0_1(i+1) < 0)
            Ch0_1(i+1) = 0;
        end
        
        
        Ch1_2(i+1) = PQa2.*exp(v(i)./PQk2).*caCh1(i).*dt;
        if(Ch1_2(i+1) > caCh1(i))
            Ch1_2(i+1) = caCh1(i);
        end
        if(Ch1_2(i+1) < 0)
            Ch1_2(i+1) = 0;
        end
        
        
        Ch2_3(i+1) = PQa3.*exp(v(i)./PQk3).*caCh2(i).*dt;
        if(Ch2_3(i+1) > caCh2(i))
            Ch2_3(i+1) = caCh2(i);
        end
        if(Ch2_3(i+1) < 0)
            Ch2_3(i+1) = 0;
        end
        
        
        Ch3_4(i+1) = PQa4.*exp(v(i)./PQk4).*caCh3(i).*dt;
        if(Ch3_4(i+1) > caCh3(i))
            Ch3_4(i+1) = caCh3(i);
        end
        if(Ch3_4(i+1) < 0)
            Ch3_4(i+1) = 0;
        end
        
        
        Ch4_O(i+1) = PQa.*caCh4(i).*dt;
        if(Ch4_O(i+1) > caCh4(i))
            Ch4_O(i+1) = caCh4(i);
        end
        if(Ch4_O(i+1) < 0)
            Ch4_O(i+1) = 0;
        end
        
        
        ChO_4(i+1) = PQb.*caChO(i).*dt;
        if(ChO_4(i+1) > caChO(i))
            ChO_4(i+1) = 0;
        end
        if(ChO_4(i+1) < 0)
            ChO_4(i+1) = 0;
        end
        
        
        Ch4_3(i+1) = PQb4.*exp(-1.*v(i)./PQk4).*caCh4(i).*dt;
        if(Ch4_3(i+1) > caCh4(i))
            Ch4_3(i+1) = caCh4(i);
        end
        if(Ch4_3(i+1) < 0)
            Ch4_3(i+1) = 0;
        end
        
        
        Ch3_2(i+1) = PQb3.*exp(-1.*v(i)./PQk3).*caCh3(i).*dt;
        if(Ch3_2(i+1) > caCh3(i))
            Ch3_2(i+1) = caCh3(i);
        end
        if(Ch3_2(i+1) < 0)
            Ch3_2(i+1) = 0;
        end
        
        
        Ch2_1(i+1) = PQb2.*exp(-1.*v(i)./PQk2).*caCh2(i).*dt;
        if(Ch2_1(i+1) > caCh2(i))
            Ch2_1(i+1) = caCh2(i);
        end
        if(Ch2_1(i+1) < 0)
            Ch2_1(i+1) = 0;
        end
        
        
        Ch1_0(i+1) = PQb1.*exp(-1.*v(i)./PQk1).*caCh1(i).*dt;
        if(Ch1_0(i+1) > caCh1(i))
            Ch1_0(i+1) = caCh1(i);
        end
        if(Ch1_0(i+1) < 0)
            Ch1_0(i+1) = 0;
        end
        
        caCh0(i+1) = caCh0(i) + (Ch1_0(i) - Ch0_1(i));
        
        caCh1(i+1) = caCh1(i) + (Ch2_1(i) - Ch1_2(i) + Ch0_1(i) - Ch1_0(i));
        
        caCh2(i+1) = caCh2(i) + (Ch3_2(i) - Ch2_3(i) + Ch1_2(i) - Ch2_1(i));
        
        caCh3(i+1) = caCh3(i) + (Ch4_3(i) - Ch3_4(i) + Ch2_3(i) - Ch3_2(i));
        
        caCh4(i+1) = caCh4(i) + (ChO_4(i) - Ch4_O(i) + Ch3_4(i) - Ch4_3(i));
        
        caChO(i+1) = caChO(i) + (Ch4_O(i) - ChO_4(i));
        
        
        %Ion channel conduction
        Gna(i) = g_Na*m(i)^3*(h(i));
        Gk(i) = g_K*n(i);
        Gca(i) = g_Ca.*caChO(i)./num_CaCh;
        
        %Ion channel Current
        Ina(i) = Gna(i)*(v(i)-E_Na);
        Ik(i) = Gk(i)*(v(i)-E_K);
        Ica(i) = Gca(i)*(v(i)-E_Ca_list(i));
        Il(i) = g_L*(v(i)-E_L);
        Ica_leakSV(i) = g_Ca_leakSV * (v(i)-E_Ca_list(i));
        Ica_leakCt(i) = g_Ca_leakCt * (v(i)-E_Ca_listCt(i));
        
        %If == 1 Makes Ca inward only
        if (Calcium_channel_inwardTF == 1)
            if(Ica(i) > 0)
                Ica(i) = 0;
            end
        end
        
        %Channel Probabilities for next time point
        n(i+1)=n(i)+dt*(n_inf(v(i))-n(i))./(n_tau(v(i)));
        m(i+1)=m(i)+dt*(m_inf(v(i))-m(i))./(m_tau(v(i)));
        h(i+1)=h(i)+dt*(h_inf(v(i))-h(i))./(h_tau(v(i)));
        
        %Free Ca diffusion
        PercentCaEq(i) = 1-(abs(Ca_freeIonSV(i)-Ca_freeIonCt(i)))/(Ca_freeIonCt(i)+Ca_freeIonSV(i));
        PercentCaEqMap(i) = round(-5*log(1-PercentCaEq(i)+0.000000001));
        RateCaDiff(i) = PercentCaEqMap_slopes(PercentCaEqMap(i)/dt+1);
        
        %Fast Buffer diffusion
        RateBfbound(i) = buffer_D(CaBound_fastSV(i), CaBound_fastCt(i), PercentBFEqMap_slopes, dt);
        RateBfnot(i) = buffer_D(FastBufferconcSV(i), FastBufferconcCt(i), PercentBFEqMap_slopes, dt);
        
        %EGTA diffusion
        RateEGTAbound(i) = EGTA_D(EGTA_BoundSV(i), EGTA_BoundCt(i), PercentEGTAMap_slopes, dt);
        RateEGTAnot(i) = EGTA_D(EGTA_freeSV(i), EGTA_freeCt(i), PercentEGTAMap_slopes, dt);
        
        %Calcium conc from channels, extrusion, and leak
        Ca_freeIonSV(i) = Ca_freeIonSV(i)-((Ica(i)+Ica_leakSV(i))*AZ_SA*1000)*dt/(2*9.648*10^13)/(AZ_vol)*1000*1000*1000*1000; % in nM
        ExtrusionSV(i) = ( Ca_extrus_rate*AZ_SA*dt*Ca_freeIonSV(i))/(AZ_vol*(Ca_freeIonSV(i)+Km_ca))/50;%/10 = tune
        Ca_freeIonSV(i) = Ca_freeIonSV(i) - ExtrusionSV(i);
        ExtrusionCt(i) = ( Ca_extrus_rate* PreSyn_SA *dt*Ca_freeIonCt(i))/( PreSyn_Vol *(Ca_freeIonCt(i)+Km_ca));
        Cyto_leak(i) = ((Ica_leakCt(i))*PreSyn_SA*1000)*dt/(2*9.648*10^13)/(PreSyn_Vol)*1000*1000*1000*1000;
        
        
        %Kinetics of Buffering
        FastBufferChangeSV(i+1) = (k_minF*CaBound_fastSV(i)-k_plusF*Ca_freeIonSV(i)*FastBufferconcSV(i))*dt;
        FastBufferChangeCt(i+1) = (k_minF*CaBound_fastCt(i)-k_plusF*Ca_freeIonCt(i)*FastBufferconcCt(i))*dt;
        
        SlowBufferChangeSV(i+1) = (k_minS*CaBound_slowSV(i)-k_plusS*Ca_freeIonSV(i)*SlowBufferconcSV(i))*dt;
        SlowBufferChangeCt(i+1) = (k_minS*CaBound_slowCt(i)-k_plusS*Ca_freeIonCt(i)*SlowBufferconcCt(i))*dt;
        
        EGTAbufferChangeSV(i+1) = (EGTA_Km*EGTA_BoundSV(i) - EGTA_Kp*Ca_freeIonSV(i)*EGTA_freeSV(i))*dt;
        EGTAbufferChangeCt(i+1) = (EGTA_Km*EGTA_BoundCt(i) - EGTA_Kp*Ca_freeIonCt(i)*EGTA_freeCt(i))*dt;
        
        %Syt Kinetics
        Syt7_conc(i+1) = Syt7_conc(i) - 6*Syt7_KonCa.*Syt7_conc(i).*Ca_freeIonSV(i)./1000.*dt   + Syt7_b.^0.*Syt7_KoffCa.*Syt7Ca1(i).*dt;
        Syt7Ca1(i+1)   = Syt7Ca1(i)   + 6*Syt7_KonCa.*Syt7_conc(i).*Ca_freeIonSV(i)./1000.*dt   + Syt7_b.^1.*Syt7_KoffCa.*Syt7Ca2(i).*dt - 5*Syt7_KonCa.*Syt7Ca1(i).*Ca_freeIonSV(i)./1000.*dt - Syt7_b.^0.*Syt7_KoffCa.*Syt7Ca1(i).*dt;
        Syt7Ca2(i+1)   = Syt7Ca2(i)   + 5*Syt7_KonCa.*Syt7Ca1(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt7_b.^2.*Syt7_KoffCa.*Syt7Ca3(i).*dt - 4*Syt7_KonCa.*Syt7Ca2(i).*Ca_freeIonSV(i)./1000.*dt - Syt7_b.^1.*Syt7_KoffCa.*Syt7Ca2(i).*dt;
        Syt7Ca3(i+1)   = Syt7Ca3(i)   + 4*Syt7_KonCa.*Syt7Ca2(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt7_b.^3.*Syt7_KoffCa.*Syt7Ca4(i).*dt - 3*Syt7_KonCa.*Syt7Ca3(i).*Ca_freeIonSV(i)./1000.*dt - Syt7_b.^2.*Syt7_KoffCa.*Syt7Ca3(i).*dt;
        Syt7Ca4(i+1)   = Syt7Ca4(i)   + 3*Syt7_KonCa.*Syt7Ca3(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt7_b.^4.*Syt7_KoffCa.*Syt7Ca5(i).*dt - 2*Syt7_KonCa.*Syt7Ca4(i).*Ca_freeIonSV(i)./1000.*dt - Syt7_b.^3.*Syt7_KoffCa.*Syt7Ca4(i).*dt;
        Syt7Ca5(i+1)   = Syt7Ca5(i)   + 2*Syt7_KonCa.*Syt7Ca4(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt7_b.^5.*Syt7_KoffCa.*Syt7Ca6(i).*dt - 1*Syt7_KonCa.*Syt7Ca5(i).*Ca_freeIonSV(i)./1000.*dt - Syt7_b.^4.*Syt7_KoffCa.*Syt7Ca5(i).*dt;
        Syt7Ca6(i+1)   = Syt7Ca6(i)   + 1*Syt7_KonCa.*Syt7Ca5(i)  .*Ca_freeIonSV(i)./1000.*dt   + 0                                      - Syt7_Krelease.*Syt7Ca6(i).*dt                       - Syt7_b.^5.*Syt7_KoffCa.*Syt7Ca6(i).*dt;
        Syt7released(i+1) = Syt7released(i)+Syt7_Krelease.*Syt7Ca6(i).*dt;
        
        
        
        Syt1_conc(i+1) = Syt1_conc(i) - 5*Syt1_KonCa.*Syt1_conc(i).*Ca_freeIonSV(i)./1000.*dt   + Syt1_b.^0.*Syt1_KoffCa.*Syt1Ca1(i).*dt;
        Syt1Ca1(i+1)   = Syt1Ca1(i)   + 5*Syt1_KonCa.*Syt1_conc(i).*Ca_freeIonSV(i)./1000.*dt   + Syt1_b.^1.*Syt1_KoffCa.*Syt1Ca2(i).*dt   - 4*Syt1_KonCa.*Syt1Ca1(i).*Ca_freeIonSV(i)./1000.*dt - Syt1_b.^0.*Syt1_KoffCa.*Syt1Ca1(i).*dt;
        Syt1Ca2(i+1)   = Syt1Ca2(i)   + 4*Syt1_KonCa.*Syt1Ca1(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt1_b.^2.*Syt1_KoffCa.*Syt1Ca3(i).*dt   - 3*Syt1_KonCa.*Syt1Ca2(i).*Ca_freeIonSV(i)./1000.*dt - Syt1_b.^1.*Syt1_KoffCa.*Syt1Ca2(i).*dt;
        Syt1Ca3(i+1)   = Syt1Ca3(i)   + 3*Syt1_KonCa.*Syt1Ca2(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt1_b.^3.*Syt1_KoffCa.*Syt1Ca4(i).*dt   - 2*Syt1_KonCa.*Syt1Ca3(i).*Ca_freeIonSV(i)./1000.*dt - Syt1_b.^2.*Syt1_KoffCa.*Syt1Ca3(i).*dt;
        Syt1Ca4(i+1)   = Syt1Ca4(i)   + 2*Syt1_KonCa.*Syt1Ca3(i)  .*Ca_freeIonSV(i)./1000.*dt   + Syt1_b.^4.*Syt1_KoffCa.*Syt1Ca5(i).*dt   - 1*Syt1_KonCa.*Syt1Ca4(i).*Ca_freeIonSV(i)./1000.*dt - Syt1_b.^3.*Syt1_KoffCa.*Syt1Ca4(i).*dt;
        Syt1Ca5(i+1)   = Syt1Ca5(i)   + 1*Syt1_KonCa.*Syt1Ca4(i)  .*Ca_freeIonSV(i)./1000.*dt   + 0                                          - Syt1_Krelease.*Syt1Ca5(i).*dt                       - Syt1_b.^4.*Syt1_KoffCa.*Syt1Ca5(i).*dt;
        Syt1released(i+1) = Syt1released(i)+Syt1_Krelease.*Syt1Ca5(i).*dt;
        
        if ((newGran_TF ==1) && (t(i) > 17)) 
            if((syt7KO_TF == 1) && (Syt1_conc(i+1) <Num_granules*0.42))
                Syt1_conc(i+1) = Syt1_conc(i+1) + 2.9*dt; 
            end
            if((syt7KO_TF == 1) && (Syt1_conc(i+1) <Num_granules*0.78))
                Syt1_conc(i+1) = Syt1_conc(i+1) + 0.090*dt; 
            end
            if((syt7KO_TF == 1) && (Syt1_conc(i+1) <Num_granules*1))
                Syt1_conc(i+1) = Syt1_conc(i+1) + 0.003*dt; 
            end
            
            
            
            if(syt7KO_TF == 0 && (Syt1_conc(i+1) < (Num_granules*Syt1percent*0.54)))
                Syt1_conc(i+1) = Syt1_conc(i+1) + 4.9.*dt; 
            end
            if(syt7KO_TF == 0 && (Syt1_conc(i+1) < (Num_granules*Syt1percent*0.89)))
                Syt1_conc(i+1) = Syt1_conc(i+1) + 0.13.*dt; 
            end
            if(syt7KO_TF == 0 && (Syt1_conc(i+1) < (Num_granules*Syt1percent*1)))
                Syt1_conc(i+1) = Syt1_conc(i+1) + 0.00005.*dt; 
            end
            
            if(syt7KO_TF == 0 && (Syt7_conc(i+1) < (Num_granules*Syt7percent*0.10)))
                Syt7_conc(i+1) = Syt7_conc(i+1) + 1.1.*dt;  
            end
            if(syt7KO_TF == 0 && (Syt7_conc(i+1) < (Num_granules*Syt7percent*0.56)))
                Syt7_conc(i+1) = Syt7_conc(i+1) + 0.27.*dt;  
            end
            if(syt7KO_TF == 0 && (Syt7_conc(i+1) < (Num_granules*Syt7percent*1)))
                Syt7_conc(i+1) = Syt7_conc(i+1) + 0.001.*dt;
            end
            
        end
        
        
        
        
        % Final Buffer conc determiniation
        if(CaBound_fastSV(i) > CaBound_fastCt(i))
            CaBound_fastSV(i) = CaBound_fastSV(i) - RateBfbound(i)*CaBound_fastSV(i);
            CaBound_fastCt(i) = CaBound_fastCt(i) + Vol_rat*RateBfbound(i)*CaBound_fastSV(i);
        else
            CaBound_fastSV(i) = CaBound_fastSV(i) + RateBfbound(i)*CaBound_fastCt(i);
            CaBound_fastCt(i) = CaBound_fastCt(i) - Vol_rat*RateBfbound(i)*CaBound_fastCt(i);
        end
        
        
        if(FastBufferconcSV(i) > FastBufferconcCt(i))
            FastBufferconcSV(i) = FastBufferconcSV(i) - RateBfnot(i)*FastBufferconcSV(i);
            FastBufferconcCt(i) = FastBufferconcCt(i) + Vol_rat*RateBfnot(i)*FastBufferconcSV(i);
        else
            FastBufferconcSV(i) = FastBufferconcSV(i) + RateBfnot(i)*FastBufferconcCt(i);
            FastBufferconcCt(i) = FastBufferconcCt(i) - Vol_rat*RateBfnot(i)*FastBufferconcCt(i);
        end
        
        if(EGTA_BoundSV(i) > EGTA_BoundCt(i))
            EGTA_BoundSV(i) = EGTA_BoundSV(i) - RateEGTAbound(i)*EGTA_BoundSV(i);
            EGTA_BoundCt(i) = EGTA_BoundCt(i) + Vol_rat*RateEGTAbound(i)*EGTA_BoundSV(i);
        else
            EGTA_BoundSV(i) = EGTA_BoundSV(i) + RateEGTAbound(i)*EGTA_BoundCt(i);
            EGTA_BoundCt(i) = EGTA_BoundCt(i) - Vol_rat*RateEGTAbound(i)*EGTA_BoundCt(i);
        end
        
        if(EGTA_freeSV(i) > EGTA_freeCt(i))
            EGTA_freeSV(i) = EGTA_freeSV(i) - RateEGTAnot(i)*EGTA_freeSV(i);
            EGTA_freeCt(i) = EGTA_freeCt(i) + Vol_rat*RateEGTAnot(i)*EGTA_freeSV(i);
        else
            EGTA_freeSV(i) = EGTA_freeSV(i) + RateEGTAnot(i)*EGTA_freeCt(i);
            EGTA_freeCt(i) = EGTA_freeCt(i) - Vol_rat*RateEGTAnot(i)*EGTA_freeCt(i);
        end
        
        FastBufferconcSV(i+1) = FastBufferconcSV(i)+FastBufferChangeSV(i+1);
        FastBufferconcCt(i+1) = FastBufferconcCt(i)+FastBufferChangeCt(i+1);
        
        SlowBufferconcSV(i+1) = SlowBufferconcSV(i)+SlowBufferChangeSV(i+1);
        SlowBufferconcCt(i+1) = SlowBufferconcCt(i)+SlowBufferChangeCt(i+1);
        
        CaBound_fastSV(i+1) = CaBound_fastSV(i)-FastBufferChangeSV(i+1);
        CaBound_fastCt(i+1) = CaBound_fastCt(i)-FastBufferChangeCt(i+1);
        
        CaBound_slowSV(i+1) = CaBound_slowSV(i)-SlowBufferChangeSV(i+1);
        CaBound_slowCt(i+1) = CaBound_slowCt(i)-SlowBufferChangeCt(i+1);
        
        EGTA_BoundSV(i+1) = EGTA_BoundSV(i) - EGTAbufferChangeSV(i+1);
        EGTA_BoundCt(i+1) = EGTA_BoundCt(i) - EGTAbufferChangeCt(i+1);
        
        EGTA_freeSV(i+1) = EGTA_freeSV(i) + EGTAbufferChangeSV(i+1);
        EGTA_freeCt(i+1) = EGTA_freeCt(i) + EGTAbufferChangeCt(i+1);
        
        
        %Final Ca Conc determinination
        if(Ca_freeIonSV(i) > Ca_freeIonCt(i))
            Ca_freeIonCt(i+1) = Ca_freeIonCt(i)+ Vol_rat*RateCaDiff(i)*Ca_freeIonSV(i) - ExtrusionCt(i) - Cyto_leak(i)+FastBufferChangeCt(i)+SlowBufferChangeCt(i)+EGTAbufferChangeCt(i);
            Ca_freeIonSV(i+1) = Ca_freeIonSV(i)-RateCaDiff(i)*Ca_freeIonSV(i)+FastBufferChangeSV(i)+SlowBufferChangeSV(i)+EGTAbufferChangeSV(i);
        else
            Ca_freeIonCt(i+1) = Ca_freeIonCt(i)- Vol_rat*RateCaDiff(i)*Ca_freeIonCt(i) - ExtrusionCt(i) - Cyto_leak(i)+FastBufferChangeCt(i)+SlowBufferChangeCt(i)+EGTAbufferChangeCt(i);
            Ca_freeIonSV(i+1) = Ca_freeIonSV(i) + RateCaDiff(i)*Ca_freeIonCt(i)+FastBufferChangeSV(i)+SlowBufferChangeSV(i)+EGTAbufferChangeSV(i);
        end
        
        %Determine voltage for next time point
        v(i+1)=v(i)+dt*(-Ina(i)-Ik(i)-Ica(i)-Il(i)+I_ext(i))/C_m;
        
        
        %Ach post-synaptic sim
        
        RateACh(i) = ACh_D(PercentAChMap_slopes, dt);
        
        AChon(i)   = AChE(i).*ACh(i).*AChE_kon1.*dt;
        if(AChon(i) < 0)
            AChon(i) = 0;
        end
        
        AChoff(i)  = AChE_koff1.*MMc1(i).*dt;
        if(AChoff(i) < 0)
            AChoff(i) = 0;
        end
        
        AChDiff(i) = RateACh(i)*ACh(i);
        if(AChDiff(i) < 0)
            AChDiff(i) = 0;
        end
        
        AChMetab_half(i) = MMc1(i).*AChE_kon2.*dt;
        if(AChMetab_half(i) < 0)
            AChMetab_half(i) = 0;
        end
        
        AChMetab_full(i) = MMc2(i).*AChE_k3.*dt;
        if(AChMetab_full(i) < 0)
            AChMetab_full(i) = 0;
        end
        
        AChMetab_off(i) = MMc2(i).*AChE_koff2.*dt;
        if(AChMetab_off(i) < 0)
            AChMetab_off(i) = 0;
        end
        
        newAChE(i) = AChE_A(i).*AChE_k4.*dt;
        if(newAChE(i) < 0)
            newAChE(i) = 0;
        end
        
        dEnzyme(i) = 0 - AChon(i) + AChoff(i) + newAChE(i);
        if(-1*dEnzyme(i) > AChE(i))
            diff2(i) = AChE(i) + dEnzyme(i);
            AChon(i) = AChon(i) + diff2(i);
            dEnzyme(i) = 0 - AChon(i) + AChoff(i) + newAChE(i);
        end
        dc1(i) = 0 + AChon(i) - AChoff(i) - AChMetab_half(i) + AChMetab_off(i);
        dc2(i) = 0 + AChMetab_half(i) - AChMetab_off(i) - AChMetab_full(i);
        de_A(i) = 0 + AChMetab_full(i) - newAChE(i);
        
        dACh(i) = 0 - AChon(i) - AChDiff(i) + AChoff(i);
        
        checkE(i) = dEnzyme(i) + dc1(i) + dc2(i) + de_A(i);
        
        
        if (t(i) > 50)
            startAvgTime(i) = i-50/dt;
        else
            startAvgTime(i) = 1;
        end
        
        injectACh(i+1) = 0; %Unused function, used to apply exogenous ACh on synapse
        
        releasedACh(i) = (quantalACh.*1.667.*(Syt7released(i+1)-Syt7released(i))+quantalACh.*(Syt1released(i+1)-Syt1released(i)))./avNum./SC_Vol*10^12;
        ACh(i+1) = ACh(i) + dACh(i) + injectACh(i) + releasedACh(i); %nM
        if(ACh(i+1) < 0)
            ACh(i+1) = 0;
        end
        
        MMc1(i+1) = MMc1(i) + dc1(i);
        MMc2(i+1) = MMc2(i) + dc2(i);
        AChE_A(i+1) = AChE_A(i) + de_A(i);
        AChE(i+1) = startingAChE_conc - MMc1(i+1) - MMc2(i+1) - AChE_A(i+1);
        if(AChE(i+1) < 0)
            AChE(i+1) = 0;
        end
        
        inacAChR(i+1)   = inacAChR(i)   - inacAChR(i).*ACh(i).*nACHR_kon.*dt + acAChR(i).*nACHR_Koff.*dt     + nACHR_Krec.*desentAChR(i).*dt     + nACHR_Krevive.*deadAChR(i).*dt;
        acAChR(i+1)     = acAChR(i)     + inacAChR(i).*ACh(i).*nACHR_kon.*dt - acAChR(i).*nACHR_Koff.*dt     - nACHR_B.*acAChR(i).*dt            + nACHR_a.*opAChR(i).*dt ;
        opAChR(i+1)     = opAChR(i)     + nACHR_B.*acAChR(i).*dt             - nACHR_a.*opAChR(i).*dt        - nACHR_kd.*opAChR(i).*dt;
        desentAChR(i+1) = desentAChR(i) + nACHR_kd.*opAChR(i).*dt            - nACHR_Krec.*desentAChR(i).*dt - nACHR_Kdead.*desentAChR(i).*dt;
        deadAChR(i+1)   = deadAChR(i)   + nACHR_Kdead.*desentAChR(i).*dt     - nACHR_Krevive.*deadAChR(i).*dt;
        
        %Nicotinic current calculation
        nicI(i+1) = opAChR(i+1).*nAChR_con.*E_Na/1000; %pA
        
    end
    
    %% Resample model output to data collection rates
    currentResamp = resample(nicI, 1,50)';
    calcium_trans = resample(Ca_freeIonSV, 1, 50)';
    currentResamp = currentResamp(1:1571).*-1;
    time = 1:0.05:78.5;
    Vm_Redone = resample(v,1,50);
    calcium_Ct = resample(Ca_freeIonCt,1,500);
    
    %% Vizualization
    
    if (Plot_TF(1) == 1)
        figure(1);
        sgtitle('Title');
        
        subplot(3,2,1);
        hold on
        plot(t,I_ext);
        title('I_s_t_i_m');
        xlabel('time (ms)');
        ylabel('Current (uA/cm^2)');
        hold off
        
        subplot(3,2,[3,4])
        hold on
        plot(t,v);
        title('V_m');
        xlabel("time (ms)");
        ylabel("Vm (mV)");
        axis([0 T_end -100 100]);
        hold off
        
        hold on
        subplot(3,2,2);
        % plot(t,m.^3.*h,t,n.^4,t,f.^2.*d.^2);
        plot(t,m.^3.*h,t,n.^4,t,caChO./num_CaCh);
        title("probabilities");
        legend("na", "k", "Ca");
        xlabel("time (ms)");
        ylabel('Probability')
        hold off;
        
        subplot(3,2,5);
        hold on
        plot(t,Gna,t,Gk, t,Gca);
        title('conduction')
        legend("gna","gk", "gca")
        xlabel("time (ms)");
        ylabel("Conductance (mS/cm^2)")
        hold off;
        
        subplot(3,2,6);
        hold on
        plot(t,Ik*AZ_SA*1000,t,Ina*AZ_SA*1000,t, Il*AZ_SA*1000, t, Ica*AZ_SA*1000);
        title("Currents");
        legend("Ik","INa", "IL" , "ICa");
        xlabel("time (ms)");
        ylabel("Current (nA)");
        hold off;
    end
    
    Vm = -150:1:150;
    ns = zeros(1,length(Vm));
    ms = zeros(1,length(Vm));
    hs = zeros(1,length(Vm));
    
    for i = 1:length(Vm)-1
        ns(i) = n_inf(Vm(i));
        ms(i) = m_inf(Vm(i));
        hs(i) = h_inf(Vm(i));
        
    end
    
    if (Plot_TF(2) == 1)
        figure(2)
        hold on
        plot(t, Ca_freeIonSV/1000, t, Ca_freeIonCt/1000)
        ylabel("[Ca] (um)");
        xlabel("time (ms)");
        legend("Free PM Ca", "Free Cytosol Ca", "location", "southeast");
        hold off
    end
    
    if (Plot_TF(3) == 1)
        figure(3)
        hold on
        plot( t, Syt1released+Syt7released, t, Syt1released, t, Syt7released);
        ylabel("Cumulitive Fusion");
        xlabel("time (ms)");
        legend("Cumulitive Fusion","Syt-1", "Syt-7", "location", "southeast");
        hold off
    end
    
    if (Plot_TF(4) == 1)
        figure(4)
        hold on
        plot(Vm,ns, Vm, ms, Vm, hs);
        legend("n_i_n_f","m_i_n_f","h_i_n_f");
        xlabel("Vm(mV)");
        hold off
    end
    
    if (Plot_TF(5) == 1)
        figure(5)
        hold on
        plot(t,v,t, caChO./num_CaCh  );
        xlabel("time (ms)");
        ylabel("Vm -- Ca opening probability");
        legend("Vm", "Ca-prob");
        hold off;
    end
    
    if (Plot_TF(6) == 1)
        figure(6);
        hold on;
        Ca_test = 0:10:300000;
        Ca_rate_test = (Ca_extrus_rate.*AZ_SA.*dt.*Ca_test)./(AZ_vol.*(Ca_test+Km_ca));
        plot(Ca_test,Ca_rate_test);
        xlabel("[Ca] nM");
        ylabel("Extrusion nM/dt");
        hold off
    end
    
    if (Plot_TF(7) == 1)
        figure(7);
        hold on;
        Ca_test = 0:10000:10000000;
        Extracellular_Ca = 1*10^6;
        E_Ca_test = 61.5./2.*log( Extracellular_Ca./Ca_test);
        plot(Ca_test, E_Ca_test);
        xlabel("[Ca] nM");
        ylabel("Ca Driving Force (mV)");
        hold off;
    end
    
    if (Plot_TF(8) == 1)
        figure(8);
        hold on;
        plot(t,v*100+15000,t, Ca_freeIonCt);
        xlabel("time (ms)");
        ylabel("Vm -- Free Ca");
        legend("Vm", "Average Cytosolic Calcium", "location", "northwest");
    end
    
    if (Plot_TF(9) == 1)
        figure(9)
        hold on
        plot( t(1:end-1), diff(Syt1released+Syt7released), t(1:end-1), diff(Syt1released), t(1:end-1), diff(Syt7released));
        ylabel("Fusion events per dt");
        xlabel("time (ms)");
        legend("Cumulitive","Syt-1", "Syt-7", "location", "southeast");
        hold off
    end
    
    
    if (Plot_TF(10) == 1)
        figure(10)
        hold on
        plot( t, SlowBufferconcCt, t, SlowBufferconcSV);
        ylabel("Slow Buffer Conc [nM]");
        xlabel("time (ms)");
        legend("Average Cytosol", "Close to PM");
        hold off
    end
    
    if (Plot_TF(11) == 1)
        figure(11)
        hold on
        plot( t, FastBufferconcCt, t, FastBufferconcSV);
        ylabel("Fast Buffer Conc [nM]");
        xlabel("time (ms)");
        legend("Average Cytosol", "Close to PM");
        hold off
    end
    
    if (Plot_TF(12) == 1)
        figure(12)
        hold on
        plot(t, -nicI,t,-ACh/100)
        ylabel("Current (pA)");
        xlabel('Time (ms)');
        legend("I_n_A_C_h_R", "[ACh]");
        
        axis([15 50 -250 0]);
        hold off
    end
    
    if (Plot_TF(13) == 1)
        
        [pks, locs, w, p] = findpeaks(nicI);
        PPR = p(p>2);
        PPR = PPR(2)/PPR(1);
        figure(13)
        hold on
        findpeaks(nicI, t, 'MinPeakDistance', 1, 'Annotate', 'extents');
        ylabel("Current (pA)");
        xlabel('Time (ms)');
        title(strcat("PPR = ", num2str(PPR)));
        hold off
        PPR
        PPR_f(freq) = PPR;
    end
    
    if (Plot_TF(14) == 1)
        figure(14)
        hold on
        plot(t, ACh)
        ylabel("ACh (nM)");
        xlabel('Time (ms)');
        hold off
    end
    
    if (Plot_TF(15) ==1)
        figure(15);
        hold on
        plot( t, acAChR, t, opAChR, t, desentAChR, t, ACh/200);
        legend("ac", "op", "desent", "ACh");
        ylabel("amount")
        hold off;
    end
    
    
    
end

%% Helpers Ion Channels

function alpha_n0 = alpha_n(v)
a = 0.02;
b = 0.002;
V_12 = 25-28;
w = v - V_12;
k = 9;
alpha_n0 = a.*w./(1-exp(-1.*w./k));
end
function beta_n0 = beta_n(v)
a = 0.02;
b = 0.002;
V_12 = 25-28;
w = v - V_12;
k = 9;
beta_n0 = -1.*b.*w./(1-exp(w./k));
end
function n_inf0 = n_inf(v)
n_inf0 = alpha_n(v)/(alpha_n(v)+beta_n(v));
end
function n_tau0 = n_tau(v)
n_tau0 = 1./(alpha_n(v)+beta_n(v))./3;
end

function alpha_m0 = alpha_m(v)
a = 0.182;
b = 0.124;
V_12 = -41;
w = v - V_12;
k = 6;
alpha_m0 = a.*w./(1-exp(-1.*w./k));
end
function beta_m0 = beta_m(v)
a = 0.182;
b = 0.124;
V_12 = -41;
w = v - V_12;
k = 6;
beta_m0 = -1.*b.*w./(1-exp(w./k));
end
function m_inf0 = m_inf(v)
m_inf0 = alpha_m(v)/(alpha_m(v)+beta_m(v));
end
function m_tau0 = m_tau(v)
m_tau0 = 1./(alpha_m(v)+beta_m(v))./3;
end

function alpha_h0 = alpha_h(v)
a = 0.024;
b = 0.0091;
V_12 = -48;
w = v - V_12;
k = 6.2;
alpha_h0 = a.*w./(1-exp(-1.*w./k));
end
function beta_h0 = beta_h(v)
a = 0.024;
b = 0.0091;
V_12 = -73;
w = v - V_12;
k = 6.2;
beta_h0 = -1.*b.*w./(1-exp(w./k));
end
function h_inf0 = h_inf(v)
h_inf0 = 1./(1+exp((v+70)./6.2));
end

function h_tau0 = h_tau(v)
a = 0.024;
b = 0.0091;
V_12a = -48;
V_12b = -73;
wa = v - V_12a;
wb = v-V_12b;
k = 5;
beta_h0 = -1.*b.*wb./(1-exp(wb./k));
alpha_h0 = a.*wa./(1-exp(-1.*wa./k));
h_tau0 = 1./(alpha_h0+beta_h0);
end


%% Helper(s) Diffusion
function buffer_D0 = buffer_D(SVconc, Ctconc, Map_slopes,dt)
percent_diff = 1-(abs(SVconc-Ctconc))/(SVconc+Ctconc+1);
percent_map = round(-600*log(1-percent_diff+0.000000001));
spot = percent_map/dt;
if(spot > length(Map_slopes)-2)
    spot = length(Map_slopes)-2;
end
buffer_D0 = Map_slopes(spot+1);
end

function ACh_D0 = ACh_D(Map_slopes,dt)
percent_diff = 0;
percent_map = round(-2*log(1-percent_diff+0.000000001));
spot = percent_map/dt;
if(spot > length(Map_slopes)-2)
    spot = length(Map_slopes)-2;
end
ACh_D0 = Map_slopes(spot+1);
end

function EGTA_D0 = EGTA_D(SVconc, Ctconc, Map_slopes,dt)
percent_diff = 1-(abs(SVconc-Ctconc))/(SVconc+Ctconc+1);
percent_map = round(-600*log(1-percent_diff+0.000000001));
spot = percent_map/dt;
if(spot > length(Map_slopes)-2)
    spot = length(Map_slopes)-2;
end
EGTA_D0 = Map_slopes(spot+1);
end
