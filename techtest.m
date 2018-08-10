%% Microwave Radio Links

% PRx [dB] -> free space propagation
% PTx [dB] -> power transmitted by the emitter
% GRx [dB] -> emitter gain
% GTx [dB] -> receiver gain
% Asis [dB] -> all losses related to equipment (OF, coaxial cables,
%receivers,...) Will be considered 3dB throughout the project.
% Ao [dB] -> Free Space attenuation
% d [km] -> Connection length
% f [GHz] -> Frequency (an array will be considered)
% lambda [m] -> Wavelength
% n_a -> Antenna efficiency
% d_a -> Antenna diameter
% Aabs [dB] -> Additional atmospheric attenuation
% gamma0_o [dB/km] -> attenuation coefficient due to oxygen
% gamma0_w [dB/km] -> attenuation coefficient due to water vapour
% gamma_r [dB/km] -> attenuation due to rain
% Ri [mm/h] -> Rain intensity
% N0 [dBW] -> Thermal noise power
% b_rf [Hz] -> Effective noise bandwidth
% b0 [Hz] -> Bandwith
% B_bw -> Bandwidth difference between b0 and b_rf
% Nf [dB] -> Noise factor due to receiver


%Variables

function [cheapest_cost,eq_ref] = techtest(d,requested_debit)


index=0;





T=23;
V=15;                                %Visibility in Km
N_fd=4;                          %Number of foggy days
D_f=3.4;                            %Duration of fog
p_time=0.1;                           %percentage of time (independent variable)
p_one=(N_fd/365.25)*(D_f/24)*100;   %percentage of time (in a year) in which the visibility does not exceed 1 km
V=(1/p_one)*p_time;                 %Visibility in Km
d_a=1;
Asis=3;
rain_frequency =readtable ('rain_k_alfa.dat');
precipitation = readtable ('precipitation.dat');

FH_equipment =readtable ('FH.dat');
nr_eq_FH=size(FH_equipment);
FSO_equipment =readtable ('FSO.dat');
nr_eq_FSO=size(FSO_equipment);
FO_equipment =readtable ('FO.dat');
nr_eq_FO=size(FO_equipment);

total_nr_eq=nr_eq_FH(1,1)+nr_eq_FSO(1,1)+nr_eq_FO(1,1);

eq_cost=ones(total_nr_eq,1);
eq_cost=eq_cost*(inf);

obs_above=1;                   %If obstacle is above the line of sight choose '1', if not, choose '-1'
h_obs_los=2;                   %height difference between tip of obstacle and line of sight in meters;





% Earth's Atmosphere
gamma0_o=0.00613408;      % pag 58 livro STVR
gamma0_w=0.000612794;     % pag 58 livro STVR
Aabs=(gamma0_o+gamma0_w)*d;
zone='K'; %Choose between "H" and "K"
if zone == 'H'
    z=1;
elseif zone == 'K'
    z=2;
end
prec_perc=0.003; %




%Comecar o for

for i=1:nr_eq_FH(1,1) 
    index=index+1;
    
    str=(table2array(FH_equipment(i,1)));
    name=str{1}; 
    f=table2array(FH_equipment(i,2));                        %choose between (1,2,4,6,7,8,10,12,15,20,25,30GHz)
    debit=table2array(FH_equipment(i,3));                    %in Mbps
    PTx=table2array(FH_equipment(i,4));                      %Transmitter Power
    GTx=table2array(FH_equipment(i,5));
    GRx=table2array(FH_equipment(i,6));                       
    SRx=table2array(FH_equipment(i,7));                      %Receiver Sensitivity
    cost=table2array(FH_equipment(i,9));

    lambda=(3*10^8)/(f*10^9);

    Ao=92.4+20*log10(d)+20*log10(f);

    %Rain attenuation
    faux=f;
    if f>30
        faux=30;
    end
    if f==13
        faux=12;
    end
    rows_rain_frequency = rain_frequency(rain_frequency.Frequency==faux,:); %loads k and alpha for assigned frequency
    k_r=(table2array(rows_rain_frequency(1,2))+(table2array(rows_rain_frequency(1,3))))/2;
    alpha_r=(table2array(rows_rain_frequency(1,4))+(table2array(rows_rain_frequency(1,5))))/2;

    rows_precipitation = precipitation(precipitation.Percentage==0.01,:); %loads Ri_001 according to the zone
    Ri_001=table2array(rows_precipitation(1,z));

    if Ri_001>100
        d0=100;
    else
        d0=35*exp(-0.015*Ri_001);
    end

    gamma_r= (k_r * (Ri_001^(alpha_r))); 
    def=d/(1+(d/d0));

    Ar_001=gamma_r*def;
    Ar_p= (Ar_001)* 0.12*prec_perc^(-0.546-0.043*log10(prec_perc)); %attenuation for given percentage

    %Obstacle attenuation

    r1e=17.32*sqrt(d/(4*f));       %radius of the first Fresnell zone in meters
    u=(obs_above)*h_obs_los*(sqrt(2)/r1e);
    Aobs=6.9+20*log10(sqrt(((u-0.1)^2)+1)+u-0.1);

    Aobs=0; %NO OBSTACLE


    % signal-to-noise ratio
    PRx=PTx+GTx+GRx-Asis-Ao-Ar_p-Aobs;          %Power detected at the receiver

    M=PRx-SRx;

    % ver pag 210

    CN_table=readtable('CN_QAM.dat');
    QAM=1024;
    CN = table2array(CN_table(CN_table.QAM==QAM,2));
    s = table2array(CN_table(CN_table.QAM==QAM,3));
    B_bw=0.3;
    b0=(debit*10^6)/log2(QAM);
    b_rf=(1+B_bw)*b0;
    N0=-204+10*log10(b_rf);
    CN_ipc=PRx-N0;


    sesr=0.00016;
    bber=0.000008;


    Kn=(1.4e-8)*f*d^(3.5);
    ms=8000/s;
    Ms=10*log10(ms);

    %SESR
    mr_sesr=Kn/sesr;
    Mr_sesr=10*log10(mr_sesr);

    mu_sesr=(mr_sesr*ms)/(ms-mr_sesr);
    Mu_sesr=10*log10(mu_sesr);

    %BBER
    alfa1=20;
    alfa2=5;
    alfa3=1;
    Nb=22500;
    rber=1e-12;
    Prber=1e-10;
    bersesr=1e-5;
    slbber=abs((log10(rber)-log10(bersesr))/(log10(Prber)-log(sesr)));

    sesr_aux=(bber-(Nb*rber)/alfa3)*(2.8*alfa2*(slbber-1)/alfa1);

    mr_bber=Kn/sesr_aux;
    Mr_bber=10*log10(mr_bber);

    mu_bber=(mr_bber*ms)/(ms-mr_bber);
    Mu_bber=10*log10(mu_bber);


     SNRmin_sesr=CN_ipc-Mu_sesr;
     SNRmin_bber=CN_ipc-Mu_bber;

     if SNRmin_bber<SNRmin_sesr
        SNRmin=SNRmin_bber;
     else
        SNRmin=SNRmin_sesr;
     end
     
     if ((SNRmin>0) && (isreal(SNRmin)) && (M>0) && (debit>=requested_debit))
%         X = sprintf('O equipamento FH %s pode ser usado. SNR= %d dB e Margem= %d dB',name,SNRmin,M);
%         disp(X)

          eq_cost(index,1)=cost;

%         if(cost<cheapest_cost)
%             cheapest_cost=cost;
%             eq_ref=index;
%             
%         end
        
%         possible_equipment(eq_aux,1)=name;
%         possible_equipment(eq_aux,2)='FH'
%         possible_equipment(eq_aux,3)=debit;
%         possible_equipment(eq_aux,4)=table2array(FH_equipment(i,9));
%         possible_equipment(eq_aux,5)=0;
%         eq_aux=eq_aux+1;
        
%      else
%         X = sprintf('O equipamento FH %s nao pode ser usado.',name);
%         disp(X)
          tgt=1;
     end
 
end
 
%% Free Space Optics

FSO_equipment =readtable ('FSO.dat');
nr_eq_FSO=size(FSO_equipment);


eff=1;
Asis=3;                                 
T=23;                                   % in Celsius Degrees
P=1018;                                 
ha=30;
c=3e8;
h_rel=0.85;
h_abs=h_rel*(-0.74+90.96*exp(T/13.67)-85.4*exp(T/13.52));
w=h_abs*d*1e-3;                         %amount of water that exists in the atmosphere along distance 'd'
Asis=3;


%Variaveis Ivo


if V<0.5
    q=0;
elseif V<1
    q=V-0.5;
elseif V<6
    q=0.16*V+0.34;
elseif V<50
    q=1.3;
else
    q=1.6;
end

%Rain Attenuation
zone='K'; %Choose between "H" and "K"

if zone == 'H'
    z=1;
elseif zone == 'K'
    z=2;
end
precipitation = readtable ('precipitation.dat');
rows_precipitation = precipitation(precipitation.Percentage==0.01,:); %loads Ri_001 according to the zone
Ri_001=table2array(rows_precipitation(1,z));

Arain=0.05556+0.00848*Ri_001-3.66e-5*Ri_001^2*d*10*log10(exp(1));



for i=1:nr_eq_FSO(1,1)
        index=index+1;
   
        str=(table2array(FSO_equipment(i,1)));
        name=str{1};
        lambda_nm=table2array(FSO_equipment(i,2));          %in nm
        lambda=lambda_nm*1e-9;                               %in mW
        debit=table2array(FSO_equipment(i,3));              %in Mbps
        PTx_mW=table2array(FSO_equipment(i,4));             %in mW
        PTx=10*log10(PTx_mW*1e-3);                          %in dBW
        GTx=table2array(FSO_equipment(i,5));                %in dB
        GRx=table2array(FSO_equipment(i,6));                %in dB
        SRx=table2array(FSO_equipment(i,7));
        cost=table2array(FSO_equipment(i,9));
        
   % if debit>=requested_debit
        
        opt_freq=c/lambda;
        B=debit*1e6;

        Agl=20*log10(4*pi*d)-20*log10(lambda*1e6);      %Geometric losses
        % Atmosphere Model
        if lambda == 785e-9
            Ai=0.0305;
            Ki=0.8;
            beta_i=0.112;
            Wi=54;
        else 
            Ai=0.211;
            Ki=0.802;
            beta_i=0.111;
            Wi=1.1;
        end
        if (w-Wi)>0             %alternative to using step function (w-Wi)
            flag=1;
        else
            flag=0;
        end
        
        gamma_abs=exp(-Ai*w^(1/2))*(1-flag)+Ki*((Wi/w)^beta_i)*flag;
        Aabs=gamma_abs*d;               %eq 3.9

        %Parte do Ivo

        c_one=(3.91/V)*(550/lambda_nm)^q;
        c_two=0.00258;
        Afog=-10*log10(exp(-(c_one+c_two*(lambda_nm^-4)*d)));

        %No Snow
        %Turbulence

        k=(2*pi)/lambda;
        C_nsq=9.8583e-18+4.9877e-16*exp(-ha/300)+2.9228e-16*exp(-ha/1200);
        cint_var=1.23*((k)^(7/6))*C_nsq*((d*1000)^(11/6));
        Aturb=sqrt(cint_var)*2;

        PRx=PTx+GTx+GRx-Asis-Agl-Aabs-Afog-Arain-Aturb;
        M=PRx-SRx;

        SNR0=sqrt((eff*(10^((PRx+Aturb)/10)))/(B*opt_freq*2*6.63e-34));

        SNR_av=SNR0/(sqrt((10^((Aturb)/10)))+10^(sqrt(cint_var)/10)*SNR0^2);

        SNR_dB=10*log10(SNR0);
    %     SNR_av_dB=10*log10(SNR_av);

        BER=0.5*erfc(0.5*sqrt((SNR0)/2));

        if ((BER<1e-6) && (M>0) && (debit>=requested_debit))
%             X = sprintf('O equipamento FSO %s pode ser usado. SNR= %d dB e Margem= %d dB',name,SNR_dB,M);
%             disp(X)
%             possible_equipment(eq_aux,1)=name;
%             possible_equipment(eq_aux,2)='FSO'
%             possible_equipment(eq_aux,3)=debit;
%             possible_equipment(eq_aux,4)=table2array(FH_equipment(i,9));
%             possible_equipment(eq_aux,5)=0;
%             eq_aux=eq_aux+1;
%             if(cost<cheapest_cost)
%             cheapest_cost=cost;
%             eq_ref=index;
%             end
        eq_cost(index,1)=cost;
        else
%             X = sprintf('O equipamento FSO %s nao pode ser usado.',name);
%             disp(X)
        tgt=1;
        end
%     else
%         X = sprintf('O equipamento FSO %s nao pode ser usado pois nao tem capacidade suficiente para a ligacao pedida.',name);
%         disp(X)
%     end
   
end

%% Fiber Optics



for i=1:nr_eq_FO(1,1)
    index=index+1;
    str=(table2array(FO_equipment(i,1)));
    name=str{1};
    debit=table2array(FO_equipment(i,2));              %in Mbps
    costperkm=(table2array(FO_equipment(i,3)));
    
    if (requested_debit<=debit)
        %if (cost<cheapest_cost)

            eq_cost(index,1)=((costperkm*d)+35000);
            
%             cheapest_cost=cost;
%             eq_ref=index;
        %end
    end

    [cheapest_cost,eq_ref] = min(eq_cost);
    
end
end
 