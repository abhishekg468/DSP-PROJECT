[y,fs] = audioread('audio1_90mins.wav');
y = y(:,1); %extract first array from stereo audio
dt = 1/fs; %sampling time

t_length=60*60; %length(y)*dt; %length of whole clip
t0=3; %length of each clip
z_cross=zeros(t_length/t0,1); %number of zero crossings
hf_energy=zeros(t_length/t0,1); %energy in higher frequencies
s_flux=zeros(t_length/t0,1); %Spectral flux
s_centroid=zeros(t_length/t0,1);
s_spread=zeros(t_length/t0,1);
s_flatness=zeros(t_length/t0,1);
s_stddev=zeros(t_length/t0,1);
state=zeros(t_length/t0,1);
uncert=zeros(t_length/t0,1);
final=[]; %store final audio file

P=10000; %Initial Sigma of State
R=2000; %Sigma of measurement
alpha=0.15;
beta=2;
gama=0.3;
K=0; %Kalman Gain
x=0; %Initial state
z=0; %sensor reading
A=0.9;
error=0;

for i = 0:t_length/t0-1
    clip=y(fs*t0*i+1:fs*t0*(i+1));
    z=find(diff(clip>0)~=0)+1;
    z_cross(i+1)=length(z);
    
    s=spectrogram(clip);
    s_flux(i+1)=sum(abs(FeatureSpectralFlux(s,fs)));
    s_centroid(i+1)=sum(abs(FeatureSpectralCentroid(s,fs)));
    s_rolloff(i+1)=sum(abs(FeatureSpectralRolloff(s,fs)));
    s_spread(i+1)=sum(abs(FeatureSpectralSpread(s, fs)));
    s_flatness(i+1)=sum(abs(FeatureSpectralFlatness(s,fs)));
    s_stddev(i+1)=std(FeatureSpectralFlux(s,fs));
    
    nf=2048; %number of point in DTFT
    FT = fft(clip,nf);
    
    f = fs/2*linspace(0,1,nf/2+1);
    hf_energy(i+1)=sum(abs(FT(81:1024)).^2);
    %z1=hf_energy(i+1)-40000;
    %z2=z_cross(i+1)-11225;
    %z=s_spread(i+1);
    %z=[0.25*alpha 1]*[z1;z2];
    
    z= [0.1453    0.0814    0.1462    0.2035    0.0151    0.1326    0.0814    0.1945]*[(hf_energy(i+1)-13000)/25000;
                                                                                        -(s_flux(i+1)-0.25)/0.5;
                                                                                        -(s_centroid(i+1)-6000)/20000;
                                                                                        (s_flatness(i+1)-1000)/3000;
                                                                                        (s_stddev(i+1)-0.03)/0.06;
                                                                                        (s_rolloff(i+1)-6000)/16000;
                                                                                        (s_spread(i+1)-10000)/30000;
                                                                                        (z_cross(i+1)-13000)/20000];

    
    %Measurement Update
    K=P/(P+R); %Kalman Gain
    x=x+K*(z-x);
    P=P-K*P;
    
    %Motion Update
    x=A*x;
    P=(1+gama)*P;
    
    state(i+1)=x;
    uncert(i+1)=P;
    
    time=t0*((i+1)/2);
    thresh=0.1;
    check=state(i+1)-thresh;
    flag=0;
    for ii=1:length(VarName1)
        if time>=VarName1(ii) && time<=VarName2(ii)
            flag=1;
            break;
        end
    end
    if (flag==1 && check<0) || (flag==0 && check>0)
        error=error+1;
    end
            
end

percentage_accuracy=100-(error/(t_length/t0))*100;
t_axis=0:t0:t_length-t0;
%{
figure;
plot(t_axis,hf_energy);xlabel('Seconds'); ylabel('HF Energy');
hold on;
line([VarName1 VarName1], [min(z_cross) max(z_cross)])
line([VarName2 VarName2], [min(z_cross) max(z_cross)])
figure;
plot(t_axis,z_cross);xlabel('Seconds'); ylabel('Zero Crossings');
hold on;
line([VarName1 VarName1], [min(z_cross) max(z_cross)])
line([VarName2 VarName2], [min(z_cross) max(z_cross)])
figure;
plot(t_axis,hf_energy);xlabel('Seconds'); ylabel('High Frequency Energy');
hold on;
line([VarName1 VarName1], [min(hf_energy) max(hf_energy)])
line([VarName2 VarName2], [min(hf_energy) max(hf_energy)])
figure;
plot(t_axis,s_flux);xlabel('Seconds'); ylabel('Spectral Flux');
hold on;
line([VarName1 VarName1], [min(s_flux) max(s_flux)])
line([VarName2 VarName2], [min(s_flux) max(s_flux)])
figure;
plot(t_axis,s_centroid);xlabel('Seconds'); ylabel('Spectral Centroid');
hold on;
line([VarName1 VarName1], [min(s_centroid) max(s_centroid)])
line([VarName2 VarName2], [min(s_centroid) max(s_centroid)])
figure;
plot(t_axis,s_rolloff);xlabel('Seconds'); ylabel('Spectral Roll Off');
hold on;
line([VarName1 VarName1], [min(s_rolloff) max(s_rolloff)])
line([VarName2 VarName2], [min(s_rolloff) max(s_rolloff)])
figure;
plot(t_axis,s_spread);xlabel('Seconds'); ylabel('Spectral Spread');
hold on;
line([VarName1 VarName1], [min(s_spread) max(s_spread)])
line([VarName2 VarName2], [min(s_spread) max(s_spread)])
%}
figure;
plot(t_axis,state);xlabel('Seconds'); ylabel('State');
hold on;
line([VarName1 VarName1], [min(state) max(state)])
line([VarName2 VarName2], [min(state) max(state)])
figure
plot(t_axis,uncert);xlabel('Seconds'); ylabel('Uncertainty');

%{
Thresholds
z_cross>13000 => commercial
hf_energy>13000 => commercial
s_flux<0.25 => commercial
s_centroid<6000 => commercial
s_rolloff>6000 => commercial

FINAL NORMALIZATION
(HF_energy-13000)/25000
-(s_flux-0.25)/0.5
-(s_centroid-6000)/20000
(s_flatness-1000)/3000
(s_flux_std_dev-0.03)/0.06
(s_rolloff-6000)/16000
(s_spread-10000)/30000
(z_cross-13000)/20000

accuracy=[66.67,54.5,58.08,61.25,50.833,57.33,54.5,60.75]
accuracy=[58.03,54.5,58.08,61.25,50.833,57.33,54.5,60.75]
accuracy_normalised = [0.1437    0.1175    0.1252    0.1320    0.1096    0.1236    0.1175
0.1310]
accuracy diff normalised=[0.2608    0.0704    0.1264    0.1760    0.0130
0.1147    0.0704    0.1682]
%}