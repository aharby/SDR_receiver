
%{
Implementatio of software defined radio receiver

Based on 
[softwarere receiver design  by C. Richard Johnson, Jr, William A.
Sethares, Andrew G. Klein]

By:
Abanoub Methany
Anas Abdulhalim

supervesion:
Prof./ Volker Strümpen

%}

%{
The work flow of the code and different processing stages of the signal 
goes as follows:

1. choose signal 
2. load sginal >> r
3. specify receiver parameters
4. anti image BPF >> x_if
5. carrier estimate and demodulation >> x_BB
6. clock recovery and interpolator >> x_down
7. matched filtering and header/ frames extractionn >> x_message, header
8. Equilizer filter :
    header ----> LMS equilizer 
    x_message -----> DD equilzer 
9. combined equilizer filter >> x_message_equ
10. decision device (quantizer) >> mprime
11.pam2letter decoder >> reconstructed_message 

%}


msg =input(['to choose message\n' ...
'mysteryA enter 1\n' ...
'mysteryB enter 2\n' ...
'mysteryC enter 3\n' ...
':: ']);

%load signal and set its specs

if msg ==1
    load("mysteryA.mat"); resultsSubDir= '\results\mysteryA\';
    srrc_length= 8;
    beta= 0.1;
    t_T= 8e-6 ;
    f_if= 1.7e6;
    fs= 800e3 ;
    

elseif msg == 2
    load("mysteryB.mat"); resultsSubDir= '\results\mysteryB\';
    srrc_length= 4;
    beta= 0.3;
    t_T= 10.2e-6 ;
    f_if= 1.95e6;
    fs= 700e3;
    
elseif msg==3
    load("mysteryC.mat");  resultsSubDir= '\results\mysteryC\';
    srrc_length= 3;
    beta= 0.5;
    t_T= 8.3e-6 ;
    f_if= 2.45e6;
    fs= 750e3;
else
    return;
end

%receiver parameters
preamble='A0Oh well whatever Nevermind';
M =fs*t_T;                                   % upsampling ratio 
f0 = f_if;
N= length(r);

%signal was sampled by a sub-Nyquist frequency
% iterate to find the signal image with lowest f0
while abs(f0-fs) < f0                     
    f0=abs(f0-fs);
end


% anit image BPF

omegaC =2*f0/fs;
f1= omegaC-0.21;
f2= omegaC-0.2;
f3= omegaC+0.2;
f4= omegaC+0.21;

fl=500; ff=[0 f1 f2 f3 f4 1];       % BPF center frequency at .4
fa=[0 0 1 1 0 0];                   % which is twice f_0
h=firpm(fl,ff,fa);                  % BPF design via firpm
x_if=filtfilt(h,1,r);                  % filter to give preprocessed r

figure(1)
subplot(2,1,1), plotspec(x_if,1/fs);
title('IF signal');



%carrier estimate and demodulation
x_BB=x_if.*carest(r',f0,fs)'; 

fl=500; ff=[0 .25 .26 1]; fa=[1 1 0 0];                  
h=firpm(fl,ff,fa);                                                                                                   
x_BB=filtfilt(h,1,x_BB);                 %Base band signal                                  

figure(1)
subplot(2,1,2), plotspec(x_BB,1/fs);
title('Base band signal');
print('IF_and_BB_signals','-dpng');


% Clock recovery

% Initialize parameters
l=srrc_length;                                  % 1/2 length of pulse shape (in symbols)                          
x_down=zeros(1,N);                              %initialize down smapled signal x_down 

%algorithm parameters
tnow=l*M+1; tau=0; 
tausave=zeros(1,N); tausave(1)=tau; i=0;                            
mu=0.1;                                                     %stepsize
delta=0.1;                                                  % dt

while tnow<N-l*M                                            % run iteration
  i=i+1;
  x_down(i)=interpsinc(x_BB,tnow+tau,l);                                   % interpolated value at tnow+tau
  x_deltap=interpsinc(x_BB,tnow+tau+delta,l);                          % get value to the right
  x_deltam=interpsinc(x_BB,tnow+tau-delta,l);                          % get value to the left
  dx=x_deltap-x_deltam;                                             % calculate numerical derivative
  qx=quantalph(x_down(i),[-3,-1,1,3]);                                  % quantize xs to nearest 4-PAM symbol
  tau=tau+mu*dx*(qx-x_down(i));                                         % alg update: DD
  tnow=tnow+M; tausave(i)=tau;                                      % save for plotting
end

% Plot results
figure(3)
subplot(2,1,1), plot(x_down(1:i-2),'b.')                                % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))                                % plot trajectory of tau
title('trajectory of tau');
ylabel('offset estimates'), xlabel('iterations')

x_down=2.7.*x_down(1:i-2);         % cut at the end of data points and upscale amplitude to fit 4-pam ranges

N= length(x_down);                  % redefine message length


% Pulse matching and preamble/frame extraction 
preamble = letters2pam(preamble);                                   % define the preamble
corrSignal = filter(flip(preamble), 1, x_down);                     % correlate preamble with the signal
mean = sum(corrSignal)/length(corrSignal);                          % get mean value
meanVector=mean*ones(1,length(corrSignal));
corrSignalSqrd = (corrSignal-meanVector).^2;                     

% Plot results
figure(4)
subplot(2,1,1); plot(1:length(corrSignal), corrSignal);             % plot correlation results 
title('correlation of preamble with the signal');
xlabel('signal data points')
subplot(2,1,2); plot(1:length(corrSignalSqrd), corrSignalSqrd);     % plot squared difference results
title('squared difference of correlation');
xlabel('signal data points')

[peak, peakIdx]=max(corrSignalSqrd);                                       % get maximum value peak 
pos=0;

peaks=zeros(1, peakIdx);                                           % initiate peaks position vector                                 
count=0;

for k=1:length(corrSignalSqrd)
    if corrSignalSqrd(k) >=peak/1.5                                 % check if the peak is over 2/3 of the max value
        if k-pos > 100                                            
            peaks(count+1)=k;
            peak=corrSignalSqrd(k);                                 % register the peak                                              
            count=count+1;
        end
    end
end                            
peaks=peaks(1:count+1);                                             % cut the excessive zeros out of peak vector
peaks(count+1)=N;                                           % define the end of the signal



%extract message and header
% header
header=x_down(peaks(1)+1-length(preamble):peaks(1));                      

%extract message frame by frame
lenMsg=length(x_down(peaks(1)+1:end))-length(preamble)*(length(peaks)-1);
x_message=zeros(1,lenMsg);                                            % initialize message array
msgIdx=1;

for k=1:length(peaks)-1                                             % iterate through peaks
    frame=x_down(peaks(k)+1:peaks(k+1)-length(header));  
    x_message(msgIdx:msgIdx+length(frame)-1)=frame;
    msgIdx = msgIdx+length(frame);
end



%Equalizer
%pass header through LMS equilizer and message through DD equilizer
%equilizer filter is the combination of both equilizers
%that seems to give best results

nTaps=7;                                                                % set the taps
muLMS=.005; delta=4;                                                % LMS stepsize and delay delta
muDD=.006;                                                          % DD stepsize
f=zeros(nTaps,1);  


for i=nTaps+1:length(header)                                             % iterate
   rr=header(i:-1:i-nTaps+1)';                                           % vector of received signal
   e=preamble(i-delta)-rr'*f;                                   % calculate error
   f=f+muLMS*e*rr;                                              % update equalizer coefficients
end
  

for i=nTaps+1:length(x_message)                                             % iterate
   rr=x_message(i:-1:i-nTaps+1)';                                           % vector of received signal
   e=quantalph(f'*rr,[-3,-1,1,3])'-f'*rr;                       % calculate error
   f=f+muDD*e*rr;                                               % update equalizer coefficients
end
 
x_message_equ=filter(f,1,x_message);                            % filter the frame with trained filter


% plot soft decisions
figure(5)
plot(x_message_equ,'b.'); 
title('sof decisions');

% decision device and symbol matching performance assessment
mprime=quantalph(x_message_equ,[-3,-1,1,3])';                             % quantize to +/-1 and +/-3 alphabet

% decode decision device output to text string
reconstructed_message=pam2letters(mprime);                           % reconstruct message

%display to console
disp(reconstructed_message);

%save results
fprintf('\nto save results in directory run saveResults in the command window\n');


    
    
