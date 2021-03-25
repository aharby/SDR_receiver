
%{

    This function uses a dual-PLL to recover carrier from received signal

    Before passing passing it through the PLL it need some pre-processing,
    namely squaring then passing it to a BPF.

    When squaring the signal its bandwidth/ center frequency doubles and 
    hence if fs<2*(2*f0), the squared signal will be distorted.

    For this reason, under that condition, the signal is upsampled by
    a factor of 2.

%}
function carest=carest(rsc,f0,fs)

upsampling= false;
if fs<2*(2*f0)
    upsampling= true;
    N= length(rsc);
    Mup=2; rup=zeros(1,N*Mup); rup(1:Mup:N*Mup)=rsc;

    fs=2*fs;


    omegaC =2*f0/fs;
    f1= omegaC-0.16;
    f2= omegaC-0.15;
    f3= omegaC+0.15;
    f4= omegaC+0.16;

    fl=500; ff=[0 f1 f2 f3 f4 1];  % BPF center frequency at .4
    fa=[0 0 1 1 0 0];                  % which is twice f_0
    h=firpm(fl,ff,fa);                 % BPF design via firpm
    r=filtfilt(h,1,rup);                  % filter to give preprocessed r
else 
    r=rsc; 
    
end   
    
q=r.^2;                            % square nonlinearity

omegaC =2*(2*f0)/fs;
f1= omegaC-0.02;
f2= omegaC-0.01;
f3= omegaC+0.01;
f4= omegaC+0.02;

fl=500; ff=[0 f1 f2 f3 f4 1];  % BPF center frequency at .4
fa=[0 0 1 1 0 0];                  % which is twice f_0
h=firpm(fl,ff,fa);                 % BPF design via firpm
rpre=filtfilt(h,1,q);                  % filter to give preprocessed r


t=1/fs:1/fs:length(r)/fs;


mu1=.01; mu2=.003;                       % algorithm stepsizes
fl=100; ff=[0 .01 .02 1]; fa=[1 1 0 0];
hlpf=firpm(fl,ff,fa);
lent=length(t); th1=zeros(1,lent);       % initialize estimates
th2=zeros(1,lent); carest=zeros(1,lent);
z1=zeros(1,fl+1);                        % initialize buffer for LPF
z2=zeros(1,fl+1);                        % initialize buffer for LPF

for k=1:lent-1
    
  z1=[z1(2:fl+1), rpre(k)*sin(4*pi*f0*t(k)+2*th1(k))];
  update1=fliplr(hlpf)*z1';
  
  z2=[z2(2:fl+1), rpre(k)*sin(4*pi*f0*t(k)+2*th1(k)+2*th2(k))];
  update2=fliplr(hlpf)*z2';
  
  th1(k+1)=th1(k)-mu1*update1;           % top PLL
  th2(k+1)=th2(k)-mu2*update2;  % bottom PLL
  carest(k)=cos(2*pi*f0*t(k)+th1(k)+th2(k));       %+  0.6283         % carrier estimate
end
figure(2)
subplot(2,1,1), plot(t,th1)              % plot first theta
title('output of first PLL')
ylabel('\theta_1')
subplot(2,1,2), plot(t,th2)              % plot second theta
title('output of second PLL')
ylabel('\theta_2')

if upsampling == true
    carest=carest(1:Mup:length(carest)/2*Mup);  
end
