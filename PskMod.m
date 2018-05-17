function y=PskMod(M,Fs,Rb,Fc,Len)
t=0:1/Fs:(Len*Fs/Rb-1)/Fs;
s=randi([0,M-1],1,Len);  
s1=pskmod(s,M);
I=real(s1);
Q=imag(s1);

rcos_Ads_i=rectpulse(I,Fs/Rb);
rcos_Ads_q=rectpulse(Q,Fs/Rb);

f0_i=cos(2*pi*Fc*t); 
f0_q=sin(2*pi*Fc*t);      
    
 
y1=rcos_Ads_i.*f0_i+rcos_Ads_q.*f0_q;  
    x1=(y1(1:Len));
y=x1./sqrt(mean(abs(x1).^2)); 
