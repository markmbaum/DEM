function [frequency,fdom] = discreteFFT(samplingfreq,tdom)
%discreteFFT takes the fourier transform of a dataset in the time domain,
%converting it into the frequency domain. The input data must have a sample
%frequecy of 1. The returned data is the magnitude of the fourier transform
%without the phasors, normalized to 1, without the DC offset removed, and 
%without the reflected part.
%[frequency,fdom] = discreteFFT(samplingfreq,tdom)
%
%   tdom - time domain signal
%   samplingfreq - samping frequency
%   
%   frequency - an array of the frequencies returned
%   fdom - frequency domain signal

fdom=fft(tdom); %take discrete fourier transform
%fdom(1)=0;%remove nonvarying offset/DC offset
fdom=abs(fdom);%remove phasors for magnitude only
fdom=fdom/max(fdom);%normalize
L=length(fdom);

if(mod(L,2)==0)
    fdom=fdom(1:(L/2+1)); %exclude reflection
else
    fdom=fdom(1:round(L/2));
end

frequency=samplingfreq*(0:1/L:.5);

end