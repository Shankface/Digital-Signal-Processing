% Ayden Shankman -  ECE310

fprintf("Upsampling Duration: \n"); 
y = srconvert(wavread("Wagner.wav"));
audiowrite("output.wav", y, 24000);

fprintf("Test Duration: \n");
y = srconvert([1 zeros(1,3000)]);
verify(y);

%%
function out = srconvert(in)
    %upF = 24000 and lowF = 11025;
    %[L,M] = rat(upF/lowF) = [147, 320]
  
    L = 320;
    M = 147;
    
    tic
    out = upsample(in, L);
    
    lpf = designfilt('lowpassfir','PassbandFrequency',1/320, ...
         'StopbandFrequency',1.2/320,'PassbandRipple',0.08, ...
         'StopbandAttenuation',72,'DesignMethod','kaiserwin');
    
    out = fftfilt(lpf, out);
    out = M*(downsample(out,M));
    toc

end