function result = DAQTimeCal(freqrange, freqstep, repN, repfreq, iter)
result = num2str(freqrange/freqstep*repN/repfreq/60*iter)+" min";
end