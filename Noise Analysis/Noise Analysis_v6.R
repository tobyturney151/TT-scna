#Set variables, load packages and clear workspace
rm(list=ls())
setwd("C:/Users/tobia/Downloads")
library(psd)

#open file
sweeps=read.csv("ChroME2f_cell4.csv");

#Useful parameters to know
#50 us per step
#200 ms lights off
#25 s lights on
#25 s lights off
#plateau current (last 10 s of each trace)
#200000/50=4001 for beginning of light (really, start at least at 4101 to allow for system lag)
#25000000/50=500000 for beginnning of dark (also give at least 100 extra for lag)
#Left-most current (should) be time, all others are current in pA.

#Make empty vectors for all your data
u_cond=vector(length=ncol(sweeps)-1);
dwell_times=vector(length=ncol(sweeps)-1);
on_currs=vector(length=ncol(sweeps)-1);
off_currs=vector(length=ncol(sweeps)-1);
#u_cond=vector(length=1); #for testing
u_cond_num=1;

#Loop through your sweeps!
for (i in c(2:(length=ncol(sweeps)))){
  #for (i in c(2:2)){
  sweep_num=i;
  print(paste("Sweep: ",sweep_num-1,sep=""))
  sub_on_mat=matrix(nrow=100001,ncol=3);
  for (j in 1:3){
    start=90000+100000*j;
    end=190000+100000*j;
    time_on=sweeps[start:end,1];
    p_curr_on=sweeps[start:end,sweep_num];
    dat_on=data.frame(x=c(time_on), y=c(p_curr_on));
    
    base_spline_on=smooth.spline(time_on,p_curr_on,nknots=5);
    model_y_on=base_spline_on$y
    sub_on_mat[,j]=p_curr_on-model_y_on;
  }
  subtracted_on=rowMeans(sub_on_mat);
  subtracted_on_df=data.frame(x=c(time_on),y=subtracted_on);
  #plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l")
  
  #Grab power spectra
  subtracted_on_limit=subtracted_on[1:ceiling(0.9*length(subtracted_on))];
  spec_on=pspectrum(subtracted_on_limit, plot=FALSE,x.frqsamp=20000,niter=1);
  PS_on_df=data.frame(x=c(spec_on$freq),y=c(spec_on$spec));
  #plot(log(y) ~ log(x), data = PS_on_df, main="Power Spectra - On")
  
  #Repeat for off cycle!
  sub_off_mat=matrix(nrow=100001,ncol=3);
  for (j in 1:3){
    start=450000+100000*j;
    end=550000+100000*j;
    time_off=sweeps[start:end,1];
    p_curr_off=sweeps[start:end,sweep_num];
    dat_off=data.frame(x=c(time_off), y=c(p_curr_off));
    
    base_spline_off=smooth.spline(time_off,p_curr_off,nknots=5);
    model_y_off=base_spline_off$y
    sub_off_mat[,j]=p_curr_off-model_y_off;
  }
  subtracted_off=rowMeans(sub_off_mat);
  #plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l")
  subtracted_off_df=data.frame(x=c(time_off),y=subtracted_off);
  subtracted_off_limit=subtracted_off[1:ceiling(0.9*length(subtracted_off))];
  spec_off=pspectrum(subtracted_off_limit, plot=FALSE,x.frqsamp=20000,niter=1);
  PS_off_df=data.frame(x=c(spec_off$freq),y=c(spec_off$spec));
  #plot(log(y) ~ log(x), data = PS_off_df, main="Power Spectra - Off",)
  
  #Get difference between PS_on and PS_off
  delta_PS=spec_on$spec-spec_off$spec;
  delta_PS_df=data.frame(x=c(spec_off$freq),y=c(delta_PS));
  #plot(log(abs(y)) ~ log(x), data = delta_PS_df, main="Power Spectra - Delta")
  
  
  #Fit everything between 2 and 1000 Hz in PS_delta to a Lorentzian function
  S0=delta_PS[1];
  lor_dat=data.frame(x=c(spec_off$freq[1:200]), y=c(delta_PS[1:200]));#[11:5000] coresponds to only fitting past 2 Hz. doesn't seem to make sense given the corner frequency is 3.8 in this example. 
  min.loren <- function(data, par) {
    with(data, sum((par[1]/(1+(x/par[2])^2) - y)^2))
  }
  (loren_result <- optim(par = c(0.06,10), fn = min.loren, data = lor_dat));
  model_l=loren_result$par[1]/(1+(spec_off$freq/loren_result$par[2])^2);
  model_dat=data.frame(x=c(spec_off$freq), y=c(model_l));
  #plot(y ~ x, data = model_dat, main="Lorentzian Model")
  
  #Get that single channel unitary conductance!!!!!!!!!!
  u_cond[u_cond_num]=loren_result$par[1]*pi*loren_result$par[2]/(-0.06*2*(mean(p_curr_on)-mean(p_curr_off)))*1000 #1000X because you're converting pS to fS
  dwell_times[u_cond_num]=2*loren_result$par[2]/pi
  on_currs[u_cond_num]=mean(p_curr_on);
  off_currs[u_cond_num]=mean(p_curr_off);
  u_cond_num=1+u_cond_num;
  if (delta_PS[1]<0){
    print("Warning: The Dark Trace has more noise than the Light Trace. Check for Instrument Noise!")
  }
}


mean_dwell=mean(dwell_times);
sd_dwell=sd(dwell_times);
mean_ucond=mean(u_cond);
sd_ucond=sd(u_cond);
sub_curr=on_currs-off_currs;
mean_sub_currs=mean(sub_curr);
sd_sub_currs=sd(sub_curr);
mean_curr_on=mean(on_currs);
sd_curr_on=sd(on_currs);
mean_curr_off=mean(off_currs);
sd_curr_off=sd(off_currs);

#Plot your results!
#Set up panels
par(mfrow=c(2,2));
#par(mfrow=c(1,1));

#Example Traces of On and Off
plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l",xlim=c(min(time_on),min(time_on)+0.01*(max(time_on)-min(time_on))),ylim=c(-20,20),xlab="Time (s)",ylab="I (pA)")
plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l",xlim=c(min(time_off),min(time_off)+0.01*(max(time_off)-min(time_off))),ylim=c(-20,20),xlab="Time (s)",ylab="I (pA)")

#Power Spectra
plot(y ~ x, data = PS_on_df, main="Power Spectra",col="red",type="l",xlim=c(1,1000),ylim=c(tail(PS_on_df$y,n=1),2*PS_on_df$y[1]),xlab="Frequency (Hz)",ylab="PSD, (pA^2/Hz)",log="xy")
lines(y ~ x, data = PS_off_df,col="black",lty=2)
legend("topright", legend=c("Light", "Dark"),lty=1:2,col=c("red", "black"))

#Lorentzian Fit
plot(y ~ x, data = delta_PS_df, main="Power Spectra - Delta", col="black",type="l",xlim=c(1,1000),ylim=c(model_dat$y[5000],2*model_dat$y[1]),xlab="Frequency (Hz)",ylab="Delta PSD, (pA^2/Hz)",log="xy")
lines(y ~ x, data = model_dat,col="red",lty=2)
legend("bottomleft", legend=c("Experiment", "Fitted"),lty=1:2,col=c("black", "red"))

