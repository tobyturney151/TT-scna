#Set variables, load packages and clear workspace
rm(list=ls())
setwd("C:/Users/tobia/Downloads")
library(psd)

#Next steps:
#Make EACH of the plots from Fig 3.
#A Subtracted traces in dark and light for each sweep
#B Overlaid PS density in light and dark log-log plots
#C Delta PS density with Fitted Lorentzian log-log plots
#D Unitary conductance plus plateau current.

#open file
sweeps=read.csv("ChroME_cell1.csv");
#50 us per step
#200 ms lights off
#25 s lights on
#25 s lights off
#plateau current (last 10 s of each trace)
#200000/50=4001 for beginning of light (really, start at least at 4101 to allow for system lag)
#25000000/50=500000 for beginnning of dark (also give at least 100 extra for lag)
#Left-most current (should) be time, all others are current in pA.
u_cond=vector(length=4);
u_cond_num=1;
for (i in c(2,4,5,7)){
  sweep_num=i;
  time_on=sweeps[390000:490000,1];
  p_curr_on=sweeps[390000:490000,sweep_num];

  dat_on=data.frame(x=c(time_on), y=c(p_curr_on));
  base_spline_on=smooth.spline(time_on,p_curr_on,nknots=5);
  model_y_on=base_spline_on$y;

  model_dat=data.frame(x=c(time_on), y=c(model_y_on));
  #plot(y ~ x, data = model_dat, main="Model", type="l")

  #Subtract exponential fit from currrent trace
  subtracted_on=p_curr_on-model_y_on;
  subtracted_on_df=data.frame(x=c(time_on),y=subtracted_on);
  #plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l")

  #Grab power spectra
  subtracted_on_limit=subtracted_on[1:ceiling(0.9*length(subtracted_on))];
  spec_on=spectrum(subtracted_on_limit,plot=FALSE);
  norm_freq_on=20000/max(spec_on$freq)*spec_on$freq;
  PS_on_df=data.frame(x=c(norm_freq_on),y=c(spec_on$spec));
  #plot(log(y) ~ log(x), data = PS_on_df, main="Power Spectra - On")

  #Repeat for off cycle!
  time_off=sweeps[850000:950000,1];
  p_curr_off=sweeps[850000:750000,sweep_num];
  dat_off=data.frame(x=c(time_off), y=c(p_curr_off));

  base_spline_off=smooth.spline(time_off,p_curr_off,nknots=5);
  model_y_off=base_spline_off$y
  subtracted_off=p_curr_off-model_y_off;
  subtracted_off_df=data.frame(x=c(time_off),y=subtracted_off);
  #plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l")

  subtracted_off_limit=subtracted_off[1:ceiling(0.9*length(subtracted_off))];
  spec_off=spectrum(subtracted_off_limit,plot=FALSE);
  norm_freq_off=20000/max(spec_off$freq)*spec_off$freq;
  norm_spec_off=
  PS_off_df=data.frame(x=c(norm_freq_off),y=c(spec_off$spec));
  #plot(log(y) ~ log(x), data = PS_off_df, main="Power Spectra - Off",)

  #Get difference between PS_on and PS_off
  delta_PS=spec_on$spec-spec_off$spec;
  delta_PS_df=data.frame(x=c(norm_freq_off),y=c(delta_PS));
  #plot(log(abs(y)) ~ log(x), data = delta_PS_df, main="Power Spectra - Delta")


  #Fit everything between 2 and 1000 Hz in PS_delta to a Lorentzian function
  S0=delta_PS[1];
  lor_dat=data.frame(x=c(spec_off$freq), y=c(delta_PS));#[11:5000] coresponds to only fitting past 2 Hz. doesn't seem to make sense given the corner frequency is 3.8 in this example. 
  min.loren <- function(data, par) {
    with(data, sum((par[1]/(1+(x/par[2])^2) - y)^2))
  }
  (loren_result <- optim(par = c(1,1), fn = min.loren, data = lor_dat));
  model_l=loren_result$par[1]/(1+(norm_freq_off/loren_result$par[2])^2);
  model_dat=data.frame(x=c(norm_freq_off), y=c(model_l));
  #plot(y ~ x, data = model_dat, main="Lorentzian Model")

  #Get that single channel unitary conductance!!!!!!!!!!
  u_cond[u_cond_num]=loren_result$par[1]*pi*loren_result$par[2]/(-60*2*mean(p_curr_on))*333333 #Figure out units!!!!!!
  u_cond_num=1+u_cond_num;
}


mean_ucond=mean(u_cond);
sd_ucond=sd(u_cond);
  
#Plot your results!
#Set up panels
par(mfrow=c(2,2));
#par(mfrow=c(1,1));

#Example Traces of On and Off
plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l",xlim=c(min(time_on),max(time_on)),ylim=c(-30,30),xlab="Time (s)",ylab="I (pA)")
plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l",xlim=c(min(time_off),max(time_off)),ylim=c(-30,30),xlab="Time (s)",ylab="I (pA)")

#Power Spectra
plot(y ~ log(x), data = PS_on_df, main="Power Spectra",col="red",type="l",xlim=c(0,3),ylim=c(0,1.2*max(spec_on$spec)),xlab="Log(Hz)",ylab="PSD")
lines(y ~ log(x), data = PS_off_df,col="black",lty=2)
legend("topright", legend=c("Light", "Dark"),lty=1:2,col=c("red", "black"))

#Lorentzian Fit
plot(y ~ log(x), data = delta_PS_df, main="Power Spectra - Delta", col="black",type="l",xlim=c(0,3),ylim=c(0,1.2*max(delta_PS)),xlab="Log(Hz)",ylab="PSD")
lines(y ~ log(x), data = model_dat,col="red",lty=2)
legend("topright", legend=c("Experiment", "Fitted"),lty=1:2,col=c("black", "red"))

