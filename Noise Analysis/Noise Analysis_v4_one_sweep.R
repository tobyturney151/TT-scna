#Set variables, load packages and clear workspace
rm(list=ls())
setwd("C:/Users/tobia/Downloads")
library(spectral)

#Next steps:
#Make EACH of the plots from Fig 3.
#A Subtracted traces in dark and light for each sweep
#B Overlaid PS density in light and dark log-log plots
#C Delta PS density with Fitted Lorentzian log-log plots
#D Unitary conductance plus plateau current.

#open file
sweeps=read.csv("PsChR_cell2.csv");
#50 us per step
#200 ms lights off
#25 s lights on
#25 s lights off
#plateau current (last 10 s of each trace)
#200000/50=4001 for beginning of light (really, start at least at 4101 to allow for system lag)
#25000000/50=500000 for beginnning of dark (also give at least 100 extra for lag)
#Left-most current (should) be time, all others are current in pA.

sweep_num=2;

time_on=sweeps[390000:490000,1];
p_curr_on=sweeps[390000:490000,sweep_num];

#full_t=sweeps[,1];
#full_curr=sweeps[,2];
#full_dat_on=data.frame(x=c(full_t), y=c(full_curr));
#plot(y ~ x, data = full_dat_on, main="Full Raw Trace - On",type="l",xlab="Time (s)",ylab="I (pA)");


#generate exponential fit
dat_on=data.frame(x=c(time_on), y=c(p_curr_on));
#plot(y ~ x, data = dat_on, main="Raw Trace - On",type="l");

#min.expo <- function(data, par) {
#  with(data, sum((par[1]*exp(par[2] * x)+par[3] - y)^2))
#}
#(result <- optim(par = c(-10, -1,mean(p_curr)), fn = min.expo, data = dat_on));

#results_lm=lm(y ~ x, data = dat_on); #linear model

#base_spline=splinefun(time_on,p_curr_on,method="natural");
base_spline_on=smooth.spline(time_on,p_curr_on,nknots=5);
model_y_on=base_spline_on$y;

#plot(y ~ x, data = dat_on, main="Least square regression");
#abline(a = result$par[1], b = result$par[2], col = "red")
#model_y=result$par[1]*exp(result$par[2]*time)+result$par[3];
#model_y=results_lm$coefficients[2]*time_on+results_lm$coefficients[1]
#for (i in 1:length(model_y))
#{
#  #print("yes")
#  model_y[i]=p_curr_on[i]-base_spline$y[ceiling(i/100)];
#}

model_dat=data.frame(x=c(time_on), y=c(model_y_on));
#plot(y ~ x, data = model_dat, main="Model", type="l")

#Subtract exponential fit from currrent trace
subtracted_on=p_curr_on-model_y_on;
subtracted_on_df=data.frame(x=c(time_on),y=subtracted_on);
#plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l")

#Grab power spectra
subtracted_on_limit=subtracted_on[1:ceiling(0.9*length(subtracted_on))];
#fft_on=spec.fft(y=subtracted_on,x=time_on);
#spec_on=fft_on$PSD;
spec_on=spectrum(subtracted_on_limit,plot=FALSE);
norm_freq_on=10000/max(spec_on$freq)*spec_on$freq;
norm_spec_on=spec_on$spec/sum(spec_on$spec);
#spec_on=pspectrum(subtracted_on_limit, plot=FALSE,x.frqsamp=20000,niter=10);
PS_on_df=data.frame(x=c(norm_freq_on),y=c(norm_spec_on));
#plot(log(y) ~ log(x), data = PS_on_df, main="Power Spectra - On",type="l",xlim=c(0,3),ylim=c(0,10));


#Repeat for off cycle!
time_off=sweeps[850000:950000,1];
p_curr_off=sweeps[850000:950000,sweep_num];
dat_off=data.frame(x=c(time_off), y=c(p_curr_off));

#min.expo <- function(data, par) {
#  with(data, sum((par[1]*exp(par[2] * x)+par[3] - y)^2))
#}
#(result <- optim(par = c(-1, -1,mean(p_curr)), fn = min.expo, data = dat_off));
#results_lm=lm(y ~ x, data = dat_off); 
#model_y=result$par[1]*exp(result$par[2]*time)+result$par[3];
#model_y=results_lm$coefficients[2]*time_off+results_lm$coefficients[1]
base_spline_off=smooth.spline(time_off,p_curr_off,nknots=5);
model_y_off=base_spline_off$y
subtracted_off=p_curr_off-model_y_off;
subtracted_off_df=data.frame(x=c(time_off),y=subtracted_off);
#plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l")

subtracted_off_limit=subtracted_off[1:ceiling(0.9*length(subtracted_off))];
spec_off=spectrum(subtracted_off_limit,plot=FALSE);
norm_freq_off=10000/max(spec_off$freq)*spec_off$freq;
norm_spec_off=spec_off$spec/sum(spec_off$spec);
#spec_on=pspectrum(subtracted_off_limit, plot=FALSE,x.frqsamp=20000,niter=10);
PS_off_df=data.frame(x=c(norm_freq_off),y=c(norm_spec_off));
#plot(log(y) ~ log(x), data = PS_off_df, main="Power Spectra - Off",type="l",xlim=c(0,3),ylim=c(0,10));

#Get difference between PS_on and PS_off
delta_PS=abs(norm_spec_on-norm_spec_off);
delta_PS_df=data.frame(x=c(norm_freq_off),y=c(delta_PS));
#plot(log(abs(y)) ~ log(x), data = delta_PS_df, main="Power Spectra - Delta",type="l",xlim=c(0,3),ylim=c(0,10))


#Fit everything between 2 and 1000 Hz in PS_delta to a Lorentzian function
S0=delta_PS[1];
lor_dat=data.frame(x=c(spec_off$freq[1:100]), y=c(delta_PS[1:100]));#[11:5000] coresponds to only fitting past 2 Hz. doesn't seem to make sense given the corner frequency is 3.8 in this example. 
min.loren <- function(data, par) {
  with(data, sum((par[1]/(1+(x/par[2])^2) - y)^2))
}
(loren_result <- optim(par = c(1,1), fn = min.loren, data = lor_dat));
model_l=loren_result$par[1]/(1+(norm_freq_off/loren_result$par[2])^2);
model_dat=data.frame(x=c(norm_freq_off), y=c(model_l));
#plot(y ~ x, data = model_dat, main="Lorentzian Model")

#Get that single channel unitary conductance!!!!!!!!!!
u_cond=loren_result$par[1]*pi*loren_result$par[2]/(-60*2*mean(p_curr_on))*333333 #Figure out units!!!!!!

#Plot your results!

#Set up panels
par(mfrow=c(2,2));
#par(mfrow=c(1,1));

#Example Traces of On and Off
plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l",xlim=c(min(time_on),min(time_on)+0.01*(max(time_on)-min(time_on))),ylim=c(-30,30),xlab="Time (s)",ylab="I (pA)")
plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l",xlim=c(min(time_off),min(time_off)+0.01*(max(time_off)-min(time_off))),ylim=c(-30,30),xlab="Time (s)",ylab="I (pA)")

#Power Spectra
plot(log(y) ~ log(x), data = PS_on_df, main="Power Spectra",col="red",type="l",xlim=c(0,3),ylim=c(-10,0),xlab="Log(Hz)",ylab="PSD")
lines(log(y) ~ log(x), data = PS_off_df,col="black",lty=2)
legend("topright", legend=c("Light", "Dark"),lty=1:2,col=c("red", "black"))

#Lorentzian Fit
plot(log(y) ~ log(x), data = delta_PS_df, main="Power Spectra - Delta", col="black",type="l",xlim=c(0,3),ylim=c(-10,0),xlab="Log(Hz)",ylab="PSD")
lines(log(y) ~ log(x), data = model_dat,col="red",lty=2)
legend("topright", legend=c("Experiment", "Fitted"),lty=1:2,col=c("black", "red"))
