#Set variables, load packages and clear workspace
rm(list=ls())
setwd("C:/Users/tobia/Downloads")
library(psd)

#Next steps:
#Come through ALL the PsChR, ChroME , ChroME2s, and ChroME2f cells (in that order)
#   until you find one with a very clear SNR! Use that as your test case (best so far is PsChR_cell2)
# Very high noise in all the traces! Try to isolate traces with minimal low freq noise.
#Get stable power spectra.
#   May need to look at other PS estimating functions (spec.prgm seemed close to working, keep trying!)
#Get stable optimizations of the Lorentzian! (seems to be working with just fitting between 0 to 200 Hz)
#   May need to brush up on other optim functions.

#open file
sweeps=read.csv("2021_10_27_0005.atf");
#50 us per step
#200 ms lights off
#25 s lights on
#25 s lights off
#plateau current (last 10 s of each trace)
#200000/50=4001 for beginning of light (really, start at least at 4101 to allow for system lag)
#25000000/50=500000 for beginnning of dark (also give at least 100 extra for lag)
#Left-most current (should) be time, all others are current in pA.

#10-22-21
#S1=736.8698, 37 on, 47 off
#S2=942.5,70 on, 90 off
#S3=312.0847, 113 on, 122 off
#S4=889.7237, 147 on, 157, off
#S5

sweep_num=2;

time_on=as.numeric(sweeps[3:8192,1]);
p_curr_on=as.numeric(sweeps[3:8192,as.integer(295/1.634)]);

full_t=sweeps[,1];
full_curr=sweeps[,2];
full_dat_on=data.frame(x=c(full_t), y=c(full_curr));
par(mfrow=c(1,1));
#plot(y ~ x, data = full_dat_on, main="Full Raw Trace",type="l",xlab="Time (s)",ylab="I (pA)");

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

#plot(y ~ x, data = dat_on, main="Raw trace",type="l");
#abline(a = result$par[1], b = result$par[2], col = "red")
#model_y=result$par[1]*exp(result$par[2]*time)+result$par[3];
#model_y=results_lm$coefficients[2]*time_on+results_lm$coefficients[1]
#for (i in 1:length(model_y))
#{
#  #print("yes")
#  model_y[i]=p_curr_on[i]-base_spline$y[ceiling(i/100)];
#}

model_dat_on=data.frame(x=c(time_on), y=c(model_y_on));
#plot(y ~ x, data = model_dat_on, main="Baseline - On", type="l")

#Subtract exponential fit from currrent trace
subtracted_on=p_curr_on-model_y_on;
subtracted_on_df=data.frame(x=c(time_on),y=subtracted_on);
#plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l")

#Grab power spectra
subtracted_on_limit=subtracted_on[1:ceiling(0.9*length(subtracted_on))];
spec_on=pspectrum(subtracted_on_limit, plot=FALSE,x.frqsamp=5000,niter=1);
PS_on_df=data.frame(x=c(spec_on$freq),y=c(spec_on$spec));
#plot(log10(y) ~ log10(x), data = PS_on_df, main="Power Spectra - On, Taper = 0.01", type="l",xlim=c(0,3),ylim=c(-3,0))


#Repeat for off cycle!
time_off=as.numeric(sweeps[3:8192,1]);
p_curr_off=as.numeric(sweeps[3:8192,as.integer(313/1.634)]);
dat_off=data.frame(x=c(time_off), y=c(p_curr_off));
#plot(y ~ x, data = dat_off, main="Raw Trace - Off",type="l",xlab="Time (s)",ylab="I (pA)");

#min.expo <- function(data, par) {
#  with(data, sum((par[1]*exp(par[2] * x)+par[3] - y)^2))
#}
#(result <- optim(par = c(-1, -1,mean(p_curr)), fn = min.expo, data = dat_off));
#results_lm=lm(y ~ x, data = dat_off); 
#model_y=result$par[1]*exp(result$par[2]*time)+result$par[3];
#model_y=results_lm$coefficients[2]*time_off+results_lm$coefficients[1]
base_spline_off=smooth.spline(time_off,p_curr_off,nknots=5);
model_y_off=base_spline_off$y
model_dat_off=data.frame(x=c(time_off), y=c(model_y_off));
#plot(y ~ x, data = model_dat_off, main="Baseline - Off", type="l")
subtracted_off=p_curr_off-model_y_off;
subtracted_off_df=data.frame(x=c(time_off),y=subtracted_off);
#plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l")

subtracted_off_limit=subtracted_off[1:ceiling(0.9*length(subtracted_off))];
spec_off=pspectrum(subtracted_off_limit, plot=FALSE,x.frqsamp=5000,niter=1);
PS_off_df=data.frame(x=c(spec_off$freq),y=c(spec_off$spec));
#plot(log(y) ~ log(x), data = PS_off_df, main="Power Spectra - Off",)

#Get difference between PS_on and PS_off
delta_PS=spec_on$spec-spec_off$spec;
delta_PS_df=data.frame(x=c(spec_off$freq),y=c(delta_PS));
#plot(log(abs(y)) ~ log(x), data = delta_PS_df, main="Power Spectra - Delta")


#Fit everything between 2 and 1000 Hz in PS_delta to a Lorentzian function
S0=delta_PS[1];
lor_dat=data.frame(x=c(spec_off$freq[1:3000]), y=c(delta_PS[1:3000]));#[11:5000] corresponds to only fitting past 2 Hz. doesn't seem to make sense given the corner frequency is 3.8 in this example. 
min.loren <- function(data, par) {
  with(data, sum((par[1]/(1+(x/par[2])^2) - y)^2))
}
(loren_result <- optim(par = c(0.06,10), fn = min.loren, data = lor_dat));
model_l=loren_result$par[1]/(1+(spec_off$freq/loren_result$par[2])^2);
model_dat_lor=data.frame(x=c(spec_off$freq), y=c(model_l));
#plot(y ~ x, data = model_dat, main="Lorentzian Model")

#Get that single channel unitary conductance!!!!!!!!!!
u_cond=loren_result$par[1]*pi*loren_result$par[2]/(-0.06*2*(mean(p_curr_on)-mean(p_curr_off)))*1000 #1000X because you're converting pS to fS

if (delta_PS[1]<0){
  print("Warning: The Dark Trace has more noise than the Light Trace. Check for Instrument Noise!")
}

#Plot your results!

#Set up panels
par(mfrow=c(2,2));
#par(mfrow=c(1,1));

#Example Traces of On and Off
plot(y ~ x, data = subtracted_on_df, main="Subtracted Current Trace - On", type="l",xlim=c(min(time_on),min(time_on)+0.04*(max(time_on)-min(time_on))),ylim=c(-10,10),xlab="Time (s)",ylab="I (pA)")
plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l",xlim=c(min(time_off),min(time_off)+0.04*(max(time_off)-min(time_off))),ylim=c(-10,10),xlab="Time (s)",ylab="I (pA)")

#Power Spectra
plot(y ~ x, data = PS_on_df, main="Power Spectra",col="red",type="l",xlim=c(1,1000),ylim=c(tail(PS_on_df$y,n=1),2*PS_on_df$y[1]),xlab="Frequency (Hz)",ylab="PSD, (pA^2/Hz)",log="xy")
lines(y ~ x, data = PS_off_df,col="black",lty=2)
legend("topright", legend=c("Light", "Dark"),lty=1:2,col=c("red", "black"))

#Lorentzian Fit
plot(y ~ x, data = delta_PS_df, main="Delta Power Spectra - Light Minus Dark", col="black",type="l",xlim=c(1,1000),ylim=c(model_dat_lor$y[3000],2*model_dat_lor$y[1]),xlab="Frequency (Hz)",ylab="Delta PSD, (pA^2/Hz)",log="xy")
lines(y ~ x, data = model_dat_lor,col="red",lty=2)
legend("bottomleft", legend=c("Experiment", "Fitted"),lty=1:2,col=c("black", "red"))