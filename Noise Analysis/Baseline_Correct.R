#Set variables, load packages and clear workspace
rm(list=ls())
setwd("C:/Users/tobia/Downloads")

#Next steps:
#Make EACH of the plots from Fig 3.
#A Subtracted traces in dark and light for each sweep
#B Overlaid PS density in light and dark log-log plots
#C Delta PS density with Fitted Lorentzian log-log plots
#D Unitary conductance plus plateau current.

#open file
sweeps=read.csv("ChroME2s_cell1.csv");
#50 us per step
#200 ms lights off
#25 s lights on
#25 s lights off
#plateau current (last 10 s of each trace)
#200000/50=4001 for beginning of light (really, start at least at 4101 to allow for system lag)
#25000000/50=500000 for beginnning of dark (also give at least 100 extra for lag)
#Left-most current (should) be time, all others are current in pA.
u_cond=vector(length=6);
u_cond_num=1;
corrected_on=matrix(nrow=100001,ncol=7);
corrected_off=matrix(nrow=100001,ncol=7);
time_on=sweeps[390000:490000,1];
time_off=sweeps[850000:950000,1];
corrected_on[,1]=time_on;
corrected_off[,1]=time_off;
for (i in c(2:6)){
  sweep_num=i;
  p_curr_on=sweeps[390000:490000,sweep_num];
  
  dat_on=data.frame(x=c(time_on), y=c(p_curr_on));
  base_spline_on=smooth.spline(time_on,p_curr_on,nknots=5);
  model_y_on=base_spline_on$y;
  
  model_dat=data.frame(x=c(time_on), y=c(model_y_on));
  #plot(y ~ x, data = model_dat, main="Model", type="l")
  
  #Subtract exponential fit from currrent trace
  subtracted_on=p_curr_on-model_y_on;
  subtracted_on_df=data.frame(x=c(time_on),y=subtracted_on);
  
  #Repeat for off cycle!
  p_curr_off=sweeps[850000:750000,sweep_num];
  dat_off=data.frame(x=c(time_off), y=c(p_curr_off));
  
  base_spline_off=smooth.spline(time_off,p_curr_off,nknots=5);
  model_y_off=base_spline_off$y
  subtracted_off=p_curr_off-model_y_off;
  subtracted_off_df=data.frame(x=c(time_off),y=subtracted_off);
  #plot(y ~ x, data = subtracted_off_df, main="Subtracted Current Trace - Off", type="l")
  
  subtracted_off_limit=subtracted_off[1:ceiling(0.9*length(subtracted_off))];
  corrected_on[,i]=subtracted_on;
  corrected_off[,i]=subtracted_off;
}
write.csv(corrected_on,"ChroME2s_cell1_corrected_on.csv")
write.csv(corrected_off,"ChroME2s_cell1_corrected_off.csv")