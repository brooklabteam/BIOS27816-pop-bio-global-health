## Building a model for COVID transmission betwee
## Paris population at large and UChicago population
## 
## By Cara Brook, 2024

rm(list=ls())                   # Clear all variables and functions


## First, write the function for your model here
## You can keep this structure and change the equations and 
## output for your own disease

# Name it whatever you want!
sir.covid.paris <- function(t,y,parms){ #this line should retain this format but you can change the name ("sir.covid.paris") if you want
  
  # The with() function gives access to the named values of parms within the
  # local environment created by the function - keep this structure as is
  with(c(as.list(y),parms),{
    # Here, you list your equations. You need a different
    # equation for every state variable in your system
    dScdt <- -BetaC*Sc*(Ip+Ic)/N #here, I divide by N (total population) to represent frequency-dependent transmission, meaning susceptibles acquire infection based on encountering infection with a hazard that scales with the proportion infected at a given time
    dIcdt <- BetaC*Sc*(Ip+Ic)/N - gamma*Ic #here susceptibles become infected and infecteds recover
    dRcdt <- gamma*Ic #recovered Chicagoans
    dSpdt <- -BetaP*Sp*Ip/N #here Parisians encounter infection based on the proportion infection at a given time. we do not count chicagoans in the population here since we think the effect is small compared to the population in Paris
    dIpdt <- BetaP*Sp*Ip/N - gamma*Ip #Parisians become infected and then recover
    dRpdt <- gamma*Ip #recovered Parisians
    # Note: Population size is constant, so don't need to specify dRdt
    return(list(c(dScdt,dIcdt,dRcdt,dSpdt,dIpdt,dRpdt)))
  })
}

N=2100000 #here is our total population size

# here, I wrote out the initial conditions for each state variable,
# which we have to feed to the model to start.
# I first wrote them as proportions of the total, such that they add
# up to 100%. below, I wrote them as total numbers - either
# approach is fine. I assumed that 1% of the total population N
# was made up of Chicagoans, so I use roughly the same distribution of S/I/R
# for Chicago and Paris (50% susceptible, 1% infected, 48% recovered)
# but I multiply the Chicago population additionally by .01 to shrink
# it. I also have to inrease the Rc class to .01*.49 (instead of .48)
# to equal 1. These are arbitrary numbers, and it does not matter much
# if there are small variations here.

inital.state.variables<- c(Sc = .5*.01*N,  
                           Ic = .01*.01*N,
                           Rc = .01*.49*N,
                           Sp = .5*N,
                           Ip = .01*N,
                           Rp = .48*N)  

# here I write them as just numbers. you could do this
# first and then just specify N = the sum of all these
# numbers if that is more intuitive.

inital.state.variables<- c(Sc = 10500,  
                           Ic = 210,
                           Rc = 10290,
                           Sp = 1050000,
                           Ip = 21000,
                           Rp = 1008000)  

#you can check that the sum here is equal to N by:
sum(inital.state.variables) #2100000

# now, make a vector of your parameters
parameters <- c(BetaC = 0.5,      # number of infectious contacts per infectious UChicago person per timestep (in days)
                gamma = 1/5,      # 1 / infectious period = 1/5 days
                N = 2100000,      # population size (constant)
                BetaP= 0.2)       # number of infectious contacts per infectious Paris person per timestep (in days)
# If you have more parameters, you will need to add lines to the above vector

# Load libary to be used for numerical integration
library(deSolve)                

## The function within deSolve that we will be using is called lsoda(). Let's
## look at the help file for this function:

?lsoda

# Now, specfiy how long to run your model for: 
# Here, we run from 0 days to 90 days (first two numbers),
# and we report the output every day (third number).
# We could report more or less frequently by decreasing or 
# increasing the third number accordingly, but the 
# model will integrate continuously and remember the outcome

time.out <- seq(0,90,1)     



# Now we run our model and store it as a dataframe
# that we call "ts.sir" for plotting (again, you could change
# the name of this dataframe if you wanted)

ts.sir <- data.frame(lsoda(
  y = inital.state.variables,   # y = the initial conditions vector you specified above
  times = time.out,             # times = the vector of how long you want to run for that you specified above
  func = sir.covid.paris,       # func = the function you want to run
  parms = parameters            # parms = the vector of parameters that you specified above
))

## Highlight above and hit "run" or "control+enter".

## Now let's look at the output:
head(ts.sir)


## Now plot it! Here are a few ways to plot it
## first, using base R and plotting all the variables on one plot:
plot(x= ts.sir$time,            # Time on the x axis
     y= ts.sir$Rp,              # Number Recovered Parisians (Rp) on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Number",           # Label the y axis
     ylim=c(0,parameters["N"]), # Here, we specify that the y axis starts at 0 and ends at the total population size - this ensures all the lines will fit on the plot
     main = "COVID in Paris",   # Plot title
     type = "l",                # Use a line plot
     col ="navy",              # Color the line navy. See here for other color options: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
     bty = "n")                 # Remove the box around the plot (stylistic only)

# You can add the other state variables one-by-one as lines of different colors
lines(ts.sir$time,  ts.sir$Ip,  col="darkred")
lines(ts.sir$time,  ts.sir$Sp,  col="darkgreen")
lines(ts.sir$time,  ts.sir$Rc,  col="royalblue")
lines(ts.sir$time,  ts.sir$Ic,  col="tomato")
lines(ts.sir$time,  ts.sir$Sc,  col="forestgreen")

# add a legend to the plot
legend(x=65,y=2100000, # xy coordinates of where to put the legend
       cex=.5, # scaling parameter for the size of the legend (it was too big at cex=1, so here I am halving it)
       col = c("navy", "darkred", "darkgreen", "royalblue", "tomato", "forestgreen"), #vector of colors of lines
       legend = c("Recovered Parisians", "Infected Parisians", "Susceptible Parisians", "Recovered Chicagoans", "Infected Chicagoans", "Susceptible Chicagoans"), #vector of what the lines represent
       lty =1)

# You can save this plot by clicking on the "Export" tab in RStudio and saving that way
# You also could type the following (without the # comment-out makring on the left)

# pdf("filenameofyourchoosing.pdf")
# put all of the text to make the plot here - so lines 112-127 above
# dev.off()

# This will save the plot to your current working directory
# You could use jpeg() or png() etc if you want a different file type
# instead of a pdf

# If you want two plots, one for Chicago and one for Paris, you could do this:

par(mfrow=c(1,2)) #this sets up the plotting window with one row and two columns

# Put Paris in the plot on the left
plot(x= ts.sir$time,            # Time on the x axis
     y= ts.sir$Rp,              # Number Recovered Parisians (Rp) on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Number",           # Label the y axis
     ylim=c(0,1.1*(inital.state.variables["Sp"])), # Here, we specify that the y axis starts at 0 and ends at 1.1*(the initial Susceptible population in Paris)
     main = "COVID in Paris",   # Plot title
     type = "l",                # Use a line plot
     col ="navy",              # Color the line navy. 
     bty = "n")                 # Remove the box around the plot (stylistic only)

# You can add the other state variables one-by-one as lines of different colors
lines(ts.sir$time,  ts.sir$Ip,  col="darkred")
lines(ts.sir$time,  ts.sir$Sp,  col="darkgreen")
legend(x=25,y=400000, # xy coordinates of where to put the legend
       cex=.5, # scaling parameter for the size of the legend (it was too big at cex=1, so here I am halving it)
       col = c("navy", "darkred", "darkgreen"), #vector of colors of lines
       legend = c("Recovered Parisians", "Infected Parisians", "Susceptible Parisians"), #vector of what the lines represent
       lty =1)


# Then, make the second plot just for UChicago
plot(x= ts.sir$time,            # Time on the x axis
     y= ts.sir$Rc,              # Number Recovered Chicagoans (Rc) on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Number",           # Label the y axis
     ylim=c(0,1.1*(inital.state.variables["Sc"])), # Here, we specify that the y axis starts at 0 and ends at 1.1*(the initial Susceptible population in Paris)
     main = "COVID in UChicago in Paris",   # Plot title
     type = "l",                # Use a line plot
     col ="royalblue",          # Color the line 
     bty = "n")                 # Remove the box around the plot (stylistic only)

lines(ts.sir$time,  ts.sir$Ic,  col="tomato")
lines(ts.sir$time,  ts.sir$Sc,  col="forestgreen")
legend(x=15,y=4000, # xy coordinates of where to put the legend
       cex=.5, # scaling parameter for the size of the legend (it was too big at cex=1, so here I am halving it)
       col = c( "royalblue", "tomato", "forestgreen"), #vector of colors of lines
       legend = c("Recovered Chicagoans", "Infected Chicagoans", "Susceptible Chicagoans"), #vector of what the lines represent
       lty =1)


dev.off()# turns off the plot


# You could also just plot the infecteds if that
# was most interesting! It is not necessary to plot all
# three state variables all the time!

# If you want to get fancy and plot in ggplot,
# you will need to reshape your data. You should be able
# to just edit my script here to do this even if you 
# do not follow all the commands:

# load libraries needed
library(reshape2)
library(ggplot2)

# Reshape data to long format:
ts.sir.long <- melt(ts.sir, id.vars = c("time"))
# Look at your data:
head(ts.sir.long)

# Plot it:
ggplot(data=ts.sir.long) + 
      theme_bw() + #sets the background black and white
      geom_line(aes(x=time, y=value, color=variable)) + #plots state variables as different colors
      ylab("number") + #label y axis
      xlab("days") # label x-axis


# If you want to divide by the population in Chicago and in Paris,
# you can add a variable to your dataframe:
ts.sir.long$pop_type <- "Paris" #sets them all to Paris
ts.sir.long$pop_type[ts.sir.long$variable=="Sc" | ts.sir.long$variable=="Ic" | ts.sir.long$variable=="Rc"] <- "UChicago" #writes over for those in states Sc, Ic, or Rc to specify they are in UChicago population

# Look at your data now:
head(ts.sir.long)
tail(ts.sir.long)

# Now plot with a "facet" to make two populations:
ggplot(data=ts.sir.long) + facet_wrap(~pop_type, scales = "free")+
  theme_bw() + #sets the background black and white
  geom_line(aes(x=time, y=value, color=variable)) + #plots state variables as different colors
  ylab("number") + #label y axis
  xlab("days") # label x-axis

## To save a ggplot, you can Export from R Studio same
## as above, or you can use the "ggsave" command
## For the latter, just type something like this after you
## make the plot (but again removing the # from the left):

# ggsave(file = "filename.png", #or .pdf or .jpeg, etc
#        units="mm", #units of plot  (can be cm or lines as well)
#        width=80,  #width in units specified
#        height=50, #heigth in units specified
#        scale=3, #scale at which it is projected in saved file - you can play with this to make the writing bigger or smaller 
#        dpi=300) #resolution of the file

## If you want to plot only part of the data (e.g. infecteds),
## you can use the "subset" function:

ggplot(data=subset(ts.sir.long, variable=="Ic" | variable=="Ip")) + facet_wrap(~pop_type, scales = "free")+
  theme_bw() + #sets the background black and white
  geom_line(aes(x=time, y=value, color=variable)) + #plots state variables as different colors
  ylab("number infected") + #label y axis
  xlab("days") # label x-axis

