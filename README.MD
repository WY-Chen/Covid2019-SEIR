Covid2019. This small project provides some intuitive guesses on which direction the outbreak will go. 
The data is from [GuangchuangYu/nCov2019](https://github.com/GuangchuangYu/nCov2019), and the model is from [Lin et al.](https://www.ijidonline.com/a…/S1201-9712(20)30117-X/fulltext). 

The main method: 
```
sim_SEIR<-function(
  region="South Korea", # province in china, or country/region of other parts of the world
  N0=51.47e6            # population of region 
  alpha=c(0,0.85),      # vector of length(breakpt)+1: levels of gov reaction, in range (0,1). 
  breakpt=c(41),        # vector of integers: indices of the change of policy
  F0=0,                 # number of zootonic cases
  I0=1,                 # initial # of cases
  offset=0,             # offset to SEIR model, use if we believe the transmission started before first reported data
  window=100,           # "forcast" window: date scale of the plot
  mu=0,                 # rate of emmigration, usually 0. 
  report_rate=c(0.3,1),    # vector of length(report_breakpt)+1: testing capacity, 0=no testing, 1=all cases gets tested
  report_breakpt=c(41),    # vector: indices of the chagne of report rate 
  explanation=NULL         # string: printed on the figure. 
)
```

For example
```
sim_SEIR("Hubei",N0=55e6,alpha=c(0,0.42,0.84),
         breakpt = c(54,64),window=100,
         offset=0,I0=5,F0=10,mu=0.02,
         report_rate = c(0.02,0.3,1),
         report_breakpt = c(54,77),
         explanation = 
           "Gov. Action Strength change on Jan. 23 & Feb 3 (level=0.42 & 0.84), 
         testing capacity change on Jan. 23 & Feb 16 (propotion=1%,30%,100%),
         emigration rate=0.02%")
```
This gives

![alt text](figures/Hubei.png)

The blue dots are real data of new cases reported; the red line is the “simulated number of new cases”, provided by the solution to the SEIR model; the green line is the “simulated number of reported new cases”, which is the solution to the SEIR model modified by the effect of insufficient testing kits.

Specifically, for the model: I used the model and parameters from Lin et al. In their paper, they simulated the SEIR (Susceptible - Exposed - Infectious - Recovered) with an additional parameter representing the strength of government/social action. The parameter modifies the transmission rate. For reference, the value of this parameter is estimated (in Lin et al.) to be around 0.42 in Wuhan during Jan. 23-29, and 0.84 afterwords. I guess most countries should have a value somewhere between 0 and 0.8, that is, between do-nothing and total shut-down.

To account for the problem of under-testing, I reported the “simulated reported cased” (in green), which is computed by multiplying the number of unreported cases with an exponentially increasing rate of “testing coverage”. I used arbitrary levels 20%, 50%. 80% to represent serious-under-testing (e.g., US before UW’s kit was produced), under-testing (e.g., Wuhan in the first few days of quarantine), and a sufficient level of testing (e.g., South Korea). This modification explains some of the lagging effect of cases that are not tested yet due to limited number of test kits.

Based on these post-hoc and not-entirely-scientific assumptions, let’s have a look at what the plots say. Of course, the more action taken, the quicker the turning point appears. Italy might see the turning point soon. And I think the US needs stricter measures (in whatever forms).

This is NOT a data analysis or a statistical analysis. The parameters are chosen post-hoc. But I think some people might find it interesting.

Other Sources:
https://rpubs.com/srijana/110753
http://www.public.asu.edu/~hnesse/classes/seir.html
https://royalsocietypublishing.org/…/full/10…/rspb.2013.1345
