outbreak in each of the countries, and provide some intuitive guesses on which direction will it go. The blue dots are real data of new cases reported; the red line is the “simulated number of new cases”, provided by the solution to the SEIR model; the green line is the “simulated number of reported new cases”, which is the solution to the SEIR model modified by the effect of insufficient testing kits.

Specifically, for the model: I used the model and parameters from Lin et al. In their paper, they simulated the SEIR (Susceptible - Exposed - Infectious - Recovered) with an additional parameter representing the strength of government/social action. The parameter modifies the transmission rate. For reference, the value of this parameter is estimated (in Lin et al.) to be around 0.42 in Wuhan during Jan. 23-29, and 0.84 afterwords. I gut guessed most countries should have a value somewhere between 0 and 0.8, that is, between do-nothing and total shut-down.

To account for the problem of under-testing, I reported the “simulated reported cased” (in green), which is computed by multiplying the number of unreported cases with an exponentially increasing rate of “testing coverage”. I used arbitrary levels 20%, 50%. 80% to represent serious-under-testing (e.g., US before UW’s kit was produced), under-testing (e.g., Wuhan in the first few days of quarantine), and a sufficient level of testing (e.g., South Korea). This modification explains some of the lagging effect of cases that are not tested yet due to limited number of test kits.

Based on these post-hoc and not-entirely-scientific assumptions, let’s have a look at what the plots say. Of course, the more action taken, the quicker the turning point appears. Italy might see the turning point soon. And I think the US needs stricter measures (in whatever forms).

This is NOT a data analysis or a statistical analysis. The parameters are chosen post-hoc. But I think some people might find it interesting.

Sources:
https://www.ijidonline.com/a…/S1201-9712(20)30117-X/fulltext
https://royalsocietypublishing.org/…/full/10…/rspb.2013.1345
