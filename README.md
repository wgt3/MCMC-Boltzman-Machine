#Using the Metropolis-Hastings Algorithm to analyze United Nations Voting data
This is a final for a Scientific Computing class in college. In this project I used the Metropolis-Hastings Algorithm, 
which is a Markov Chain Monte Carlo method, to sample from a Boltzmann Machine/Ising Model. The Boltzmann Machine was
fit to a United Nations voting data set using steepest ascent. Results are discussed in attached word doc.
# Basic Explantion of the Problem and the Math
REAL-WORLD PROBLEM: We want to predict how countries will vote in the United Nations. If each country is a node, then the
links between these nodes can tell us which countries tend to vote together. For each issue brought before the UN, the voting
choice is binary: yes or no. The Ising Model, orgininally conceived of to deal with magnetic moments (also binary) on a lattice, 
can be used to describe this UN voting data. By fitting the data to an Ising Model, we can examine interactions between nodes (countries) 
to find out who usually votes with who.

THE MATH: Fitting the data and sampling from the ditribution are interconnected (I assume by design). Sampling from the distribution is 
very difficult, so MCMC is necessary. The Metropolis Hastings algorithm is an MCMC method that essentially reduces the problem to sampling from a uniform distribution (which the computer does for us). The derivatives needed for steepest descent contain expected values for the MH sampler. So at each step of the steepest ascent, we are running an MH sampler and updating the weights. At the maximum of the likliehood
function we have some weights which we can examine to see the connection between nodes.
