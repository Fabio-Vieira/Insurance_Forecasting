# Forecasting Insurance Reserves

I have used a Bayesian State-Space model, assuming the errors follow a student's t distribution, to forecast future reserves of an insurance portfolio. The data comes from Choy (2016), I have used it in my master dissertation. The original data was a 18 x 18 triangle, I, however, extracted a 9 x 9 matrix in order to have out of sample data to compare with the predictions of my model. Also, the estimation is done via Markov chain Monte Carlo methods. In the case of the degrees of freedom of the student's t I have used the Shaby and Wells (2010): the metropolis-hastings with log adaptive proposal. The code ends with a calculation of DIC, because in my work I have compared this models with other two, using the same data.

# References

```
Choy, S. B., Chan, J. S., & Makov, U. E. (2016). Robust bayesian analysis of loss reserving data using scale mixtures distributions. Journal of Applied Statistics , 43 (3), 396411.

Shaby, B and  Wells, M. T. Exploring an adaptive metropolis algorithm. Currently under review, 1:117, 2010.
