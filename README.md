# Calcualation of Integrated Resilience Index
identifying quasi-equilibrium potential states in long-term NDVI time series

This includes the following steps:
1. Decompose the NDVI time series data into signal, noise and trend
2. It is decomposed in such a way that the number of minima and maxima in the trend data matches that of the annual NDVI data
3. Identify the maxima and minima points
4. If the first point in the time series in maxima, then it is considered as the first quasi-equiibrium point and subsequent maximas are considered as next quasi-equibrium points
5. If the first point in the time series in minima, then the same above procedure is followed
6. After identification of cycles with the help of quasi-equlibrium states, then the return rate indices from the spring-damper mechanical system is applied (Todman et al. 2016)
7. Take the mean and standard deviation of the return rate indices
8. Calculate a single composite return rate index from the return rate indices
9. Integrated resilience index =  exp(-composite return rate index)


