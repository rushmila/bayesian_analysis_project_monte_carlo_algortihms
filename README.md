Assessment 2 

Monte Carlo algorithms project: DNA methylation of cytosines

In biology, epigenetics is the study of heritable phenotype changes that do not involve alterations in the DNA sequence. One important epigenetic mechanism is DNA methylation of cytosines. This mechanism has been shown to be related to complex diseases, such as heart diseases, schizophrenia and different cancers, for example at early stages of cancer. For this reason, using large-scale sequencing to detect genome-wide methylation can be used to early detection of some forms of cancer. 

The data usually come as relative frequencies of methylation. A relative frequency is a variable defined in 
[
0
,
1
]
[0,1], then a beta distributions can be considered a reasonable model. However data can come from different tissues, then we have different levels of methylation (usually categorised into three different broad categories - low, median, and high methylation), then a "mixture" of beta distributions can be more suitable. 

A reference to mixture of beta distributions to model methylation data is A Beta-mixture model for dimensionality reduction, sample classification and analysis by Laurila et al. (2011), BMC Bioinformatics.

When analysing a methylation data, it may be needed to simulate values from mixtures of beta distributions. In this assessment, you are requested to simulate values from a mixture of beta distributions using accept/reject sampling and importance sampling. 

Define 
𝑋
X a random variable describing the frequency of methylation. 
𝑋
1
,
…
,
𝑋
𝑛
X 
1
​
 ,…,X 
n
​
  are i.i.d. distributed according to 
𝑓
(
𝑥
)
f(x). Define a mixture of three beta distributions as following:
f(x)= 
3
1
​
 Be(α 
1
​
 =1,β 
1
​
 =5)+ 
3
1
​
 Be(α 
2
​
 =3,β 
2
​
 =5)+ 
3
1
​
 Be(α 
3
​
 =10,β 
3
​
 =5)

Note: 2 points will be granted if the submitted 
R
R code works (runs) with no errors, independently from its correctness.   

-  Plot the density 
𝑓
(
𝑥
)
f(x).


Implement an accept/reject algorithm, with 
𝑈
𝑛
𝑖
𝑓
(
0
,
1
)
Unif(0,1) proposal distribution.  Include the relevant code in your report. 

-  Set the seed in 
R
R at 1234.

-  Define the value 
𝐾
K, i.e. the maximum of the density and associate that with the corresponding 
𝑥
x value.

-  Simulate 10,000 values from a uniform distribution. 

-  Define the acceptance probability.

-  Accept or reject the simulated values according to the corresponding acceptance probability.


-  Compute the observed acceptance rate and compare it with the theoretical one.  Provide the expression for the theoretical acceptance rate in your report, along with the code for observed acceptance rate.  


Implement an importance sampling algorithm. Include the relevant code in your report. 

- Use the same 10,000 values generated from the 
𝑈
𝑛
𝑖
𝑓
(
0
,
1
)
Unif(0,1) for the implementation of the accept/reject algorithm. 

-  Compute the importance weights. 


-  Compare on the same plot the target density, the density of the accepted values for the accept/reject algorithm, and the density of the values weighted by importance sampling. Notice that 
R
R has an option in the 
density()
density() function where you can add a vector of weights. 


-  Which of the two algorithms seems to approximate the target distribution better? Which one you would use and why?
