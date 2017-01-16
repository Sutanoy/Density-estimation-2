# Density-estimation-2
Geometric framework for nonparametric conditional density estimation

The file for conditional density estimaion(cde) with 1 predictor is univariatecde.m. It calls locallinearreg.m for employing local linear regression to compute the mean function for the initial estimate
The file for cde with 2 predictors is bivariatecde.m. It calls locallinearregbivariate.m to compute bivariate local linear regression to compute initial mean function.
The file for cde with more than 2 predictors(and irrelevant predictors) is multivariatecde.m. It does not call any other function to compute the initial mean function, but computes a Nadaraya Watson estimate itself.

Each of the functions for cde calls the function formwteddensityregressionMh.m to compute the unpenalized log likelihood objective function used to obtain the conditional density estimate at a prespecified point. The point has to be specified in the functions univariatecde/bivariatecde/multivariatecde, as the case. The function formwteddensityregressionMh.m itself calls the function FormGammaFromC.m to obtain the \gamma from its corresponding tangent space representation.

laprnd.m is a file downloaded from Fileexchange online. It is used to generate samples from Laplace distribution. The author is Elvis Chen. The link to the file is https://www.mathworks.com/matlabcentral/fileexchange/13705-laplacian-random-number-generator

gestationdata.mat is the data studied in Longnecker(2001), testgest.m performs the data analysis with 1 predictor. testgestbivariate.m performs the analysis with 2 predictors. 
