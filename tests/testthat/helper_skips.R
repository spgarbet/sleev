library(testthat)
library(sleev)

# should we skip tests longer than a few seconds?
SKIP_CRAN_TESTS <- TRUE

# define default values for linearMEXY
{
	rho = -.3
	p = 0.3
	hn_scale = 1
	nsieve = 20

	NSIM = 20
	n = 1000
	n2 = 400
	alpha = 0.3
	beta = 0.4
}

# generate an instance of data from default values for linearMEXY
generate_data <- function(basis = "cubic")
{
    if (!require(splines) || !require(MASS))
    {
        return(list(data=NULL, Bspline= NULL))
    }
	simX = rnorm(n)
    epsilon = rnorm(n)
    simY = alpha+beta*simX+epsilon
    error = MASS::mvrnorm(n, mu=c(0,0), Sigma=matrix(c(1, rho, rho, 1), nrow=2))

    simS = rbinom(n, 1, p)
    simU = simS*error[,2]
    simW = simS*error[,1]
    simY_tilde = simY+simW
    simX_tilde = simX+simU

    id_phase2 = sample(n, n2)

    simY[-id_phase2] = NA
    simX[-id_phase2] = NA

    # histogram basis
    if (basis == "histogram")
    {
        Bspline = matrix(NA, nrow=n, ncol=nsieve)
        cut_x_tilde = cut(simX_tilde, breaks=quantile(simX_tilde, probs=seq(0, 1, 1/nsieve)), include.lowest = TRUE)
        for (i in 1:nsieve) {
            Bspline[,i] = as.numeric(cut_x_tilde == names(table(cut_x_tilde))[i])
        }
        colnames(Bspline) = paste("bs", 1:nsieve, sep="")
    }
    # histogram basis

    # linear basis
    if (basis == "linear")
    {
        Bspline = splines::bs(simX_tilde, df=nsieve, degree=1, Boundary.knots=range(simX_tilde), intercept=TRUE)
        colnames(Bspline) = paste("bs", 1:nsieve, sep="")
    }
    # linear basis

    # quadratic basis
    if (basis == "quadratic")
    {
        Bspline = splines::bs(simX_tilde, df=nsieve, degree=2, Boundary.knots=range(simX_tilde), intercept=TRUE)
        colnames(Bspline) = paste("bs", 1:nsieve, sep="")
    }
    # quadratic basis

    # cubic basis
    if (basis == "cubic")
    {
        Bspline = splines::bs(simX_tilde, df=nsieve, degree=3, Boundary.knots=range(simX_tilde), intercept=TRUE)
        colnames(Bspline) = paste("bs", 1:nsieve, sep="")
    }
    # cubic basis

    data = data.frame(Y_tilde=simY_tilde, X_tilde=simX_tilde, Y=simY, X=simX, Bspline)

    return(list(data=data, Bspline= Bspline))
}
