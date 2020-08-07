\name{pxlogistic}
\alias{pxlogistic}
\title{
Parameter-expanded ECME algorithm for logistic regression
}
\description{
This is an implementation of a parameter-expanded ECME algorithm for estimating the regression coefficients in a logistic regression model.
}
\usage{
pxlogistic(par, X, y, n.trials=rep(1, length(y)), control=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
A vector of initial parameter values.
}
\item{X}{
The design matrix.
}
\item{y}{
The vector of outcomes.
}
\item{n.trials}{
The vector containing the "number of trials" for each outcome.
}
\item{control}{
A list of control parameters specifying any changes to default values of algorithm control parameters. Full names of control list elements must be specified, otherwise, user-specifications are ignored. See *Details*.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following components
  \item{coef}{The vector of regression coefficient estimates.}
  \item{iter}{The number of PX-ECME iterations required to reach convergence.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zhongzhe Ouyang and Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 500
beta.true <- c(1, 0, 0.3, -0.3)
X <- cbind(1, matrix(rnorm(n*3, mean=0), ncol=3))
phi <- X \%*\% beta.true
probs <- 1/(1+exp(-phi))
y <- rbinom(n, size=1, prob=probs)

beta.init <- runif(4)
logist.fit <- pxlogistic(par=beta.init, X=X, y=y) 
}
\keyword{ regression }
\keyword{ optimize }
