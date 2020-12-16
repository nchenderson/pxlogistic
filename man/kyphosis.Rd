\name{kyphosis}
\alias{kyphosis}
\docType{data}
\title{
Data on Children who have had Corrective Spinal Surgery
}
\description{
The \code{kyphosis} data has 81 rows and 4 columns.
This contains data on children who have had a corrective spinal operation.
}
\usage{data("kyphosis")}
\format{
  A data frame with 81 observations on the following 4 variables.
  \describe{
    \item{\code{Kyphosis}}{
             A vector of binary outcomes: 1 denotes that kyphosis was
             present after the operation; 0 denotes that kyphosis was
             absent after the operation.
          }
     \item{\code{Age}}{
            age of patient in months
       }
      \item{\code{Number}}{
           the number of vertebrae involved in the operation 
       }
       \item{\code{Start}}{
           the number of the first (topmost) vertebra operated on.
      }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
John M. Chambers and Trevor J. Hastie eds. (1992)
  \emph{Statistical Models in S}, Wadsworth and Brooks/Cole,
  Pacific Grove, CA.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data("kyphosis")
X <- model.matrix(Kyphosis ~ ., data=kyphosis)
kyph.fit <- pxlogistic(par=runif(4), X=X, y=kyphosis$Kyphosis)
kyph.fit$coef
}
\keyword{datasets}
