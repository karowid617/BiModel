% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BernoulliEM.R
\name{BernoulliEM}
\alias{BernoulliEM}
\title{Function to perform binomial mixture model (BMM)}
\usage{
BernoulliEM(
  X,
  K = 2,
  start_ini = 20,
  ini = "random",
  m_iter = 3000,
  eps = 1e-40
)
}
\arguments{
\item{X}{Matrix of binary data to cluster where in rows there are features and in columns observations.}

\item{K}{Number of clusters to divide your data. Default K=2.}

\item{start_ini}{Number of starting initialization.}

\item{ini}{Method for parameters initialization. There are three options 'random' (default), 'kmeans' or "kmeanspp".}

\item{m_iter}{Maximum number of iteration for EM algorithm. Default value 3000.}

\item{eps}{Minimum delta of model LL to stop EM algorithm.Default value 1e-40.}
}
\value{
Function returns a \code{list} which contains: \describe{
 \item{probs}{Matrix of finall probabilities for each observation and cluster.}
 \item{alphas}{Finall vector of each cluster weight in Bernoulli mixture model.}
 \item{clusters}{Cluster assigment of each observation.}
 \item{iter}{Number of EM interations.}
 \item{delta}{Final EM step delta fo algorithm stop.}
 \item{bic}{Bayesian information criterion (BIC) value for fitted model.}
 \item{ll}{Log-Likelihood of fitted model in clustering.}
}
}
\description{
Function description
}
\examples{
\dontrun{
data(example)
res<-BernoulliEM(example$nouli_data, 4, start_ini = 20, ini = "random", m_iter=3000, eps=1e-40)
}

}
\seealso{
\code{\link{bernoulliEM_ini}}
}
