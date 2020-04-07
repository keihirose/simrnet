\name{nagano2019}
\alias{nagano2019}
\title{Data}
\description{
 The \code{nagano2019} data frame is a empirical dataset about wild \emph{Arabidopsis halleri} subsp. \emph{gemmifera}.  
}
\usage{data(nagano2019)
}

\format{
  A data frame.
  \describe{
    \item{X}{
    Design matrix of transcriptome data and has 835 observations (rows) and 100 variables (columns).
    }
    \item{Y}{
    Response vector containing 835 observations.
    }
  }
}

\details{
The \code{nagano2019} data frame is a part of the data which are described in Nagano, A. J. \emph{et al.} (2019). \code{X} is a data frame of transcriptome data and \code{Y} is a data frame of temperature data. Although the original \code{X} has 835 rows and 17205 columns, we only use the first 100 columns to reduce the computational loads; we need several days to complete the simulation with the original data.
}
\source{
Nagano, A. J. \emph{et al.} (2019). Annual transcriptome dynamics in natural environments reveals plant seasonal adaption. \emph{Nature plants}, \bold{5}(1), 74-83.
}
\references{
Okinaga, Y., Kyogoku, D., Kondo, S., Nagano, J. A. and Hirose, K.  The Accuracy of Lasso Estimation of Scale-Free Network, \emph{manuscript.}
}
\keyword{datasets}
