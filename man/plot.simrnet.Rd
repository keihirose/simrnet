\name{plot.simrnet}
\alias{plot.simrnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box and whisker plot of RMSE
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Box and whisker plot of RMSE obtained by simrnet function.
}
\usage{
\method{plot}{simrnet}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The list of RMSE obtained by simrnet function. 
}
  \item{...}{
Additional arguments of details in boxplot.
}
%  \item{\dots}{Other graphical parameters to plot}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%}
\references{
Okinaga, Y., Kyogoku, D., Kondo, S., Nagano, J. A. and Hirose, K.  The Accuracy of Lasso Estimation of Scale-Free Network, \emph{manuscript.}
}
\author{
Kei Hirose and Yuichi Okinaga\cr
\email{mail@keihirose.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{simrnet} object.}

\examples{
data(nagano2019)
attach(nagano2019)
rho <- seq(0.1, 0.9, by = 0.1)
pars <- genpar(X,Y,rho)
result <- simrnet(pars,times.sim=5)
plot(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
