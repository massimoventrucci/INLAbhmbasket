
%%%%%% TO DO

\name{sim_basket}
\alias{sim_basket}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bla%%  ~~function to do ... ~~
}
\description{
MANUAL IN PROGRESS; Simulate a basket trial under several study design choices and priors on the arm-level variance}
\usage{
sim_basket(nsim, m, N, p_true, p_null, p_target, ia1_fraction = 0.4, step = 0.5, futility_threshold = 0.05, efficacy_threshold = 0.9, prior, parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{
number of simulated trials
}
  \item{m}{
number of arms
}
  \item{N}{
max number of patients enrolled per arm  (same for each arm)
}
  \item{p_true}{
a vector of length \texttt{m} with response rate for generating the data in each arm
}
  \item{p_null}{
the null hypothesis uninteresting threshold (efficacy rate under H0)
}
  \item{p_target}{
the alternative hypothesis target threshold (efficacy rate under H1)
}
  \item{ia1_fraction}{
a fraction of N at which the 1st interim analysis (IA) takes place; default 0.4
}
  \item{step}{
define subsequent IAs' timing as an increment of the patients enrolled at the 1st IA; default 0.5
}
  \item{futility_threshold}{
response rate at which the arm can be closed based on futility; default 0.05
}
  \item{efficacy_threshold}{
response rate at which the arm can be closed based on efficacy; default 0.9
}
  \item{prior}{
prior distribution for the arm-level variance, choices are:
\begin{itemize}
\item Gamma(a, b) on the precision;
\item half-t(gamma, ni) on the standard deviation (\texttt{half-t}); \item uniform (a, b) on the standard deviation;
\item PC(sd_x) prior (\texttt{PC});
\item EPC (please use the EPC_sd function to compute the desired sd_x value, TODO...);
\end{itemize}
}
\item{parameters}{
parameters of the prior distrinution on the arm-level variance}
}
\details{
bla%%  ~~ If necessary, more details than the description above ~~
}
\value{
bla%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
bla%% ~put references to the literature/web site here ~
}
\author{
bla%%  ~~who you are~~
}
\note{
bla%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
bla%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
