\name{pps}
\alias{pps_draw}
\alias{pps_probabilities}
\alias{print.prob_pps}
\title{Sampling with probability proportional to size (pps without replacement)}
\usage{
pps_probabilities(size, n)
pps_draw(x, method = "brewer", sort = TRUE)

\method{print}{prob_pps}(x, ...)
}
\arguments{
    \item{size}{\code{[numeric vector]} measure of size.}
    \item{n}{\code{[integer]} sample size.}
    \item{x}{object of class \code{prob_pps}.}
    \item{method}{\code{[character]} currently only method \code{"brewer"}
        is implemented.}
    \item{sort}{\code{[logical]} indicating whether the sampled indices
       are sorted in ascending order (default: \code{TRUE}).}
    \item{\dots}{additional arguments.}
}
\description{
    Methods to compute the first-order sample inclusion probabilities (given
    a measure of size) and sampling mechanisms to draw samples with
    probabilities proportional to size (pps).
}
\details{
    Function \code{pps_probabilities} computes the first-order sample inclusion
    probabilities for a given sample size \code{n}; see e.g., Särndal et
    al., 1992 (p. 90). The probabilities (and additional attributes) are
    returned as a vector, more precisely as an object of class \code{prob_pps}.
    To get the probabilities, use the function \code{as.numeric()} on the
    object.

    For an object of class \code{prob_pps} (inclusion probabilities and
    additional attributes), function \code{pps_draw} draws a pps sample
    without replacement and returns the indexes of the population elements.
    Only the method of Brewer (1963, 1975) is currently implemented.
}
\value{
    Function \code{pps_probabilities} returns the probabilities (an object
    of class \code{prob_pps}).

    Function \code{pps_draw} returns a pps sample of indexes from the
    population elements.
}
\references{
Brewer, K. W. R. (1963). A Model of Systematic Sampling with Unequal
    Probabilities. \emph{Australian Journal of Statistics} \bold{5}, 93--105.
    \doi{10.1111/j.1467-842X.1963.tb00132.x}

Brewer, K. W. R. (1975). A simple procedure for \eqn{\pi}pswor,
    \emph{Australian Journal of Statistics} \bold{17}, 166--172.
    \doi{10.1111/j.1467-842X.1975.tb00954.x}

Särndal, C.-E., Swensson, B., Wretman, J. (1992).
    \emph{Model Assisted Survey Sampling}, New York: Springer-Verlag.
}
\examples{
# We are going to pretend that the workplace sample is our population.
head(workplace)

# The population size is N = 142. We want to draw a pps sample (without
# replacement) of size n = 10, where the variable employment is the measure of
# size. The first-order sample inclusion probabilities are calculated as
# follows
p <- pps_probabilities(workplace$employment, n = 10)

# Extract the probabilities
as.numeric(p)

# Now, we draw a pps sample using Brewer's method.
pps_draw(p, method = "brewer")
}
