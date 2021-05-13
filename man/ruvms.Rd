\name{ruvms}
\alias{ruvms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Removing unwanted variations in mass spectrometry data
}
\description{
The ruvms algorithm that removes unwanted variations from high dimensional, omics style data such as mass spectrometry data.
}
\usage{
ruvms(Y, M, ctl, k1, lambda, ..., inputcheck)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The input raw data for normalisation. A m by n matrix, where m is the number of samples and n is the number of measurements (e.g., genes, proteins, metabolites or lipids).
}
  \item{M}{
The replication structure. See RUVIII from the ruv R package for details.
}
  \item{ctl}{
A binary vector to sepecify the negative control measurements. If the j-th measurement is a negative control feature, then the j-th element of this vector is TRUE, and FALSE otherwise.
}
  \item{k1}{
A tuning parameter. Default is set to 0.
}
  \item{lambda}{
The regularisation parameter. The default is set to 10^-5.
}
  \item{eta}{
Measurement-wise (as opposed to sample-wise) covariates. See RUVIII and RUVI for details.
}
  \item{include.intercept}{
A required input parameter of RUVI. See RUVI for full details.
}
  \item{average}{
Logical; whether to average replications after normalisation.
}
  \item{return.info}{
Logical; if TRUE, the sepctural decomposition results on the normalisation matrix will also be returned. If FALSE, only normalised data matrix will be returned.
}
  \item{inputcheck}{
Logical; a basic check on the data inputs.
}

}
\details{

}
\value{
If return.info = FALSE, only normalised data matrix is returned. If return.info = TRUE, a list of outputs are returned:

   \item{newY}{The normalised data matrix. }
   \item{di}{Eigenvalues of the normalisation matrix N, ordered decreasingly. }
   \item{hi}{Eigenvectors of the normalisation matrix N. }

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mengbo Li <mengbo.li@outlook.com>
}

\seealso{
 \code{\link[ruv]{RUVIII}};
 \code{\link[ruv]{RUVI}}.
}

\examples{

normalised <- ruvms(Y = rawData, M, ctl)
