\name{sampling_generator}
\alias{sampling_generator}

\title{
Generator function
}
\description{
This is a generator function that yields batches of training data then pass the function to the "fit_generator" function.
}


\usage{
sampling_generator(X_data, Y_data, batch_size)
}


\arguments{
  \item{X_data}{
inputs
}
  \item{Y_data}{
targets
}
  \item{batch_size}{
batch size
}
}


\value{
generator for "fit_generator"
}


\author{
Dongmin Jung
}


\examples{
X_data <- matrix(rnorm(200), ncol = 2)
Y_data <- matrix(rnorm(100), ncol = 1)
sampling_generator(X_data, Y_data, 32)
}


