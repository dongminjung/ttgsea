\name{bi_gru}
\alias{bi_gru}

\title{
Bidirectional GRU with embedding layer
}
\description{
A predefined function that is used as a model in "ttgsea". This is a simple model, but you can define your own model. The loss function is "mean_squared_error" and the optimizer is "adam". Pearson correlation is used as a metric.
}
\usage{
bi_gru(num_tokens, embedding_dim, length_seq, num_units)
}

\arguments{
  \item{num_tokens}{
maximum number of tokens
}
  \item{embedding_dim}{
a non-negative integer for dimension of the dense embedding
}
  \item{length_seq}{
length of input sequences, input length of "layer_embedding""
}
  \item{num_units}{
dimensionality of the output space in the GRU layer
}
}


\value{
model
}


\author{
Dongmin Jung
}


\seealso{
keras::keras_model_sequentia, keras::layer_embedding, keras::layer_gru,
keras::bidirectional, keras::layer_dense, keras::compile
}


\examples{
library(reticulate)
if (keras::is_keras_available() & reticulate::py_available()) {
  num_tokens <- 1000
  length_seq <- 30
  embedding_dim <- 50
  num_units <- 32
  model <- bi_gru(num_tokens, embedding_dim, length_seq, num_units)
  
  # stacked gru
  num_units_1 <- 32
  num_units_2 <- 16
  stacked_gru <- function(num_tokens, embedding_dim, length_seq,
                          num_units_1, num_units_2)
  {
    model <- keras::keras_model_sequential() \%>\%
      keras::layer_embedding(input_dim = num_tokens,
                             output_dim = embedding_dim,
                             input_length = length_seq,
                             mask_zero = TRUE) \%>\%
      keras::layer_gru(units = num_units_1,
                       activation = "relu",
                       return_sequences = TRUE) \%>\%
      keras::layer_gru(units = num_units_2,
                       activation = "relu") \%>\%
      keras::layer_dense(1)
      
    model \%>\%
      keras::compile(loss = "mean_squared_error",
                     optimizer = "adam",
                     metrics = custom_metric("pearson_correlation",
                                             metric_pearson_correlation))
  }
}
}
