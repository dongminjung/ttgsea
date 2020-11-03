\name{metric_pearson_correlation}
\alias{metric_pearson_correlation}

\title{
Pearson correlation coefficient
}
\description{
Pearson correlation coefficient can be seen as one of the model performance metrics. This is a measure of how close the predicted value is to the true value. If it is close to 1, the model is considered a good fit. If it is close to 0, the model is not good. A value of 0 corresponds to a random prediction.
}



\author{
Dongmin Jung
}


\seealso{
keras::k_mean, keras::sum, keras::k_square, keras::k_sqrt
}


\examples{
library(reticulate)
if (keras::is_keras_available() & reticulate::py_available()) {
  num_tokens <- 1000
  length_token <- 30
  embedding_dims <- 50
  num_units_1 <- 32
  num_units_2 <- 16
  
  stacked_gru <- function(num_tokens, embedding_dims, length_seq,
                          num_units_1, num_units_2)
  {
    model <- keras::keras_model_sequential() \%>\%
      keras::layer_embedding(input_dim = num_tokens,
                             output_dim = embedding_dims,
                             input_length = length_seq) \%>\%
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
