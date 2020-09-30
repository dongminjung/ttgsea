\name{bi_lstm}
\alias{bi_lstm}

\title{
Bidirectional LSTM
}
\description{
A predefined function that is used as a model in "ttgsea". This is a simple model, but you can define your own model.
}
\usage{
bi_lstm(num_tokens, embedding_dims, length_seq, num_units)
}

\arguments{
  \item{num_tokens}{
total number of tokens
}
  \item{embedding_dims}{
a non-negative integer for dimension of the dense embedding
}
  \item{length_seq}{
length of input sequences, input length of "layer_embedding""
}
  \item{num_units}{
dimensionality of the output space in the LSTM layer
}
}


\value{
model
}


\author{
Dongmin Jung
}


\seealso{
keras::keras_model_sequentia, keras::layer_embedding, keras::layer_lstm,
keras::bidirectional, keras::layer_dense, keras::compile
}


\examples{
num_tokens <- 1000
length_seq <- 30
embedding_dims <- 50
num_units <- 32
model <- bi_lstm(num_tokens, embedding_dims, length_seq, num_units)


# stacked lstm
num_units_1 <- 32
num_units_2 <- 16
stacked_lstm <- function(num_tokens, embedding_dims, length_seq,
                         num_units_1, num_units_2)
{
  model <- keras::keras_model_sequential() \%>\%
    keras::layer_embedding(input_dim = num_tokens,
                           output_dim = embedding_dims,
                           input_length = length_seq) \%>\%
    keras::layer_lstm(units = num_units_1,
                      activation = "relu",
                      return_sequences = TRUE) \%>\%
    keras::layer_lstm(units = num_units_2,
                      activation = "relu") \%>\%
    keras::layer_dense(1)
    
  model \%>\%
    keras::compile(loss = "mean_squared_error",
                   optimizer = "adam",
                   metrics = metric_pearson_correlation)
}
}

