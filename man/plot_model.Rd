\name{plot_model}
\alias{plot_model}

\title{
visualization of the model architecture
}
\description{
You are allowed to create a visualization of your model architecture.
}
\usage{
plot_model(x)
}

\arguments{
  \item{x}{
deep learning model
}
}


\value{
plot for the model architecture
}


\author{
Dongmin Jung
}


\seealso{
purrr::map, purrr::map_chr, purrr::pluck, purrr::imap_dfr, DiagrammeR::grViz
}

\examples{
inputs1 <- layer_input(shape = c(1000))
inputs2 <- layer_input(shape = c(1000))

predictions1 <- inputs1 \%>\%
  layer_dense(units = 128, activation = 'relu') \%>\%
  layer_dense(units = 64, activation = 'relu') \%>\%
  layer_dense(units = 32, activation = 'softmax')

predictions2 <- inputs2 \%>\%
  layer_dense(units = 128, activation = 'relu') \%>\%
  layer_dense(units = 64, activation = 'relu') \%>\%
  layer_dense(units = 32, activation = 'softmax')

combined <- layer_concatenate(c(predictions1, predictions2)) \%>\%
  layer_dense(units = 16, activation = 'softmax')


model <- keras_model(inputs = c(inputs1, inputs2),
                     outputs = combined)
model \%>\%
  compile(optimizer = "rmsprop",
          loss = "categorical_crossentropy",
          metrics = c("accuracy"))
plot_model(model)
}

