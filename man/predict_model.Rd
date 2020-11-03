\name{predict_model}
\alias{predict_model}

\title{
Model prediction
}
\description{
From the result of the function "ttgsea", we can predict enrichment scores. For each new term, lemmatized text, predicted enrichment score, Monte Carlo p-value and adjusted p-value are provided. The function "token_vector" is used for encoding as we did for training. Of course, mapping from tokens to integers should be the same.
}
\usage{
predict_model(object, new_text, num_simulations = 1000,
              adj_p_method = "fdr")
}

\arguments{
  \item{object}{
result of "ttgsea"
}
  \item{new_text}{
new text data
}
  \item{num_simulations}{
number of simulations for Monte Carlo p-value (default: 1000)
}
  \item{adj_p_method}{
correction method (default: "fdr")
}
}


\value{
table for lemmatized text, predicted enrichment score, MC p-value and adjusted p-value
}


\author{
Dongmin Jung
}


\seealso{
stats::p.adjust
}

\examples{
library(reticulate)
if (keras::is_keras_available() & reticulate::py_available()) {
  library(fgsea)
  data(examplePathways)
  data(exampleRanks)
  names(examplePathways) <- gsub("_", " ",
                            substr(names(examplePathways), 9, 1000))
  set.seed(1)
  fgseaRes <- fgsea(examplePathways, exampleRanks)
  
  num_tokens <- 1000
  length_seq <- 30
  batch_size <- 32
  embedding_dims <- 50
  num_units <- 32
  epochs <- 1
  
  ttgseaRes <- fit_model(fgseaRes, "pathway", "NES",
                         model = bi_gru(num_tokens,
                                        embedding_dims,
                                        length_seq,
                                        num_units),
                         num_tokens = num_tokens,
                         length_seq = length_seq,
                         epochs = epochs,
                         batch_size = batch_size)
  
  set.seed(1)
  predict_model(ttgseaRes, "Cell Cycle")
}
}
