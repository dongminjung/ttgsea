\name{fit_model}
\alias{fit_model}

\title{
Deep learning model fitting
}
\description{
From the result of GSEA, we can predict enrichment scores for unique tokens or words from text in names of gene sets by using deep learning. The function "text_token" is used for tokenizing text and the function "token_vector" is used for encoding. Then the encoded sequence is fed to the embedding layer of the model.
}
\usage{
fit_model(gseaRes, text, score, model, ngram_min = 1, ngram_max = 2,
          num_tokens, length_seq, epochs, batch_size, ...)
}

\arguments{
  \item{gseaRes}{
a table with GSEA result having rows for gene sets and columns for text and scores  
}
  \item{text}{
column name for text data
}
  \item{score}{
column name for enrichment score
}
  \item{model}{
deep learning model, input dimension and length for the embedding layer must be same to the "num_token" and "length_seq", respectively
}
  \item{ngram_min}{
minimum size of an n-gram (default: 1)
}
  \item{ngram_max}{
maximum size of an n-gram (default: 2)
}
  \item{num_tokens}{
total number of tokens, it must be equal to the input dimension of "layer_embedding" in the "model"
}
  \item{length_seq}{
length of input sequences, it must be equal to the input length of "layer_embedding" in the "model"
}
  \item{epochs}{
number of epochs
}
  \item{batch_size}{
batch size
}
  \item{...}{
additional parameters for the "fit_generator"
}
}


\value{
  \item{model}{
trained model
}
  \item{tokens}{
information for tokens
}
  \item{token_pred}{
prediction for every token, each row has a token and its predicted score
}
  \item{token_gsea}{
list of the GSEA result only for the corresponding token
}
  \item{num_tokens}{
total number of tokens
}
  \item{length_seq}{
length of input sequences
}
}


\author{
Dongmin Jung
}


\seealso{
keras::fit_generator, keras::layer_embedding, keras::pad_sequences, textstem::lemmatize_strings, text2vec::create_vocabulary, text2vec::prune_vocabulary
}

\examples{
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
}

