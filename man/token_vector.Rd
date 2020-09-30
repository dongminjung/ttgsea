\name{token_vector}
\alias{token_vector}

\title{
Vectorization of tokens
}
\description{
A vectorization of words or tokens of text is necessary for machine learning. Vectorized sequences are padded or truncated.
}
\usage{
token_vector(text, token, length_seq)
}


\arguments{
  \item{text}{
text data
}
  \item{token}{
result of tokenization (output of "text_token")
}
  \item{length_seq}{
length of input sequences
}
}


\value{
sequences of integers
}


\author{
Dongmin Jung
}


\seealso{
tm::removeWords, stopwords::stopwords, textstem::lemmatize_strings, tokenizers::tokenize_ngrams, keras::pad_sequences
}

\examples{
library(fgsea)
data(examplePathways)
data(exampleRanks)
names(examplePathways) <- gsub("_", " ",
                          substr(names(examplePathways), 9, 1000))
set.seed(1)
fgseaRes <- fgsea(examplePathways, exampleRanks)
tokens <- text_token(data.frame(fgseaRes)[,"pathway"],
          num_tokens = 1000)
sequences <- token_vector("Cell Cycle", tokens, 10)
}


