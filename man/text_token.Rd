\name{text_token}
\alias{text_token}

\title{
Tokenizing text
}
\description{
An n-gram is used for tokenization. This function can also be used to limit the total number of tokens.
}
\usage{
text_token(text, ngram_min = 1, ngram_max = 1, num_tokens)
}

\arguments{
  \item{text}{
text data
}
  \item{ngram_min}{
minimum size of an n-gram (default: 1)
}
  \item{ngram_max}{
maximum size of an n-gram (default: 1)
}
  \item{num_tokens}{
total number of tokens
}
}


\value{
  \item{token}{
result of tokenizing text
}
  \item{ngram_min}{
minimum size of an n-gram
}
  \item{ngram_max}{
maximum size of an n-gram
}
}


\author{
Dongmin Jung
}


\seealso{
tm::removeWords, stopwords::stopwords, textstem::lemmatize_strings, text2vec::create_vocabulary, text2vec::prune_vocabulary
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
}

