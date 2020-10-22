library(fgsea)
data(examplePathways)
data(exampleRanks)
names(examplePathways) <- gsub("_", " ", substr(names(examplePathways), 9, 1000))
set.seed(1)
fgseaRes <- fgseaSimple(examplePathways, exampleRanks, nperm = 1000)

num_tokens <- 100
length_seq <- 30
batch_size <- 32
embedding_dims <- 50
num_units <- 32
epochs <- 1


library(reticulate)
if (keras::is_keras_available() & reticulate::py_available()) {
  ttgseaRes <- fit_model(fgseaRes, "pathway", "NES",
                         model = bi_lstm(num_tokens, embedding_dims,
                                         length_seq, num_units),
                         num_tokens = num_tokens,
                         length_seq = length_seq,
                         epochs = epochs,
                         batch_size = batch_size)
}



check_python_available <- function() {
  if (!reticulate::py_available()) {
    skip("Python is not available on this system")
  }
}


check_keras_available <- function() {
  if (!keras::is_keras_available()) {
    skip("Keras is not available on this system")
  }
}



test_that("fit_model: the output of fit_model is list", {
  check_python_available()
  check_keras_available()
  expect_type(ttgseaRes, "list")
})


test_that("fit_model: fit_model yields a model", {
  check_python_available()
  check_keras_available()
  expect_type(ttgseaRes$model, "closure")
})


test_that("fit_model: num_tokens in the model", {
  check_python_available()
  check_keras_available()
  expect_equal(ttgseaRes$num_tokens, num_tokens)
})


test_that("fit_model: length_seq in the model", {
  check_python_available()
  check_keras_available()
  expect_equal(ttgseaRes$length_seq, length_seq)
})


test_that("fit_model: swap values of text and score", {
  check_python_available()
  check_keras_available()
  expect_error(fit_model(fgseaRes, "NES", "pathway",
                         model = bi_lstm(num_tokens, embedding_dims,
                                         length_seq, num_units),
                         num_tokens = num_tokens,
                         length_seq = length_seq,
                         epochs = epochs,
                         batch_size = batch_size))
})


test_that("fit_model: miss model", {
  check_python_available()
  check_keras_available()
  expect_error(fit_model(fgseaRes, "pathway", "NES",
                         num_tokens = num_tokens,
                         length_seq = length_seq,
                         epochs = epochs,
                         batch_size = batch_size))
})


test_that("fit_model: miss epochs", {
  check_python_available()
  check_keras_available()
  expect_error(fit_model(fgseaRes, "pathway", "NES",
                         model = bi_lstm(num_tokens, embedding_dims,
                                         length_seq, num_units),
                         num_tokens = num_tokens,
                         length_seq = length_seq,
                         batch_size = batch_size))
})


test_that("fit_model: miss batch_size", {
  check_python_available()
  check_keras_available()
  expect_error(fit_model(fgseaRes, "pathway", "NES",
                         model = bi_lstm(num_tokens, embedding_dims,
                                         length_seq, num_units),
                         num_tokens = num_tokens,
                         length_seq = length_seq,
                         epochs = epochs))
})



test_that("predict_model: the output of predict_model is list", {
  check_python_available()
  check_keras_available()
  expect_type(predict_model(ttgseaRes, "Cell Cycle"), "list")
})


test_that("predict_model: the number of rows of the result is the length of new_text", {
  check_python_available()
  check_keras_available()
  expect_equal(nrow(predict_model(ttgseaRes, c("Cell Cycle", "DNA Replication"))), 2)
})


test_that("predict_model: miss new_text", {
  check_python_available()
  check_keras_available()
  expect_error(predict_model(ttgseaRes))
})



test_that("plot_model: output", {
  check_python_available()
  check_keras_available()
  expect_visible(plot_model(ttgseaRes$model))
})


test_that("plot_model: only model is required", {
  check_python_available()
  check_keras_available()
  expect_error(plot_model(ttgseaRes))
})


