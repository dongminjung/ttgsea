text_token <- function(text, ngram_min = 1, ngram_max = 1, num_tokens) {
    text <- tm::removeWords(text, stopwords::stopwords("en"))
    text <- textstem::lemmatize_strings(text)
    it <- text2vec::itoken(text, preprocess_function = identity,
                           tokenizer = text2vec::word_tokenizer,
                           progressbar = FALSE)
    token <- text2vec::create_vocabulary(it, ngram = c(ngram_min, ngram_max),
                                         sep_ngram = " ")
    token <- text2vec::prune_vocabulary(token, vocab_term_max = num_tokens)

    result <- list()
    result$token <- token
    result$ngram_min <- ngram_min
    result$ngram_max <- ngram_max
    result
}




token_vector <- function(text, token, length_seq) {
    text <- tm::removeWords(text, stopwords::stopwords("en"))
    text <- textstem::lemmatize_strings(text)
    text_ngrams <- tokenizers::tokenize_ngrams(text, n = token$ngram_max,
                                               n_min = token$ngram_min,
                                               lowercase = FALSE,
                                               stopwords = stopwords::stopwords("en"))
    text_seqs <- lapply(text_ngrams, function(x) na.omit(match(x, token$token$term))-1)
    text_seqs %>% keras::pad_sequences(maxlen = length_seq)
}




metric_pearson_correlation <- function(y_true, y_pred) {
    y_true_dev <- y_true - keras::k_mean(y_true)
    y_pred_dev <- y_pred - keras::k_mean(y_pred)
    r_num <- keras::k_sum(y_true_dev * y_pred_dev)
    r_den <- keras::k_sqrt(keras::k_sum(keras::k_square(y_true_dev)) * 
                             keras::k_sum(keras::k_square(y_pred_dev)))
    r_num / r_den 
}




bi_lstm <- function(num_tokens, embedding_dims, length_seq, num_units) {
    model <- keras::keras_model_sequential() %>% 
      keras::layer_embedding(input_dim = num_tokens,
                             output_dim = embedding_dims,
                             input_length = length_seq) %>% 
      keras::bidirectional(keras::layer_lstm(units = num_units,
                                             activation = "relu")) %>%
      keras::layer_dense(1)
    model %>%
      keras::compile(loss = "mean_squared_error", optimizer = "adam",
                     metrics = custom_metric("pearson_correlation",
                                             metric_pearson_correlation))
}




bi_gru <- function(num_tokens, embedding_dims, length_seq, num_units) {
    model <- keras::keras_model_sequential() %>%
      keras::layer_embedding(input_dim = num_tokens,
                             output_dim = embedding_dims,
                             input_length = length_seq) %>%
      keras::bidirectional(keras::layer_gru(units = num_units,
                                            activation = "relu")) %>%
      keras::layer_dense(1)
    model %>%
      keras::compile(loss = "mean_squared_error", optimizer = "adam",
                     metrics = custom_metric("pearson_correlation",
                                             metric_pearson_correlation))
}




sampling_generator <- function(X_data, Y_data, batch_size) {
    function() {
        rows <- sample(seq_len(nrow(X_data)), batch_size, replace = TRUE)
        list(X_data[rows,], Y_data[rows,])
    }
}




fit_model <- function(gseaRes, text, score, model, ngram_min = 1, ngram_max = 2,
                      num_tokens, length_seq, epochs, batch_size,
                      use_generator = TRUE, ...) {
    gseaRes <- data.frame(gseaRes)
    
    message("pre-processing...")
    # text tokenization
    tokens <- text_token(gseaRes[,text], ngram_min, ngram_max, num_tokens)
    x_train <- token_vector(gseaRes[,text], tokens, length_seq)
    y_train <- gseaRes[,score]
    
    message("model fitting...")
    if(use_generator) {
      model %>% keras::fit_generator(sampling_generator(as.matrix(x_train), 
                                                        as.matrix(y_train),
                                                        batch_size = batch_size),
                                     steps_per_epoch = nrow(x_train)/batch_size, 
                                     epochs, ...)
    } else {
      model %>% keras::fit(x_train, y_train, batch_size, epochs, ...)
    }
    
    # prediction for every token
    token_term <- textstem::lemmatize_strings(tokens$token$term)
    token_term_vector <- token_vector(token_term, tokens, length_seq)
    pred <- as.vector(predict(model, token_term_vector))
    token_pred <- data.frame(token_term, pred)
    token_pred <- token_pred[order(token_pred$pred, decreasing = TRUE),]
    
    message("post-processing...")
    # text for each token
    gsea_text <- gsub("[[:punct:]]", " ",
                      tm::removeWords(gseaRes[,text],
                                      stopwords::stopwords("en")))
    gsea_text <- textstem::lemmatize_strings(gsea_text)
    token_gsea <- lapply(token_term,
                         function(x) gseaRes[which(grepl(x, gsea_text)),])
    names(token_gsea) <- token_term
    
    result <- list()
    result$model <- model
    result$token_pred <- data.table::data.table(token_pred)
    result$token_gsea <- token_gsea
    result$tokens <- tokens
    result$num_tokens <- num_tokens
    result$length_seq <- length_seq
    result
}




predict_model <- function(object, new_text, num_simulations = 1000,
                          adj_p_method = "fdr") {
    model <- object$model
    num_tokens <- object$num_tokens
    length_seq <- object$length_seq
    token <- object$tokens
    new_text <- textstem::lemmatize_strings(new_text)
    x_test <- token_vector(new_text, token, length_seq)
    test_value <- as.vector(predict(model, x_test))
    
    m <- num_simulations
    n <- rowSums(x_test != 0)
    MC_p_value <- ifelse(n, mapply(function(x, y) {
      simulation_matrix <- matrix(sample(0:(num_tokens-1),
                                         m*x, replace = TRUE), m, x)
      temp <- lapply(seq_len(nrow(simulation_matrix)),
                     function(z) simulation_matrix[z,])
      x_test_temp <- temp %>% keras::pad_sequences(maxlen = length_seq)
      2*min(mean(as.vector(predict(model, x_test_temp)) > abs(y)),
            mean(as.vector(predict(model, x_test_temp)) < -abs(y)))
    }, n, test_value), NA)
    
    adj_p_value <- p.adjust(MC_p_value, adj_p_method, length(MC_p_value))
    
    data.frame(new_text, test_value, MC_p_value, adj_p_value)
}




plot_model <- function(x) {
    # layer information
    model_layers <- x$get_config()$layers
    if(length(x$layers) == (length(model_layers)-1)) {
      model_layers <- list()
      for(i in 1:length(x$layers)) {
        model_layers[[i]] <- x$get_config()$layers[[i+1]]
      }
    }
    
      # node information
    layer_name <- model_layers %>% 
      purrr::map_chr(~(purrr::`%||%`(purrr::pluck(., "config", "name"), "")))
    layer_name_sub <- model_layers %>% 
      purrr::map_chr(~(purrr::`%||%`(purrr::pluck(., "config", "layer", "config", "name"), "")))
    layer_type <- model_layers %>% 
      purrr::map_chr("class_name")
    layer_type_sub <- model_layers %>% 
      purrr::map_chr(~(purrr::`%||%`(purrr::pluck(., "config", "layer", "class_name"), "")))
    layer_input_shape <- unlist(lapply(x$layers, function(x) {
        ifelse(length(purrr::pluck(x$input_shape, 1)) > 0,
               paste("[", paste(unlist(lapply(x$input_shape,
                   function(x) paste("(", toString(paste(x)), ")", sep = ""))),
                   collapse = ", "), "]", sep = ""),
               paste("(", toString(paste(x$input_shape)), ")", sep = ""))
        }))
    layer_input_shape <- gsub("NULL", "None", layer_input_shape)
    layer_output_shape <- unlist(lapply(x$layers, function(x) {
        ifelse(length(purrr::pluck(x$output_shape, 1)) > 0,
               paste("[", paste(unlist(lapply(x$output_shape,
                   function(x) paste("(", toString(paste(x)), ")", sep = ""))),
                   collapse = ", "), "]", sep = ""),
               paste("(", toString(paste(x$output_shape)), ")", sep = ""))
        }))
    layer_output_shape <- gsub("NULL", "None", layer_output_shape)
    node_info <- data.frame(layer_name, layer_name_sub, layer_type,
                            layer_type_sub, layer_input_shape, layer_output_shape)
    
    # edge information
    inbound <- lapply(model_layers,
                      function(x) switch(length(x$inbound_nodes),
                                         x$inbound_nodes[[1]] %>% 
                                           purrr::map_chr(c(1, 1))))
    
    if (length(Filter(Negate(is.null), inbound)) == 0) {
      edge_info <- embed(rownames(node_info), dimension = 2)
      from <- edge_info[, 2]
      to <- edge_info[, 1]
      edge_info <- data.frame(from, to, stringsAsFactors = FALSE)
    } else {
      names(inbound) <- purrr::map(model_layers, "name")
      inbound <- Filter(Negate(is.null), inbound)
      edge_info <- purrr::imap_dfr(inbound, ~ data.frame(from = .x, to = .y, stringsAsFactors = FALSE))
      edge_info$from <- rownames(node_info)[match(edge_info$from, node_info$layer_name)]
      edge_info$to <- rownames(node_info)[match(edge_info$to, node_info$layer_name)]
    }
  
    # plot
    nodes <- paste(unlist(lapply(seq_len(nrow(node_info)),
                          function(x) paste(toString(x), " [label = '@@", toString(x), "']", sep = ""))),
                   collapse = "")
    names <- paste(unlist(lapply(seq_len(nrow(node_info)), function(x) {
        paste(" [", toString(x), "]: ", "'", node_info$layer_name[x], " ",
              ifelse(node_info$layer_name_sub[x] != "",
                     paste("(", node_info$layer_name_sub[x], ")", sep = ""), ""),
              " : ",
              node_info$layer_type[x], " ",
              ifelse(node_info$layer_type_sub[x] != "",
                     paste("(", node_info$layer_type_sub[x], ")", sep = ""), ""),
              "|{input: | output:}", "|{",
              node_info$layer_input_shape[x], "|", node_info$layer_output_shape[x], "}", "'", sep = "")
        })), collapse = "\n")
    edges <- gsub(",", "->", paste(apply(edge_info, 1, toString), collapse = " "))
    
    DiagrammeR::grViz(paste("digraph{
                             graph [layout = dot]
                             node [shape = record]",
                             nodes, edges, "} \n", names))
}
