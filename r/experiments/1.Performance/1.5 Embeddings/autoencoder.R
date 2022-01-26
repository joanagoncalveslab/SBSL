library(keras)

autoencode <- function(X.train, X.test, hidden = NULL) {
  input_size <- dim(X.train)[2]
  latent_size <- hidden
  
  encoder_input <- layer_input(shape = input_size)
  encoder_output <- 
    encoder_input %>%
    layer_dense(units = latent_size, activation = "tanh")
  encoder <- keras_model(encoder_input, encoder_output)
  
  decoder_input <- layer_input(shape = latent_size)
  decoder_output <- 
    decoder_input %>%
    layer_dense(units = input_size)
  decoder <- keras_model(decoder_input, decoder_output)
  
  aen_input = layer_input(shape = input_size)
  aen_output = aen_input %>% 
    encoder() %>% 
    decoder()
  
  autoencoder_model = keras_model(aen_input, aen_output)
  
  autoencoder_model %>% compile(
    loss='mean_squared_error',
    optimizer='adam',
    metrics = c('accuracy')
  )
  
  summary(autoencoder_model)
  
  history <-
    autoencoder_model %>%
    keras::fit(X.train,
               X.train,
               epochs=100,
               shuffle=TRUE,
               validation_data= list(X.test, X.test),
               callbacks = list(
                 callback_early_stopping(monitor = "val_loss")
               )
    )
  
  plot(history)
  
  encoder
}

binary_autoencoder <- function(X.train, X.test, hidden = NULL) {
  input_size <- dim(X.train)[2]
  latent_size <- hidden
  
  encoder_input <- layer_input(shape = input_size)
  encoder_output <- 
    encoder_input %>%
    layer_dense(units = latent_size, activation = "tanh")
  encoder <- keras_model(encoder_input, encoder_output)
  
  decoder_input <- layer_input(shape = latent_size)
  decoder_output <- 
    decoder_input %>%
    layer_dense(units = input_size, activation = "sigmoid")
  decoder <- keras_model(decoder_input, decoder_output)
  
  aen_input = layer_input(shape = input_size)
  aen_output = aen_input %>% 
    encoder() %>% 
    decoder()
  
  autoencoder_model = keras_model(aen_input, aen_output)
  
  autoencoder_model %>% compile(
    loss='binary_crossentropy',
    optimizer='adam'
  )
  
  summary(autoencoder_model)
  
  history <-
    autoencoder_model %>%
    keras::fit(X.train,
               X.train,
               epochs=100,
               shuffle=TRUE,
               validation_data= list(X.test, X.test),
               callbacks = list(
                 callback_early_stopping(monitor = "val_loss")
               )
    )
  
  plot(history)
  
  encoder
}


