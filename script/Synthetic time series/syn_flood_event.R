syn1 <- structure(list(t = 1:20, 
                       discharge = c(0, 0.445364188838454, 
                                     0.67648582965612, 
                                     0.941847713557885, 
                                     1.2928102051699, 
                                     1.8378818077613,
                                     2.38849798386105, 
                                     2.99721341646724, 
                                     7.25822144471062, 
                                     25.276198249854, 
                                     32.4590403546072, 
                                     33.0677557872134, 
                                     11.032257126869, 
                                     5.67556131993451, 
                                     3.86082843647728, 
                                     3.24069958950972, 
                                     2.48265865234232, 
                                     1.92435246649882, 
                                     1.44594018118489, 
                                     0.927580945606176)), 
                  class = "data.frame", 
                  row.names = c(NA, -20L))

syn2 <- data.frame(t = 1:40,
                   discharge = c(syn1$discharge, syn1$discharge))
