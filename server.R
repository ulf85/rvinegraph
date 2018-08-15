
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(networkD3)
library(VineCopula)
library(xtable)
#library(ggplot2)

shinyServer(function(input, output) {
  
  ## function to load the RVM object
  f.loadData <- reactive({
    if (!is.null(input$file)){
      myData = new.env()
      load(input$file$datapath, envir=myData)
      return(myData)
    } else {
      return(NULL)
    }
  })
  
  output$RVineTree <- renderSimpleNetwork({
    plotInput()
  })
  
  plotInput <- reactive({
    
    out2 <- f.loadData()
    
    if(!is.null(out2) & input$RVMname!=""){
      
      ## TODO: Name of the RVM-object may be different
      out <- eval(parse(text=paste0("out2$",input$RVMname)))
      
      if(!is.null(out) & !is.null(input$tree) & !is.na(input$tree)){
      
        Matrix <- out$Matrix
        family <- out$family
        par <- out$par
        par2 <- out$par2
        names <- out$names
        d <- dim(Matrix)[1]
        if(is.null(names)) names <- paste0("V", 1:d)
      
#     Matrix <- c(5, 2, 3, 1, 4,
#                 0, 2, 3, 4, 1,
#                 0, 0, 3, 4, 1,
#                 0, 0, 0, 4, 1,
#                 0, 0, 0, 0, 1)
#     Matrix <- matrix(Matrix, 5, 5)
#     
#     # define R-vine pair-copula family matrix
#     family <- c(0, 1, 3, 4, 4,
#                 0, 0, 3, 4, 1,
#                 0, 0, 0, 4, 1,
#                 0, 0, 0, 0, 3,
#                 0, 0, 0, 0, 0)
#     family <- matrix(family, 5, 5)
#     
#     # define R-vine pair-copula parameter matrix
#     par <- c(0, 0.2, 0.9, 1.5, 3.9,
#              0, 0, 1.1, 1.6, 0.9,
#              0, 0, 0, 1.9, 0.5,
#              0, 0, 0, 0, 4.8,
#              0, 0, 0, 0, 0)
#     par <- matrix(par, 5, 5)
#     
#     # define second R-vine pair-copula parameter matrix
#     par2 <- matrix(0, 5, 5)
#     
      # define RVineMatrix object
      RVM <- RVineMatrix(Matrix = Matrix, family = family,
                         par = par, par2 = par2,
                         names = names)
    
    
    M <- RVM$Matrix
    edge.labels <- c("family")
    legend <- FALSE
    empTauMat <- matrix(NA, d, d)
    
    edges <- list()
    for (j in 1:(d - 1)) edges[[j]] <- array(NA, dim = c(d - j, 2, j))
    
    weight <- list()
    for (j in 1:(d - 1)) weight[[j]] <- rep(NA, d - j)
    
    if (edge.labels[1] != FALSE) {
      numlabels <- length(edge.labels)
      elabels <- list()
      for (j in 1:(d - 1)) elabels[[j]] <- matrix(NA, d - j, numlabels)
    }
    
    # initial edge
    edges[[1]][1, , ] <- sort(c(M[d - 1, d - 1], M[d, d - 1]))
    weight[[1]][1] <- ifelse(is.null(data), theoTauMat[d, d - 1], empTauMat[d, d - 1])
    if (edge.labels[1] != FALSE) {
      for (jj in 1:numlabels) {
        if (edge.labels[jj] == "family") 
          elabels[[1]][1, jj] <- BiCopName(RVM$family[d, d - 1],
                                           short = TRUE)
        if (edge.labels[jj] == "par") 
          elabels[[1]][1, jj] <- parMat[d, d - 1]
        if (edge.labels[jj] == "par2") 
          elabels[[1]][1, jj] <- parMat2[d, d - 1]
        if (edge.labels[jj] == "theotau") 
          elabels[[1]][1, jj] <- theoTauMat[d, d - 1]
        if (edge.labels[jj] == "emptau") 
          elabels[[1]][1, jj] <- empTauMat[d, d - 1]
        if (edge.labels[jj] == "pair") 
          if (legend == TRUE) {
            elabels[[1]][1, jj] <- paste(RVM$Matrix[d - 1, d - 1],
                                         RVM$Matrix[d, d - 1],
                                         sep = ",")
          } else {
            elabels[[1]][1, jj] <- paste(RVM$names[RVM$Matrix[d - 1, d - 1]],
                                         RVM$names[RVM$Matrix[d, d - 1]],
                                         sep = ",")
          }
      }
    }
    
    for (i in (d - 2):1) {
      
      # new edge in first tree
      ee <- sort(c(M[i, i], M[d, i]))
      edges[[1]][d - i, , ] <- ee
      weight[[1]][d - i] <- ifelse(is.null(data), theoTauMat[d, i], empTauMat[d, i])
      if (edge.labels[1] != FALSE) {
        for (jj in 1:numlabels) {
          if (edge.labels[jj] == "family") 
            elabels[[1]][d - i, jj] <- BiCopName(RVM$family[d, i], 
                                                 short = TRUE)
          if (edge.labels[jj] == "par") 
            elabels[[1]][d - i, jj] <- parMat[d, i]
          if (edge.labels[jj] == "par2") 
            elabels[[1]][d - i, jj] <- parMat2[d, i]
          if (edge.labels[jj] == "theotau") 
            elabels[[1]][d - i, jj] <- theoTauMat[d, i]
          if (edge.labels[jj] == "emptau") 
            elabels[[1]][d - i, jj] <- empTauMat[d, i]
          if (edge.labels[jj] == "pair") 
            if (legend == TRUE) {
              elabels[[1]][d - i, jj] <- paste(RVM$Matrix[i, i],
                                               RVM$Matrix[d, i],
                                               sep = ",")
            } else {
              elabels[[1]][d - i, jj] <- paste(RVM$names[RVM$Matrix[i, i]],
                                               RVM$names[RVM$Matrix[d, i]],
                                               sep = ",")
            }
        }
      }
      # edges in further trees
      for (k in 1:(d - i - 1)) {
        edges[[k + 1]][d - i - k, 1, ] <- ee
        
        # identify conditioned and conditioning sets
        if (length(M[(d - k):d, i]) >= 3) {
          if (setequal(M[(d - k):d, i], ee_old)) {
            edges[[k + 1]][d - i - k, 2, ] <- ee_old
          } else {
            for (j in 1:(d - i - k)) {
              if (setequal(M[(d - k):d, i], edges[[k + 1]][j, 1, ])) 
                edges[[k + 1]][d - i - k, 2, ] <- edges[[k + 1]][j, 1, ]
              if (setequal(M[(d - k):d, i], edges[[k + 1]][j, 2, ])) 
                edges[[k + 1]][d - i - k, 2, ] <- edges[[k + 1]][j, 2, ]
            }
          }
        } else {
          edges[[k + 1]][d - i - k, 2, ] <- sort(M[(d - k):d, i])
        }
        
        # create edge lables
        weight[[k + 1]][d - i - k] <- ifelse(is.null(data), theoTauMat[d - k, i], empTauMat[d - k, i])
        if (edge.labels[1] != FALSE) {
          for (jj in 1:numlabels) {
            if (edge.labels[jj] == "family") 
              elabels[[k + 1]][d - i - k, jj] <- BiCopName(RVM$family[d - k, i], short = TRUE)
            if (edge.labels[jj] == "par") 
              elabels[[k + 1]][d - i - k, jj] <- parMat[d - k, i]
            if (edge.labels[jj] == "par2") 
              elabels[[k + 1]][d - i - k, jj] <- parMat2[d - k, i]
            if (edge.labels[jj] == "theotau") 
              elabels[[k + 1]][d - i - k, jj] <- theoTauMat[d - k, i]
            if (edge.labels[jj] == "emptau") 
              elabels[[k + 1]][d - i - k, jj] <- empTauMat[d - k, i]
            if (edge.labels[jj] == "pair") {
              if (legend == TRUE) {
                handle1 <- paste(RVM$Matrix[i, i], 
                                 RVM$Matrix[d - k, i],
                                 sep = ",")
                handle2 <- paste(RVM$Matrix[(d - k + 1):d, i],
                                 collapse = ",")
                handle3 <- paste(handle1, 
                                 handle2, 
                                 sep = ";")
              } else {
                handle1 <- paste(RVM$names[RVM$Matrix[i, i]], 
                                 RVM$names[RVM$Matrix[d - k, i]],
                                 sep = ",")
                handle2 <- paste(RVM$names[RVM$Matrix[(d - k + 1):d, i]],
                                 collapse = ",")
                handle3 <- paste(handle1, 
                                 handle2, 
                                 sep = ";")
              }
              elabels[[k + 1]][d - i - k, jj] <- handle3  #paste(handle1,handle2,sep=';')
            }
          }
        }
        
        # identify conditioned and conditioning sets
        ee <- c(sort(c(setdiff(ee, M[(d - k):d, i]), 
                       setdiff(M[(d - k):d, i], ee))),
                sort(intersect(ee, M[(d - k):d, i])))
      }
      
      ee_old <- ee
      
    }
    
    # label the nodes
    if (legend == FALSE) {
      for (j in 1:(d - 1)) for (k in 1:d) edges[[j]][edges[[j]] == k] <- RVM$names[k]
    }
    
    # convert to edge lists
    edgelist <- list()
    for (j in 1:(d - 1)) edgelist[[j]] <- matrix(NA, d - j, 2)
    
    edgelist[[1]] <- matrix(as.character(edges[[1]][, , 1]), d - 1, 2)
    
    for (j in 1:(d - 2)) edgelist[[2]][j, ] <- c(paste(edges[[2]][j, 1, ], collapse = ","), 
                                                 paste(edges[[2]][j, 2, ], collapse = ","))
    
    # separate conditioned and conditioning sets
    if (d > 3) {
      for (i in 3:(d - 1)) {
        for (j in 1:(d - i)) {
          edgelist[[i]][j, 1] <- paste(paste(edges[[i]][j, 1, 1:2], collapse = ","),
                                       paste(edges[[i]][j, 1, 3:i], collapse = ","), sep = ";")
          edgelist[[i]][j, 2] <- paste(paste(edges[[i]][j, 2, 1:2], collapse = ","),
                                       paste(edges[[i]][j, 2, 3:i], collapse = ","), sep = ";")
        }
      }
    }
    
      # combine edge lables
      if (edge.labels[1] != FALSE) {
        elabels2 <- list()
        for (j in 1:(d - 1)) {
          elabels2[[j]] <- rep(NA, d - j)
          for (i in 1:(d - j)) elabels2[[j]][i] <- paste(elabels[[j]][i, ], collapse = ",")
        }
      }
    
    
      # Create data.frame for plot
      tree <- 1  # for test here the first tree
      src <- c(edgelist[[input$tree]][,1])
      target <- c(edgelist[[input$tree]][,2])

      networkData <- data.frame(src, target)
    
      # Plot
      simpleNetwork(Data=networkData, fontSize = 15)
      }
    }
    
  })
  
  
  # for test
  output$table <- renderTable({
    out2 <- f.loadData()
    #if(!is.null(out)) return(out$RVM$family)
    if(!is.null(out2) ){
      #out <- eval(parse(text=paste0("out2$",input$RVMname)))
      content <- ls(out2)
      if(length(content)%%5 !=0){
        mis <- 5-length(content)%%5
        content <- c(content, rep("",mis))
      }
      out <- data.frame(matrix(content, ncol=5, byrow = TRUE))
      colnames(out)<-NULL
      return(out)
    }
    #data.frame(text)
    #data.frame(a=c("bla", "bla"))
  }, include.rownames = FALSE, 
  caption="List of objects in the loaded .RData file:",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL)
  )
  
  output$text <- renderUI({
    out2 <- f.loadData()
    if(!is.null(out2) ){
      out <- "List of objects in the loaded .RData file:<br>"
      content <- ls(out2)
      for(i in 1:length(content)){
        out <- paste(out, content[i])
        if(i%%5==0) out <- paste(out, "<br>")
      }
      HTML(out)
    }
  })
  
  
  
  output$downloadPlot <- downloadHandler(
    filename = "Shinyplot.pdf",
    content = function(file) {
      file.copy("Shinyplot.pdf", file)
    }
#       filename = function() {
#         paste("test.png",sep="")
#       },
#       content = function(file) {
#         p1 <- plotInput()
#         p1$save(file, standalone = TRUE)
#       }
    )  

})
