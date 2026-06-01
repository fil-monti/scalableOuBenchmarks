# Project: AdditiveGaussianProcesses
# This files stores some functions to plot the output of the MCMC's
# ---- plots ----
library(R6)
plots <-
  R6Class("plots",
          public = list(
            name = NULL,
            initialize = function(name = NULL) {
              self$name <- name
            }
          )
  )

normalizeRates = function(logRates) {
  n_rates <- length(logRates)
  n_states <- 0.5 * (1 + sqrt(1 + 4 * n_rates))
  Q <- exp(logRates) |> getQ()
  normalizationConstant <- - sum(diag(Q)) / n_states
  # cat("The normalization constant is: ", normalizationConstant, "\n")
  logRates <- logRates - log(normalizationConstant)
  return(logRates)
}
              
# ---- plotLogRatesQuantiles ----
plotLogRatesQuantiles <- function(summStat, trueLogRates, predictor, jobName) {
  summaryLogRates <- summStat |> row.names() |> 
    sapply(function(x) grepl("logRates", x)) %>% summStat[.,] |> as_tibble()
  lmod <- lm(median ~ predictor, summaryLogRates)
  trueLogRates <- normalizeRates(exp(trueLogRates))
  summaryLogRatesPred <- cbind(trueLogRates, predictor, summaryLogRates) |> as_tibble()
  slope <- coef(lmod)[2] %>% round(2) %>% as.character()
  slope <- paste0("\nSlopes\nTrue: 1\n", "Est: ", slope)
  ytextSlope <- Inf
  ## max(max(summaryLogRates$q75), max(predictor - mean(predictor)))
  plotLogRates <- ggplot(summaryLogRatesPred, aes(x = predictor)) + 
    geom_errorbar(aes(ymin = q25, ymax = q75, color = "darkgreen"), width = 0.05) +
    # geom_errorbar(aes(ymin = q05, ymax = q95, color = "green"), width = 0.1) +
    # geom_point(aes(y = median, color = "blue"), size = 1.5) + 
    # geom_smooth(aes(y = median), method = "lm", se = FALSE, color = "blue", linewidth = 0.4) +
    geom_text(aes(x = min(predictor), y = ytextSlope), label = slope, hjust = -0.1, vjust = 1.5, size = 5, color = "blue") +
    # geom_point(aes(y = trueLogRates, color = "red"), size = 1.5) +
    #geom_smooth(aes(y = logRates), method = "lm", se = FALSE, color = "red", linewidth = 0.4) +
    #geom_line(aes(y = logRates, color = "red"), linewidth = 0.4) +
    labs(title = paste0(jobName), x = "Predictor", y = "Log rates") + 
    scale_color_manual(name = "Legend",  # Adding manual color scale with legend
                       values = c("blue" = "blue", "darkgreen" = "darkgreen", "red" = "red"),
                       labels = c("Median", "25%-75% quantiles", "True log rates")) +
    theme(legend.text = element_text(size = 20), 
          legend.title = element_text(size = 20), 
          legend.key.size = unit(2, 'cm')) +
    theme_minimal() 
  plotLogRates
  #geom_pointrange(aes(y = median, ymin = q25, ymax = q75), color = "green") +
  # annotate("text", x = -Inf, y = Inf, label = slop, size = 5, color = "blue") +
}

boxPlotLogRates <- function(data, trueLogRates, predictor, jobName) {
  nrows <- nrow(data)
  parameters <- data[, grepl("^logRates", names(data))]
  colnames(parameters) <- gsub("^logRates", "", colnames(parameters))
  parametersLong <- parameters |> 
    pivot_longer(cols = everything(), names_to = "estimates", values_to = "y")
  parametersLong$estimates <- factor(parametersLong$estimates,
                                     levels = unique(parametersLong$estimates))
  parametersLong$trueLogRates <- rep(trueLogRates, nrows)
  boxPlot <- ggplot(parametersLong, aes(x = estimates, y = y)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_point(aes(y = trueLogRates), color = "red") +
    #position = position_jitter(width = 0.2)) +
    labs(title = paste0(jobName), x = "") +
    scale_y_continuous(limits = quantile(parametersLong$y, c(0.01, 0.99)))
  boxPlot
}

boxPlotHyper <- function(data, jobName) {
  parameters <- data[, grepl("^logRates", names(data))]
  colnames(parameters) <- gsub("^logRates", "", colnames(parameters))
  parametersLong <- parameters |> 
    pivot_longer(cols = everything(), names_to = "estimates", values_to = "y")
  parametersLong$estimates <- factor(parametersLong$estimates,
                                     levels = unique(parametersLong$estimates))
  
  boxPlot <- ggplot(parametersLong, aes(x = estimates, y = y)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste0(jobName), x = "") +
    scale_y_continuous(limits = quantile(parametersLong$y, c(0.01, 0.99)))
  boxPlot
}

# ---- plotLogRatesRabies ----
plottingLogRates <- 
  R6Class("plottingLogRates",
          inherit = plots,
          public = list(
             summStat = NULL,
             otherData = NULL, 
             predictor = NULL,
             predictorName = NULL,
             jobName = NULL,
             data = list(),
             dataRates = list(),
             other = list(),
             attributes = list(),
             actions = list(),
             
             initialize = function(summStat, otherData, predictor, predictorName, jobName, attributes, actions = list()) {
               super$initialize(name = "plottingLogRates")
               self$summStat <- summStat
               self$predictor <- predictor
               self$predictorName <- predictorName
               self$jobName <- jobName
               self$attributes <- attributes
               self$actions <- actions
               
               # self$preProcessing()
               self$data <- summStat
               self$other <- otherData
             }, 
             
             
             preProcessing = function() {
               summaryLogRates <- self$summStat |> row.names() |> 
                 sapply(function(x) 
                   grepl("logRates", x)) %>% 
                 self$summStat[.,] |> as_tibble() 
               # summaryRates <- self$summStat |> row.names() |> 
               #   sapply(function(x) 
               #     grepl("rates", x)) %>% 
               #   self$summStat[.,] |> as_tibble()
               # |> arrange(median) # TODO delete arrange
               # summaryLogRates <- summaryLogRates %>% mutate(across(everything(), ~ . - 3))
               cat("Stats and predictor lengths are ", length(summaryLogRates$median), "and", length(self$predictor), "\n")
               lmod <- lm(median ~ self$predictor, summaryLogRates)
               # self$predictor <- sort(self$predictor)
               summaryLogRatesPred <- cbind(self$predictor, summaryLogRates) |> as_tibble()
               # summaryRatesPred <- cbind(self$predictor, summaryRates) |> as_tibble()
               slope <- coef(lmod)[2] %>% round(4) %>% as.character()
               slope <- paste0("\nSlopes\nTrue: 1\n", "Est: ", slope)
               ytextSlope <- Inf
               ymin <- min(summaryLogRates$q05)
               ymax <- max(summaryLogRates$q95)
               # trueLine <- self$predictor*0.23 - 2.80
               # trueLine <- -self$predictor*0.5 - 2.50
               if (self$attributes$predictorName == "Body size") {
                 trueLine <- -0.14581 * self$predictor - 2.8
               } else if (self$attributes$predictorName == "Host genetic distance") {
                 # trueLine <- -0.6171 * self$predictor - 3.1
                 trueLine <- normalizeRates(-0.5694 * self$predictor) + log(0.007693286)
                 # - 3.1
               } else if (self$attributes$predictorName == "Range overlap") {
                 trueLine <- normalizeRates(0.428 * self$predictor) 
                   # normalizeRates(0.2957 * self$predictor - 2.8)
               } else if (self$attributes$predictorName == "SmallHost") {
                 trueLine <- normalizeRates(self$predictor^2)
                   # as.numeric(strsplit(predictor, " ")[[1]]) "0 -0.6931472 -1.098612 -1.098612 -0.6931472 0 0 -0.6931472 -1.098612 -1.098612 -0.6931472 0 -0.6931472 -1.098612 -1.098612 0 -0.6931472 -1.098612 0 -0.6931472 0 0 -0.6931472 -1.098612 -1.098612 -0.6931472 0 0 -0.6931472 -1.098612 -1.098612 -0.6931472 0 -0.6931472 -1.098612 -1.098612 0 -0.6931472 -1.098612 0 -0.6931472 0"
               } else if (self$attributes$predictorName == "Total seats") {
                 cat("Choosing slope for total seats \n")
                 trueLine <-  
                   normalizeRates(1.3113*self$predictor) #1.3113*self$predictor - 3
               } else if (self$attributes$predictorName == "Origin population density") {
                 cat("Choosing slope for donor population density \n")
                 trueLine <-  
                   normalizeRates(0.30716*self$predictor) #0.30716*self$predictor - 3
               } else {
                 cat("Choosing default slope \n")
                 trueLine <-  normalizeRates(0.5*self$predictor)
               }
               
               
               if (!is.null(self$attributes$shape) && nzchar(self$attributes$shape)) {
                 cat("Choosing shape \n")
                 if (self$attributes$shape == "Concave") {
                   trueLine <- normalizeRates(- 0.3 * self$predictor^2 - 1.7 * self$predictor - 1)
                 } else if  (self$attributes$shape == "Convex") {
                   trueLine <- normalizeRates(0.2 * self$predictor^2 + 0.3 * self$predictor - 3.5) + log(0.2)
                 }
               }
               
               # cat("trueLine is ", trueLine, "\n")
               # cat("predictor is ", self$predictor, "\n")
               # cat("medians are ", summaryLogRatesPred$median, "\n")
               
               
               self$data <- summaryLogRatesPred
               # self$dataRates <- summaryRatesPred
               self$other <- list(slope = slope, ytextSlope = ytextSlope, 
                                  ymin = ymin, ymax = ymax, 
                                  trueLine = trueLine)
               # print(names(self$data))
             }, 
             
             publicationQualitySettings = function(myplot) {
               library(sysfonts)
               sysfonts::font_add("CMU Serif", regular = "/Users/filippomonti/Library/Fonts/cmunrm.otf")

               fontSize <- 12
               fontSizeSmall <- 10
               fontType <- "CMU Serif"
               
               if (self$attributes$predictorName == "Host genetic distance") {
                 myplot <- myplot +
                   theme( legend.position = c(0.95, 0.9),
                          legend.justification = c(1, 1))
               } else {
                 myplot <- myplot +
                   theme( legend.position = c(0, 1),
                          legend.justification = c(0, 1)) 
                 #  xlim(-0.3, 1.5) +
                 # yliim(-4, -1)
               }
               
               myplot <- myplot +
                 theme_classic() +
                 theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # frame around plot)
                   panel.grid.major.y = element_line(color = "grey90", linetype = "33", size = 0.5),
                   panel.grid.minor.y = element_blank(),  # optional: remove minor grid lines
                   panel.grid.major.x = element_blank(),  # optional: remove vertical lines if you want
                    axis.ticks.x = element_line(size = 0.3),
                    axis.ticks.y = element_line(size = 0.3),
                    axis.line = element_blank()
                 ) +
                 theme(
                   panel.grid.major.y = element_line(color = "grey90", linetype = "33", size = 0.5),
                   # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                   legend.box.just = "right",
                   legend.key.width = unit(0.5, "lines"),
                   legend.spacing.y = unit(0, "cm"),          # reduce space between legend items
                   legend.key.height = unit(0.3, "cm"),        # reduce height of each legend key (line)
                   # legend.margin = margin(0, 0, 0, 0),
                   # legend.margin = margin(6, 6, 6, 6),
                   legend.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),
                   # legend.title = element_text(size = fontSizeSmall, family = fontType, color = "black", face = "bold"),
                   legend.title=element_blank(),
                   legend.background = element_rect(
                     fill = alpha("white", 0))  # Semi-transparent white background
                   # color = "black",             # Border color
                   # size = 0.5)                   # Border thickness
                 ) +
                 theme(text = element_text(size = fontSize, family = fontType, color = "black"),  # Set base text size
                       axis.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),  # Axis tick labels
                       axis.title = element_text(size = fontSize, family = fontType, color = "black")
                 )
               return(myplot)
             },
             
             plottingDerivatives = function() {
               summaryLogRates <- self$summStat |> row.names() |> 
                 sapply(function(x) 
                   grepl("prediction", x)) %>% 
                 self$summStat[.,] |> as_tibble() 
               summaryLogRates |> names() |>  print()
               cat("creating Derivative plot \n")
               summaryLogRatesPred <- cbind(self$predictor, summaryLogRates) |> as_tibble()
               myplotDer <- ggplot(summaryLogRatesPred, aes(x = self$predictor)) + 
                 geom_ribbon(aes(ymin = q05, ymax = q95), fill="darkgoldenrod1", alpha = 0.5) +
                 geom_line(aes(y=median, , color = "b")) +
                 geom_line(aes(y = (0.4 * self$predictor + 0.3), color = "a"), size = 0.5) + 
                 labs(x = paste0(self$predictorName), y = "Log rates derivatives") + 
                 scale_color_manual(name = "",  # Adding manual color scale with legend
                                    values = c("a" = "blue", "b" = "#F8766D"),
                                    # "c" = "darkgoldenrod1"),
                                    labels = c( "Truth", "GP estimates")) 
               
               plotLogRates <- ggplot(self$data, aes(x = self$predictor)) + 
                 # geom_point(aes(y = median, color = "GP"), size = 2.5) +
                 # geom_point(aes(y = self$other$trueLine, color = "TrueValues"), size = 2.5) +
                 geom_ribbon(aes(ymin = q05, ymax = q95), fill="darkgoldenrod1", alpha = 0.5)  #darkgoldenrod1
               
               return(self$publicationQualitySettings(myplotDer))
             },
             
             plottingSimulation = function() {
               cat("Plotting simulation \n")
               self$data <- summ_data_GP %>%
                 mutate(
                   lowHPD = predict(gam(q05 ~ s(predictor, bs = "cs")), newdata = self$data),
                   highHPD = predict(gam(q95 ~ s(predictor, bs = "cs")), newdata = self$data)
                 )
               plotLogRates <- ggplot(self$data, aes(x = self$predictor)) + 
                 # geom_point(aes(y = median, color = "GP"), size = 2.5) +
                 # geom_point(aes(y = self$other$trueLine, color = "TrueValues"), size = 2.5) +
                 geom_ribbon(aes(ymin = lowHPD, ymax = highHPD), fill="green", alpha = 0.5) + #darkgoldenrod1
                 
                 # geom_line(aes(y = median, color = "b"), size = 0.5) + 
                 geom_smooth(aes(y = median, color = "b"), se = FALSE, size = 0.0, alpha=0.0)
                 geom_line(aes(y = self$other$trueLine, color = "a"), size = 0.5) + 
                 # geom_line(aes(y = q05, color = "c"), size = 1.5, alpha = 0) +
                 labs(x = paste0(self$predictorName), y = "Log rates") + 
                 theme_minimal() +
                 scale_color_manual(name = "",  # Adding manual color scale with legend
                                    values = c("a" = "blue", "b" = "#F8766D"),
                                    # "c" = "darkgoldenrod1"),
                                    labels = c( "Truth", "GP estimates")) 
               
               if (!self$actions$publicationQuality) {
                 plotLogRates <- plotLogRates +
                   labs(title = paste0(self$jobName))
               }
               
               return(plotLogRates)
             },
             
             plotting = function() {
               library(mgcv)
               names(self$data)[1] <- "x"
               self$data <- self$data %>%
                 mutate(
                   q05_smooth = predict(gam(q05 ~ s(x, bs = "cs")), newdata = self$data),
                   q95_smooth = predict(gam(q95 ~ s(x, bs = "cs")), newdata = self$data)
                 )
               plotLogRates <- ggplot(self$data, aes(x = self$predictor)) +
                 geom_ribbon(aes(ymin = q05_smooth, ymax = q95_smooth), fill = "darkgoldenrod1", alpha = 0.5) +
                 # geom_ribbon(aes(ymin = q05, ymax = q95), fill="darkgoldenrod1", alpha = 0.5) +
                 # geom_point(aes(y = median, color = "b"), size = 1.5) + #####################################################################
                 # geom_line(aes(y = q05, color = "c", linetype = "c"), size = 1) +
                
                 # geom_smooth(aes(y = q05, color = "c", linetype = "c"), se = FALSE, size = 0.0, alpha=0.0) +
                 # geom_smooth(aes(y = q95, color = "c", linetype = "c"), se = FALSE,size = 1) +
                 # geom_line(aes(y = q95, color = "c", linetype = "c"), size = 0) +
                 
                 # geom_line(aes(y = median, color = "b", linetype = "b"), size = 1) +
                 geom_smooth(aes(y = median, color = "b", linetype = "b"), method = "gam", se = FALSE,size = 0.5) +
                 geom_line(aes(y = self$other$trueLine, color = "a", linetype = "a"), size = 0.5) +
                 
                 # Legend controls
                 scale_color_manual(
                   name = "Legend",
                   values = c("a" = "blue", "b" = "#F8766D"), 
                   labels = c("LL", "GP")
                 ) +
                 scale_linetype_manual(
                   name = "Legend",
                   values = c("a" = "solid", "b" = "solid"), 
                   labels = c("LL", "GP")
                 ) + 
                 labs(x = paste0(self$predictorName), y = "Log rates") + 
                 theme_classic() +
               theme(
                 panel.grid.major.y = element_line(color = "grey90", linetype = "33", size = 0.5),
                 panel.grid.minor.y = element_blank(),  # optional: remove minor grid lines
                 panel.grid.major.x = element_blank(),  # optional: remove vertical lines if you want
                 axis.ticks.x = element_line(size = 0.3),
                 axis.ticks.y = element_line(size = 0.3),
                 axis.line = element_blank(),
                 panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # frame around plot
               )
               
               if (self$actions$publicationQuality) {
                 fontSize <- 12
                 fontSizeSmall <- 10
                 fontType <- "CMU Serif"
                 
                 if (self$attributes$predictorName == "Host genetic distance") {
                   plotLogRates <- plotLogRates +
                     theme( legend.position = c(0.95, 0.9),
                            legend.justification = c(1, 1))
                 } else {
                   plotLogRates <- plotLogRates +
                     theme( legend.position = c(0, 1),
                            legend.justification = c(0, 1)) 
                  #  xlim(-0.3, 1.5) +
                  # yliim(-4, -1)
                 }
                 
                 
                 plotLogRates <- plotLogRates +
                   theme(
                     # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                     legend.box.just = "right",
                     legend.key.width = unit(0.5, "lines"),
                     legend.spacing.y = unit(0, "cm"),          # reduce space between legend items
                     legend.key.height = unit(0.3, "cm"),        # reduce height of each legend key (line)
                     # legend.margin = margin(0, 0, 0, 0),
                     # legend.margin = margin(6, 6, 6, 6),
                     legend.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),
                     # legend.title = element_text(size = fontSizeSmall, family = fontType, color = "black", face = "bold"),
                     legend.title=element_blank(),
                     legend.background = element_rect(
                       fill = alpha("white", 0))  # Semi-transparent white background
                     # color = "black",             # Border color
                     # size = 0.5)                   # Border thickness
                   ) +
                   theme(text = element_text(size = fontSize, family = fontType, color = "black"),  # Set base text size
                         axis.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),  # Axis tick labels
                         axis.title = element_text(size = fontSize, family = fontType, color = "black")
                   )
               }
               

               return(plotLogRates)
             }, 
             
             plottingFull = function() {
               scale_factor= max(self$dataRates$median)/max(self$data$median) 
               plotLogRates <- ggplot() + 
                 # geom_point(data=self$data, aes(x = self$predictor, y = median*scale_factor, color = "LogRates"), size = 1.5) + 
                 geom_ribbon(data=self$data, aes(x = self$predictor, ymin = q05*scale_factor, ymax = q95*scale_factor), fill = "darkslategray1", alpha = 0.5) +
                 # geom_point(data=self$dataRates, aes(x = self$predictor, y = median, color = "Rates"), size = 1.5) + 
                 geom_ribbon(data=self$dataRates,aes(x = self$predictor, ymin = q05, ymax = q95), fill = "darkgoldenrod1", alpha = 0.5) +
                 # geom_errorbar(aes(ymin = q25, ymax = q75, color = "darkgreen"), width = 0.05) +
                 # geom_errorbar(aes(ymin = q05, ymax = q95, color = "green"), width = 0.1) +
                 #geom_smooth(aes(y = median), method = "lm", se = FALSE, color = "blue", linewidth = 0.4) +
                 #geom_line(aes(y = self$other$trueLine), color = "orange" ) +
                 #geom_text(aes(x = min(self$predictor), y = self$other$ytextSlope), label = self$other$slope, hjust = -0.1, vjust = 1.5, size = 3, color = "blue") +
                 # geom_point(aes(y = trueLogRates, color = "red"), size = 1.5) +
                 #geom_smooth(aes(y = logRates), method = "lm", se = FALSE, color = "red", linewidth = 0.4) +
                 #geom_line(aes(y = logRates, color = "red"), linewidth = 0.4) +
                 labs(x = paste0(self$predictorName), y = "Rates") + 
                 scale_y_continuous(
                   name = "rates",
                   sec.axis = sec_axis(trans=~./scale_factor, name="LogRates")
                 ) + 
                 scale_color_manual(name = "Legend",  # Adding manual color scale with legend
                                    values = c("LogRates" = "#2171B5", "Rates" = "#F8766D"),
                                    labels = c("LogRates", "Rates")) +
                 theme_minimal() +
                 theme(legend.position = "none") 
               # +ylim(self$other$ymin, self$other$ymax)
               # theme(legend.position = "right",
               #       legend.justification = c(1,1),
               #       legend.text = element_text(size = 10), 
               #       legend.title = element_text(size = 10), 
               #       legend.key.size = unit(0.5, 'cm')) +
               
               # cat(self$other$ymin, self$other$ymax, "end")
               if (self$actions$publicationQuality) {
                 fontSize <- 30
                 fontType <- "serif"
                 plotLogRates <- plotLogRates +
                   # theme(legend.position = "bottom", 
                   #       legend.direction = "horizontal",
                   #       legend.text = element_text(size = fontSize, family = fontType, color = "black"),
                   #       legend.title = element_text(size = fontSize, family = fontType, color = "black"),
                   #       legend.key.size = unit(2, 'cm'),
                   #       plot.margin = margin(1, 1, 2, 1, "cm"),  # Increase bottom margin for legend
                   #       legend.box.spacing = unit(1, "cm")
                   # ) +
                   theme(text = element_text(size = fontSize, family = fontType, color = "black"),  # Set base text size
                         axis.text = element_text(size = fontSize, family = fontType, color = "black"),  # Axis tick labels
                         axis.title = element_text(size = fontSize, family = fontType, color = "black")
                   )
                 
               }
               
               return(plotLogRates)
             }
           )
)


plotLogRatesMeansSds <- function(summStat, predictors, jobName) {
  summaryLogRates <- summStat |> row.names() |> 
    sapply(function(x) grepl("logRates", x)) %>% summStat[.,]
  summaryLogRatesPred <- cbind(predictors, summaryLogRates) |> as_tibble()
  plotLogRates <- ggplot(summaryLogRatesPred, aes(x = predictors, y = mean)) + 
    # geom_point(aes(color = "blue")) + 
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = "darkgreen"), width = 0.1) +  # Add error bars
    geom_line(aes(y = predictors, color = "red"), linewidth = 0.4) +
    labs(title = paste0(jobName), x = "Predictors", y = "Log rates (mean/sd)") + 
    scale_color_manual(name = "Legend",  # Adding manual color scale with legend
                       values = c("blue" = "blue", "darkgreen" = "darkgreen", "red" = "red"),
                       labels = c("Mean", "Standard deviation", "True Log rates")) +
    theme_minimal() 
  plotLogRates
}


# boxPlotGlmCoefficients <- function(data, jobName) {
#   glmCoefficients <- data |> names() |> 
#     sapply(function(x) grepl("glmCoefficients", x)) %>% data[.,]
#   if (ncol(glmCoefficients) == 1) {
#     glmCoefficientsLong <- glmCoefficients |> rename(y = glmCoefficients) |> 
#       mutate(glmCoefficients = 1)
#   } else {
#     names(glmCoefficients) <- gsub("glmCoefficients", "", names(glmCoefficients))
#     glmCoefficientsLong <- pivot_longer(glmCoefficients, cols = everything(), 
#                                         names_to = "glmCoefficients", values_to = "y")
#     glmCoefficientsLong$glmCoefficients <- factor(glmCoefficientsLong$glmCoefficients, 
#                                                   levels = unique(glmCoefficientsLong$glmCoefficients))
#   }
#   boxPlot <- ggplot(glmCoefficientsLong, aes(x = glmCoefficients, y = y)) +
#     geom_boxplot(outlier.shape = NA) +
#     labs(title = paste0(jobName), x = "") +
#     scale_y_continuous(limits = quantile(glmCoefficientsLong$y, c(0.01, 0.99)))
#   boxPlot
# }





# --- Final Plots ----
library(sysfonts)
publicationQualitySettings = function(myplot) {
  fontSize <- 12
  fontSizeSmall <- 10
  fontType <- "CMU Serif"
  
  myplot <- myplot +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # frame around plot)
          panel.grid.major.y = element_line(color = "grey90", linetype = "33", size = 0.5),
          panel.grid.minor.y = element_blank(),  # optional: remove minor grid lines
          panel.grid.major.x = element_blank(),  # optional: remove vertical lines if you want
          axis.ticks.x = element_line(size = 0.3),
          axis.ticks.y = element_line(size = 0.3),
          axis.line = element_blank()
    ) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linetype = "33", size = 0.5),
      # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      legend.box.just = "right",
      legend.key.width = unit(0.6, "lines"),
      legend.spacing.y = unit(0, "cm"),          # reduce space between legend items
      legend.key.height = unit(0.3, "cm"),        # reduce height of each legend key (line)
      # legend.margin = margin(0, 0, 0, 0),
      # legend.margin = margin(6, 6, 6, 6),
      legend.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),
      # legend.title = element_text(size = fontSizeSmall, family = fontType, color = "black", face = "bold"),
      legend.title=element_blank(),
      legend.background = element_rect(
        fill = alpha("white", 0))  # Semi-transparent white background
      # color = "black",             # Border color
      # size = 0.5)                   # Border thickness
    ) +
    theme(text = element_text(size = fontSize, family = fontType, color = "black"),  # Set base text size
          axis.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),  # Axis tick labels
          axis.title = element_text(size = fontSize, family = fontType, color = "black")
    )
  return(myplot)
}

otherPublicationQualitySettings = function(myplot) {
  sysfonts::font_add("CMU Serif", regular = "/Users/filippomonti/Library/Fonts/cmunrm.otf")
  
  fontSize <- 12
  fontSizeSmall <- 10
  fontType <- "CMU Serif"
  myplot <- myplot +
    theme_classic() +
    theme(legend.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),
          axis.text = element_text(size = fontSizeSmall, family = fontType, color = "black"),  # Axis tick labels
          axis.title = element_text(size = fontSize, family = fontType, color = "black")
    )
  return(myplot)
}

createPlottingDataLogRates <- function(object, addPredictor = TRUE, addTrueline=FALSE) {
  raw_data <- object$data
  data <- object$summStatCalculator(raw_data)
  summ_data <- data |> row.names() |> 
    sapply(function(x) 
      grepl("logRates", x)) %>% 
    data[.,] |> as_tibble()
  if (addPredictor) {
    predictor <- object$getPredictor()
    summ_data <- cbind(summ_data, predictor)
  }
  
  if (addTrueline) {
    trueLine <- object$otherData$trueLine
    summ_data <- cbind(summ_data, trueLine)
  }
  
  return(summ_data)
}


# --- Rewards Functions ----

reward_clean_robust <- function(single_tree_obj) {
  tree_data <- as_tibble(single_tree_obj)
  phylo_obj <- as.phylo(single_tree_obj)
  
  trunk_nodes <- nodepath(phy = phylo_obj, from = 821, to = 975)
  
  tree_data_filtered <- tree_data %>%
    filter(node %in% trunk_nodes)
  
  ac_columns <- grep("^AC[0-9]+_R$", names(tree_data_filtered), value = TRUE)
  if(length(ac_columns) == 0) {
    warning("No AC_R columns found in filtered tree data. Returning empty tibble.")
    return(tibble())
  }
  
  tree_data_selected <- tree_data_filtered %>%
    dplyr::select(all_of(ac_columns))
  
  names(tree_data_selected) <- sub("_R$", "", names(tree_data_selected))
  
  rewards <- as_tibble(as.list(colSums(tree_data_selected, na.rm = TRUE)))
  return(rewards)
}

rewardSummaryCreator <- function(full_beast_trees) {
  num_trees <- length(full_beast_trees)
  burnIn <- 0
  
  trees_to_process <- full_beast_trees[(burnIn + 1):num_trees]
  
  trunk_rewards <- map_dfr(trees_to_process, reward_clean_robust)
  
  # print(dim(trunk_rewards))
  # print(head(trunk_rewards))
  
  a <- trunk_rewards
  a[] <- lapply(a, as.numeric)
  summary_stats <- a %>%
    pivot_longer(
      cols = everything(), # Select all columns
      names_to = "variable",
      values_to = "value"
    ) %>%
    group_by(variable) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      lower_HPD = quantile(value, 0.01, na.rm = TRUE),
      upper_HPD = quantile(value, 0.99, na.rm = TRUE),
      .groups = 'drop' # Drop the grouping after summarising
    )
  summary_stats <- summary_stats |> rename(aircommunity = variable) 
  if (any(summary_stats$aircommunity %in% names(ac_to_letter))) { #ac_to_letter defined by the dictionary above
    summary_stats$aircommunity <- ac_to_letter[summary_stats$aircommunity]
  }
  summary_stats$aircommunity <- factor(summary_stats$aircommunity, levels = sort(ac_to_letter, decreasing = FALSE))
  return(summary_stats)
}

