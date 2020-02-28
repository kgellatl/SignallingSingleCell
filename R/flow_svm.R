#' Flow Support Vector Machine
#'
#' This function will take an expression set class after flow sorting, and train to identify cells that failed to pass the gates.
#'
#' @param input the input data ex_sc
#' @export
#' @details
#' This will take an expression set class and classify cells via SVM that did not pass the gates.
#' Note that the SVM parameters are all defaults. You are encourages to modify this function to suit your purposes.
#' This will do SVM training and classification on the Components that were fed INTO tSNE.
#' @examples
#' flow_svm(input = sc_dat)

flow_svm <- function(input, pcnames){
  training <- rownames(pData(input))[which(!is.na(pData(input)$"Pass_Gate"))]
  training <- pData(input)[training,c(grep("Pass_Gate", colnames(pData(input))), grep(pcnames, colnames(pData(input))))]
  trctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 5)
  set.seed(3233)
  svm_Linear <- caret::train(Pass_Gate ~., data = training, method = "svmLinear",
                      trControl=trctrl,
                      tuneLength = 10)
  prediction <- predict(svm_Linear, newdata = pData(input)[,grep(pcnames, colnames(pData(input)))])
  pData(input)$SVM_Classify <- as.character(prediction)
  return(input)
}
