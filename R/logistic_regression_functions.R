#logistic regression for a specific variable
log_reg = function(data, index) {
  logreg=glm(formula=data[,index]~.,data = data[,-index], family=binomial(link="logit"))
  coefs = logreg$coefficients
  coefs[-1] = coefs[-1]/2
  return(coefs)
}

bootstrap_var = function(data, index) {
  boot_data = data[sample(nrow(data), 1e3, replace=T), ]
  logreg = glm(formula=boot_data[, index]~., data=boot_data[,-index], family=binomial(link="logit"))
  ses = summary(logreg)$coefficients[,2]
  vars = ses ^ 2
  return(vars)
}