#logistic regression for a specific variable
log_reg = function(data, index) {
  logreg=glm(formula=data[,index]~.,data = data[,-index], family=binomial(link="logit"))
  coefs = logreg$coefficients
  coefs[-1] = coefs[-1]/2
  return(coefs)
}

#helper function for the bootstrap standard errors of the logistic regression
bootstrap_lr = function(data, indices) {
  d = data[indices,]
  logreg=glm(formula=d[,1]~., data=d[,-1], family=binomial(link="logit"))
  coefs = logreg$coefficients
  coefs[-1] = coefs[-1]/2
  return(coefs)
}