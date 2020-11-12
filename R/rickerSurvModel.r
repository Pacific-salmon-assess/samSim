## Temporary file.  This function will eventually be incoprated into samSim package via the stockRecruitModels.r file.


rickerSurvModel <- function(S, a, b, error, ppnAges, gamma, mSurvAtAge) {

  R2 <- ppnAges[1] * S * exp(a-b*S + gamma * mSurvAtAge[1])  * exp(error)
  R3 <- ppnAges[2] * S * exp(a-b*S + gamma * mSurvAtAge[2]) * exp(error)
  R4 <- ppnAges[3] * S * exp(a-b*S + gamma * mSurvAtAge[3]) * exp(error)
  R5 <- ppnAges[4] * S * exp(a-b*S + gamma * mSurvAtAge[4]) * exp(error)
  R6 <- ppnAges[5] * S * exp(a-b*S + gamma * mSurvAtAge[5]) * exp(error)

  RecBY<-sum(R2,R3,R4,R5,R6)

  return(list(R2,R3,R4,R5,R6,RecBY))
}
