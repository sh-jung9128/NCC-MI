# Functions for ncc-ipw
risk_set_f <- function(data, X, Delta){
  
  exit_time <- data$X
  FT <- unique(sort(exit_time[data$Delta == 1])) # Store the unique values of the time points where events (failures) occurred, sorted in ascending order, in FT (removing duplicate event times)
  
  risk_set <- data.frame(cbind(FT, t(sapply(FT, function(x){c(sum(exit_time == x), sum(exit_time >= x))}))))
  colnames(risk_set)<-c("failure_time","cases","at_risk")
  return(risk_set)
}

inclusion_prob_f <- function(data, X, Delta, control_size, risk_set){ 
  
  exit_time <- data$X
  
  CS <- risk_set$cases; AR <- risk_set$at_risk; FT <- risk_set$failure_time
  inclusion.prob <- sapply(exit_time, function(x){
    p <- pmin(1, control_size*CS/(AR-CS))
    p[FT > x] <- 0 # zero when not at risk
    1-prod(1-p) #inclusion prob
  })
  inclusion.prob[data$Delta == 1] <- 1 # Data where events have occurred have an inclusion probability of 1
  
  return(inclusion.prob)
}

ncc_data_f <- function(data, control_size = 1) {
  risk_set <- risk_set_f(data = data, X = "X", Delta = "Delta")  # Calculate the risk set
  data$inc.prob <- inclusion_prob_f(data = data, X = "X", Delta = "Delta", control_size = control_size, risk_set = risk_set)  # Calculate the Inclusion probability
  
  # Control sampling and NCC data generation
  case_data <- data[data$Delta == 1, ]
  case_data$weight <- 1 # Cases are always included, so the weight = 1
  ncc_data <- case_data # The initial NCC data starts with case data
  
  for (i in 1:nrow(case_data)) {
    # Event occurrence time
    failure_time <- case_data$X[i]
    
    # Define the risk set for the current event occurrence time
    risk_set_current <- data[data$X >= failure_time & data$Delta == 0, ]
    
    if (nrow(risk_set_current) >= control_size) {
      sampled_controls <- risk_set_current[sample(1:nrow(risk_set_current), control_size), ]
      sampled_controls$weight <- 1 / sampled_controls$inc.prob  # calculate the weights
      ncc_data <- rbind(ncc_data, sampled_controls)
    }
  }
  
  return(ncc_data)
}