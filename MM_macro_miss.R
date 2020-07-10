# Estimate MM model in misspecified case


#-------------------------------------------------------------
# Initialize quantities
#-------------------------------------------------------------

set.seed(123)

# multinomial parameters initialization at true values
theta = array(0, dim = c(data$p1, 2, 4))

for(j in 1:5){
    for(h in 1:2){
        for(d in 1:4){
         theta[j,h,d] = true_parameter$kern[d,j,h]
        }
    }
}
# observed data
obs           = array(0, dim = c(data$n, data$p1, 1, 1))
obs[, , 1, 1] = data$data[,1:data$p1] -1

initial = mixedMemModel(Total = data$n, J = data$p1, Rj = rep(1,data$p1),
Nijr = array(1, dim = c(data$n, data$p1,1)), K = 2, Vj = rep(4,data$p1), alpha = rep(0.2,2),
theta = theta, dist = rep('multinomial',data$p1), obs = obs)

#-------------------------------------------------------------
## fitting the model
#-------------------------------------------------------------

out1 =  mmVarFit(initial, printStatus = 1, printMod = 25)


lambda.point1 = out1$phi / rowSums(out1$phi)


##########################################
##########################################
## model 2
##########################################
##########################################

set.seed(123)

# multinomial parameters initialization at true values
theta = array(0, dim = c(data$p2, 2, 4))

for(j in 6:10){
    for(h in 1:2){
        for(d in 1:4){
         theta[j-5,h,d] = true_parameter$kern[d,j,h]
        }
    }
}
# observed data
obs = array(0, dim = c(data$n, 5, 1, 1))
obs[, , 1, 1] = data$data[,6:10] -1

initial = mixedMemModel(Total = data$n, J = data$p2, Rj = rep(1,data$p2),
Nijr = array(1, dim = c(data$n, data$p2,1)), K = 2, Vj = rep(4,data$p2), alpha = rep(0.2,2),
theta = theta, dist = rep('multinomial',data$p2), obs = obs)

#-------------------------------------------------------------
## fitting the model
#-------------------------------------------------------------

out2  = mmVarFit(initial, printStatus = 1, printMod = 25)
lambda.point2  = out2$phi / rowSums(out2$phi)




##########################################
##########################################
## MODEL 1 with 4 pure types --- (correctly specified) 
##########################################
##########################################



#-------------------------------------------------------------
# Initialize quantities
#-------------------------------------------------------------

set.seed(123)

# multinomial parameters initialization at true values
theta = array(0, dim = c(data$p1, 4, 4))

for(j in 1:5){
    for(h in 1:4){
        for(d in 1:4){
         theta[j,h,d] = true_parameter$kern[d,j,h]
        }
    }
}
# observed data
obs           = array(0, dim = c(data$n, data$p1, 1, 1))
obs[, , 1, 1] = data$data[,1:data$p1] -1

initial = mixedMemModel(Total = data$n, J = data$p1, Rj = rep(1,data$p1),
Nijr = array(1, dim = c(data$n, data$p1,1)), K = 4, Vj = rep(4,data$p1), alpha = rep(0.5,4),
theta = theta, dist = rep('multinomial',data$p1), obs = obs)

#-------------------------------------------------------------
## fitting the model
#-------------------------------------------------------------

out4 =  mmVarFit(initial, printStatus = 1, printMod = 25)
lambda.point4 = out4$phi / rowSums(out4$phi)
