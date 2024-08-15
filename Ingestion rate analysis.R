food = read.csv("C:/RData/areej_data.csv")
food = food[order(food$food),]
head(food)

library("bbmle")
library("deSolve")

Ingestion_model =function(t, y, parameters) { 
  Food = y[1]
  with(as.list(parameters),{                 # The with() to keep equations tidy

    # Dynamics
    dFooddt = -d*Food - a*L^2*Food/(1 + a*L^2*h*Food)
    
    result = c(dFooddt)          # Store the result in order of y-vector
    return(list(result))
  }
  ) 
} 

timespan=0:1
inits = c(Food = 0.033)
parms = c(d = 0.01, a = 0.0001, h = 0.1, L = 10)
ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda")

Ingestion_predictions = function(df=food, d=d, a=a, h=h){
  predictions = numeric()
  timespan = 0:3
  for(i in 1:dim(df)[1]){
    inits = c(Food = df[i, "initial_weight"])
    parms = c(d = d, a = a, h = h, L = df[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda")[4,2])
  }
  predictions
}

plot(Ingestion_predictions(d = 0.07, a = 0.01, h = 50), food$final_weight)
abline(0, 1)

Ingestion_NLL = function(df=food, d, a, h, sd){
  predictions = numeric()
  timespan = 0:3
  for(i in 1:dim(df)[1]){
    inits = c(Food = df[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(a), h = abs(h), L = df[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  #print(predictions)
  #print((df$final_weight+0.0001)*1000)
  NLL = -sum(dnorm(x=log((df$final_weight+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  #print(NLL)
  NLL
}

Ingestion_NLL(d = 0.07, a = 0.01, h = 0.05, sd=0.1)

m1 = mle2(minuslogl = Ingestion_NLL, start=list(d = 0.07, a = 0.01, h = 0.05, sd=0.1),
          control=list(parscale = c(d = 0.07, a = 0.01, h = 0.05, sd=0.1), trace=TRUE, maxit=1e3))

m1
saveRDS(m1, file="C:/RData/areej_m1.RDA")
coef(m1)



Ingestion_predictions = function(df=food, d=d, a=a, h=h){
  predictions = numeric()
  timespan = 0:3
  for(i in 1:dim(df)[1]){
    inits = c(Food = df[i, "initial_weight"])*1000
    parms = c(d = d, a = a, h = h, L = df[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda")[4,2])
  }
  predictions
}

plot(Ingestion_predictions(d=as.numeric(coef(m1)[1]), a=as.numeric(coef(m1)[2]), h=as.numeric(coef(m1)[3])), food$final_weight*1000)
abline(0, 1)





Ingestion_NLL_d = function(df=food, dA, dB, dC, dD, a, h, sd){
  predictions = numeric()
  timespan = 0:3
  df_A=subset(df, food=="A")
  df_B=subset(df, food=="B")
  df_C=subset(df, food=="C")
  df_D=subset(df, food=="D")
  predictions = numeric()
  for(i in 1:dim(df_A)[1]){
    inits = c(Food = df_A[i, "initial_weight"])*1000
    parms = c(d = abs(dA), a = abs(a), h = abs(h), L = df_A[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_A = -sum(dnorm(x=log((df_A[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))

  predictions = numeric()
  for(i in 1:dim(df_B)[1]){
    inits = c(Food = df_B[i, "initial_weight"])*1000
    parms = c(d = abs(dB), a = abs(a), h = abs(h), L = df_B[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_B = -sum(dnorm(x=log((df_B[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))

  predictions = numeric()
  for(i in 1:dim(df_C)[1]){
    inits = c(Food = df_C[i, "initial_weight"])*1000
    parms = c(d = abs(dC), a = abs(a), h = abs(h), L = df_C[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_C = -sum(dnorm(x=log((df_C[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))

    predictions = numeric()
  for(i in 1:dim(df_D)[1]){
    inits = c(Food = df_D[i, "initial_weight"])*1000
    parms = c(d = abs(dD), a = abs(a), h = abs(h), L = df_D[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_D = -sum(dnorm(x=log((df_D[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))

  NLL_A + NLL_B + NLL_C + NLL_D
}

Ingestion_NLL_d(dA = 0.07, dB = 0.07, dC = 0.07, dD = 0.07, a = 0.01, h = 0.05, sd=0.1)
Ingestion_NLL(d = 0.07, a = 0.01, h = 0.05, sd=0.1)



m1_d = mle2(minuslogl = Ingestion_NLL_d, start=list(dA = 0.14, dB = 0.14, dC = 0.14, dD = 0.14, a = 0.02, h = 0.15, sd=0.5),
           control=list(parscale = c(dA = 0.14, dB = 0.14, dC = 0.14, dD = 0.14, a = 0.02, h = 0.15, sd=0.5), trace=TRUE, maxit=1e3))

m1_d
saveRDS(m1_d, file="C:/RData/areej_m1_d.RDA")
coef(m1_d)

#### Need to do m1_a and m1_ad ####

Ingestion_NLL_a = function(df=food, d, aA, aB, aC, aD, h, sd){
  predictions = numeric()
  timespan = 0:3
  df_A=subset(df, food=="A")
  df_B=subset(df, food=="B")
  df_C=subset(df, food=="C")
  df_D=subset(df, food=="D")
  predictions = numeric()
  for(i in 1:dim(df_A)[1]){
    inits = c(Food = df_A[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(aA), h = abs(h), L = df_A[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_A = -sum(dnorm(x=log((df_A[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_B)[1]){
    inits = c(Food = df_B[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(aB), h = abs(h), L = df_B[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_B = -sum(dnorm(x=log((df_B[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_C)[1]){
    inits = c(Food = df_C[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(aC), h = abs(h), L = df_C[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_C = -sum(dnorm(x=log((df_C[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_D)[1]){
    inits = c(Food = df_D[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(aD), h = abs(h), L = df_D[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_D = -sum(dnorm(x=log((df_D[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  NLL_A + NLL_B + NLL_C + NLL_D
}

Ingestion_NLL_a(d = 0.14, aA = 0.03, aB = 0.03, aC = 0.03, aD = 0.03, h = 0.15, sd=0.5)



m1_a = mle2(minuslogl = Ingestion_NLL_a, start=list(d = 0.14, aA = 0.03, aB = 0.03, aC = 0.03, aD = 0.03, h = 0.15, sd=0.5),
            control=list(parscale = c(d = 0.14, aA = 0.03, aB = 0.03, aC = 0.03, aD = 0.03, h = 0.15, sd=0.5), trace=TRUE, maxit=1e3))

m1_a
saveRDS(m1_a, file="C:/RData/areej_m1_a.RDA")
coef(m1_a)
warnings()

############### Actually, I think we want m1_h #########################


Ingestion_NLL_h = function(df=food, d, a, hA, hB, hC, hD, sd){
  predictions = numeric()
  timespan = 0:3
  df_A=subset(df, food=="A")
  df_B=subset(df, food=="B")
  df_C=subset(df, food=="C")
  df_D=subset(df, food=="D")
  predictions = numeric()
  for(i in 1:dim(df_A)[1]){
    inits = c(Food = df_A[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(a), h = abs(hA), L = df_A[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_A = -sum(dnorm(x=log((df_A[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_B)[1]){
    inits = c(Food = df_B[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(a), h = abs(hB), L = df_B[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_B = -sum(dnorm(x=log((df_B[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_C)[1]){
    inits = c(Food = df_C[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(a), h = abs(hC), L = df_C[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_C = -sum(dnorm(x=log((df_C[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_D)[1]){
    inits = c(Food = df_D[i, "initial_weight"])*1000
    parms = c(d = abs(d), a = abs(a), h = abs(hD), L = df_D[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_D = -sum(dnorm(x=log((df_D[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  NLL_A + NLL_B + NLL_C + NLL_D
}

Ingestion_NLL_h(d = 0.14, a = 0.03, hA = 0.12, hB = 0.12, hC = 0.12, hD = 0.12, sd=0.5)



m1_h = mle2(minuslogl = Ingestion_NLL_h, start=list(d = 0.14, a = 0.03, hA = 0.12, hB = 0.12, hC = 0.12, hD = 0.12, sd=0.5),
            control=list(parscale = c(d = 0.14, a = 0.03, hA = 0.12, hB = 0.12, hC = 0.12, hD = 0.12, sd=0.5), trace=TRUE, maxit=1e3))

m1_h
saveRDS(m1_h, file="C:/RData/areej_m1_h.RDA")
coef(m1_h)

AICtab(m1, m1_d, m1_h, sort=T, delta=T, weights=T)




Ingestion_NLL_hd = function(df=food, dA, dB, dC, dD, a, hA, hB, hC, hD, sd){
  predictions = numeric()
  timespan = 0:3
  df_A=subset(df, food=="A")
  df_B=subset(df, food=="B")
  df_C=subset(df, food=="C")
  df_D=subset(df, food=="D")
  predictions = numeric()
  for(i in 1:dim(df_A)[1]){
    inits = c(Food = df_A[i, "initial_weight"])*1000
    parms = c(d = abs(dA), a = abs(a), h = abs(hA), L = df_A[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_A = -sum(dnorm(x=log((df_A[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_B)[1]){
    inits = c(Food = df_B[i, "initial_weight"])*1000
    parms = c(d = abs(dB), a = abs(a), h = abs(hB), L = df_B[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_B = -sum(dnorm(x=log((df_B[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_C)[1]){
    inits = c(Food = df_C[i, "initial_weight"])*1000
    parms = c(d = abs(dC), a = abs(a), h = abs(hC), L = df_C[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_C = -sum(dnorm(x=log((df_C[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  predictions = numeric()
  for(i in 1:dim(df_D)[1]){
    inits = c(Food = df_D[i, "initial_weight"])*1000
    parms = c(d = abs(dD), a = abs(a), h = abs(hD), L = df_D[i, "major"])
    predictions[i] = as.numeric(ode(y=inits, times = timespan, parms = parms, func = Ingestion_model, method="lsoda",
                                    maxsteps=1e6)[4,2])
  }
  NLL_D = -sum(dnorm(x=log((df_D[,"final_weight"]+0.0001)*1000), mean=log(predictions), sd=sd, log=T))
  
  NLL_A + NLL_B + NLL_C + NLL_D
}

Ingestion_NLL_hd(dA = 0.14, dB = 0.14, dC = 0.14, dD = 0.14, a = 0.03, hA = 0.12, hB = 0.12, hC = 0.12, hD = 0.12, sd=0.5)



m1_hd = mle2(minuslogl = Ingestion_NLL_hd, start=list(dA = 0.14, dB = 0.14, dC = 0.14, dD = 0.14, a = 0.03, hA = 0.12, hB = 0.12, hC = 0.12, hD = 0.12, sd=0.5),
            control=list(parscale = c(dA = 0.14, dB = 0.14, dC = 0.14, dD = 0.14, a = 0.03, hA = 0.12, hB = 0.12, hC = 0.12, hD = 0.12, sd=0.5), trace=TRUE, maxit=1e3))

m1_hd
saveRDS(m1_hd, file="C:/RData/areej_m1_hd.RDA")
coef(m1_hd)


AICtab(m1, m1_d, m1_h,m1_hd, sort=T, delta=T, weights=T)



anova(m1_d, m1_hd)

# profiles
p1_d = profile(m1_d)
p1_hd = profile(m1_hd)

saveRDS(p1_d, file="C:/RData/areej_p1_d.RDA")
saveRDS(p1_hd, file="C:/RData/areej_p1_hd.RDA")
