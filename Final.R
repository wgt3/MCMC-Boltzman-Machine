##2 construct MH sampler with know eta
dd = un_data_small##need to import
#eta = rep(1,28)
indicator= function(i,j){ ##this function computes the product of Xi's
  if(i==j)
    return(1)
  else
    return(-1)
}

gnuf = function(e,eta){
  
  theta_sum1 = 0
  for(i in 1:2){ ##for US and Russia
    for (j in 1:6){##for other countries
     
      theta_sum1 = theta_sum1 +indicator(e[i],e[2+j])*eta[((i-1)*6)+j]
      
    }
  
  theta_sum1 = theta_sum1 + indicator(e[1],e[2])*eta[13]##for us and russia edge
  theta_sum2 = 0
  for(i in 1:2){##for issue edges to us and russia
    for(j in 1:2){
      theta_sum2 = theta_sum2 + indicator(e[i],e[j+8])*eta[13+i+j]
    }
  }
  theta_sum = theta_sum1 +theta_sum2
  alpha_sum = 0
  
  for(i in 1:10){
    alpha_sum = alpha_sum + e[i]*eta[17+i]
  }
  
  return(exp(-theta_sum-alpha_sum))##z will cancel in MH algorithm 
  
  
  }
}

##Discrete Sampler function
DiscreteSampler = function(p,q){
  u = runif(1)
  out = 0
  
  
  
  lowbound = 0
  highbound = 0
  
  for(i in 1: length(p)){
    
    highbound = highbound +p[i]
    #print(highbound)
    lowbound = highbound -p[i] 
    #print(lowbound)
    
    # print(p[i])
    
    
    if (u>lowbound && u<highbound){
      out = q[i]}
    
    
  }
  
  return(out)
}

get_proposal = function(e){
  
  digit = DiscreteSampler(rep(1/10,10),seq(1,10))
  e[digit] = e[digit]*-1 ##flips 1 to -1 and vice versa
  if(e[9]==e[10]){
    e[9]=e[9]*-1 ##makes sure we only have one issue at a time
  }
  return(e)
  
}



etest = c(1,1,-1,1,-1,1,-1,1,-1,1)
gnuf(etest,eta)
##MH sampler
mh_sampler = function(e,t,eta){
  ##initial state
  alpha = 0
  
  for (i in 1:t){
   
    proposal = get_proposal(e)
    
    choice = exp(log(gnuf(proposal,eta))-log(gnuf(e,eta)))##symmetric proposal so q's cancel
    ##use exponential to deal with large numbers and floating point error
    alpha = min(1,choice)
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
  }
  return(e)
}

print(mh_sampler(etest,1000,eta))

##2a
##we need a function that computes x_us*x_ru
##will use the LLN for markov chains to compute prob of this combo
sample_path <- function(e,N,eta)
{
  
  stat <- rep(0, N)
  for (i in 1:N){
    
    proposal = get_proposal(e)
    
    choice = (gnuf(proposal,eta)/gnuf(e,eta))##symmetric proposal so q's cancel
    
    alpha = min(1,choice)
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
    stat[i] = e[1]*e[2]
  }
  
  #plot(1:N, stat, type="l")
  #return (stat)
  return(sum(stat)/N)
}
#print(sample_path(etest,1000))
##try burn in time of about 700 iterations
Expected_value = function(e,eta){
  burn = mh_sampler(e,1000,eta)
  return(sample_path(burn,10000,eta))
  
}
print(Expected_value(etest,eta))

##2b
##modify proposal function so that x_eco is always 1
get_proposal_modified = function(e){
  
  digit = DiscreteSampler(rep(1/8,8),seq(1,8))##x_eco is last index on e
  e[digit] = e[digit]*-1 ##flips 1 to -1 and vice versa
  return(e)
  
}
mh_sampler_modified = function(e,t,eta){
  ##initial state
  alpha = 0
  
  for (i in 1:t){
    
    proposal = get_proposal_modified(e)
    
    choice = (gnuf(proposal,eta)/gnuf(e,eta))##symmetric proposal so q's cancel
    
    alpha = min(1,choice)
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
  }
  return(e)
}
etest = c(rep(-1,9),1) ##x_eco will always be 1 and x_hu will always be -1
#print(mh_sampler_modified(etest,1000))
sample_path_modified <- function(e,N,eta)
{
  
  stat <- rep(0, N)
  for (i in 1:N){
    
    proposal = get_proposal_modified(e)
    
    choice = (gnuf(proposal,eta)/gnuf(e,eta))##symmetric proposal so q's cancel
    
    alpha = min(1,choice)
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
    if(e[1]==e[2])
      stat[i]=1
    else
      stat[i]=0
  }
  
  #plot(1:N, stat, type="l")
  #return (stat)
  return(sum(stat)/N)
}
#print(sample_path_modified(etest,1000))
prob_us_ru = function(e,eta){
  burn = mh_sampler_modified(etest,1000,eta)
  return(sample_path_modified(burn,10000,eta))
  
}
#print(prob_us_ru(etest))


####5 compute gradient
##need a function that outputs a vector of expected values of combos given eta in order
expected_value_vector = function(eta){
  #browser()
  etest = c(1,-1,1,-1,1,-1,1,-1,1,-1)
  burn = mh_sampler(etest,1000,eta)
  e = burn
  theta_list = vector(mode="list", length=27) ##this will store all of the values that need summming
  for(i in 1:27){
    theta_list[[i]] = rep(0,10000)##creates a vector at each place in the list
  }
  for (x in 1:10000){##this loop will fill vectors with datto be averaged
    #browser()
    gnu_e = gnuf(e,eta)
    proposal = get_proposal(e)
    gnu_prop = gnuf(proposal,eta)
    choice = exp(log(gnu_prop)-log(gnu_e))
    
   
    while(is.finite(choice)==FALSE){ ##so alpha doesnt get messed up....
      gnu_e = gnuf(e,eta)
      proposal = get_proposal(e)
      gnu_prop = gnuf(proposal,eta)
      choice = exp(log(gnu_prop)-log(gnu_e))
    }##symmetric proposal so q's cancel
    ##exponential deals with floating point error
    alpha = min(1,choice)
    
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
    ##compute xi*xj
    for (i in 1:2) {##for US and Russia
      for (j in 1:6) {##for other countries
      #browser()
        theta_list[[(i-1)*6+j]][x] = e[i]*e[2+j]
        

      }
      #browser()
    }
    #browser()
    theta_list[[13]][x] = e[1]*e[2]
    for(i in 1:2){##for issue edges to us and russia
      for(j in 1:2){
        #browser()
        theta_list[[13+(i-1)*2+j]][x] = e[i]*e[j+8]
      }
    }
    #browser()
    
    ##compute xi
    for(i in 18:27){
      theta_list[[i]][x]=e[i-17]
    }

    #browser()
    
  }
  expected_values = rep(0,27)
  for (y in 1:27){
     expected_values[y] = sum(theta_list[[y]])/10000
  }
  return(expected_values)
  
  
  
  
}
##lets test this...##at this point makes sense that most would be 0
#print(expected_value_vector(eta))
yes_no = function(a){
  if(a=="yes")
    return(1)
  if(a=="no")
    return(-1)
}

####need to create a function that turns this data into something useful
create_data = function(){
  vote_id = unique(dd$rcid)##creates vector that has all vote id's
  country_id = unique(dd$country_code)
  #browser()
  new_data = matrix(rep(0,10190),nrow = 1019)
  for (i in 1:length(vote_id)){
    
    rcid_index = which(dd$rcid==vote_id[i])
    ##stores lines of data in which data for this vote is
    ##if statements to identify issue
    issue = as.integer(dd[rcid_index[1],6])
    if(issue==4){
      new_data[i,9]=1
      new_data[i,10]=-1
    }
    else{
      new_data[i,10]=1
      new_data[i,9]=-1
    }
    
    
      for(j in 1:length(rcid_index)){##for country's votes
        #going to do seperate if statements cause only 8 countries
        if(dd[rcid_index[j],3]=="US"){
          new_data[i,1] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="RU"){
          new_data[i,2] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="CA"){
          new_data[i,3] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="CH"){
          new_data[i,4] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="ME"){
          new_data[i,5] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="SD"){
          new_data[i,6] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="ss"){
          new_data[i,7] = yes_no(dd[rcid_index[j],4])
        }
        if(dd[rcid_index[j],3]=="TZ"){
          new_data[i,8] = yes_no(dd[rcid_index[j],4])
        }
      
        
    }
      
  }
  ##now we need to run through the matrix and fill in zeroes 
  ##with -1 and 1 with equal probability
  for(i in 1:1019){
    for(j in 1:10){
      if(new_data[i,j]==0)
        new_data[i,j]=DiscreteSampler(c(.5,.5),c(-1,1))
    }
      
  }
  return(new_data)
}
#vote_matrix = create_data()

grad_Log_eta = function(eta){
  expected = 10000*expected_value_vector(eta)
  theta_partial = rep(0,17)
  alpha_partial = rep(0,10)
  #for US edges
    for(k in 3:8){
      sum1 = 0
      for(j in 1: nrow(vote_matrix)){
        sum1 = sum1 + (vote_matrix[j,1]*vote_matrix[j,k])
      }
      theta_partial[k-2] = -sum1 + expected[k-2]
    }
  ##for RU edges
  for(k in 3:8){
    sum1 = 0
    for(j in 1: nrow(vote_matrix)){
      sum1 = sum1 + (vote_matrix[j,2]*vote_matrix[j,k])
    }
    theta_partial[4+k] = -sum1 + expected[4+k]
  }
  # ##for US-RU edges
  sum1 = 0
  for(j in 1: nrow(vote_matrix)){
    sum1 = sum1 + (vote_matrix[j,1]*vote_matrix[j,2])
  }
  theta_partial[13] = -sum1 +expected[13]
  ##for issue edges
  for(i in 1:2){
    for(k in 9:10){
      sum1 = 0
      for(j in 1: nrow(vote_matrix)){
        sum1 = sum1 + (vote_matrix[j,i]*vote_matrix[j,k])
      }
      theta_partial[k+3+2*i] = -sum1 +expected[k+3+2*i]
    }
  #   
  }
  for(i in 1:10){
    sum1 = 0
    for(j in 1: nrow(vote_matrix)){
      sum1 = sum1 + vote_matrix[j,i]
    }
    alpha_partial[i] = -sum1 + expected[i+17]
  }
  return(c(theta_partial,alpha_partial))
}
#print(grad_Log_eta(eta))

####numero sixo, el steepest ascensor
mag = function(a){
  sum1 = 0
  for(i in length(a)){
    sum1 = sum1 + (a[i])^2
  }
  return(sqrt(sum1))
}

log_eta_estimate = function(eta){
  sum1 = 0
  for(i in 1:nrow(vote_matrix)){
    sum1 = sum1 + log(gnuf(vote_matrix[i, ],eta))
  }##note- not going to compute z-just want to test steepest ascent
  return(sum1/1000)#to deal with its growth
}
steepest_ascent = function(t){
  s = .5
  eta = rep(0,27)
  #etavec = rep(0,t)
  for(i in 1:t){
    gradient = grad_Log_eta(eta)
    gradient = gradient/mag(gradient)##normalizes
    eta = eta + s*gradient
    if(i%%50==0){##tracks it every 50 iterations
      est = log_eta_estimate(eta)
      print(est)
      
    }
    
    #etavec[i] = log_eta_estimate(eta)  
    }
  #plot(seq(1,t),etavec)
  return(eta)
    
}
eta = steepest_ascent(1100)
print(eta)


##write some functions that test how good eta is-compare to dataset
##this one will compute the expected value of X_us..this sucks, try somthin else
sample_path_test <- function(e,N,eta)
{
  
  stat <- rep(0, N)
  for (i in 1:N){
    #browser()
    gnu_e = gnuf(e,eta)
    proposal = get_proposal(e)
    gnu_prop = gnuf(proposal,eta)
    
    choice = exp(log(gnu_prop)-log(gnu_e))##symmetric proposal so q's cancel
    
    alpha = min(1,choice)
    #print(alpha)
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
    stat[i] = e[10]
    #print(e)
  }
  
  #plot(1:N, stat, type="l")
  #print(stat)
  return(sum(stat)/N)
}

Expected_value_test = function(e,eta){
  
  burn = mh_sampler(e,1000,eta)
  return(sample_path_test(burn,10000,eta))
  
}
etest = c(rep(1,9),-1)
expected_x_us = function(eta){
  from_data = 0
  from_estimate = 0
  from_data = sum(vote_matrix[ ,10 ])/nrow(vote_matrix)
  from_estimate = Expected_value_test(etest,eta)
  return(cat("From estimate of eta, the expected value of X_us is", from_estimate, ".", 
         "From the data set, the expected value of X_us is", from_data))
  
}
print(expected_x_us(eta))
##this one will compute the expected value of x_us*x_ca
sample_path_test2 <- function(e,N,eta)
{
  
  stat <- rep(0, N)
  for (i in 1:N){
    gnu_e = gnuf(e,eta)
    proposal = get_proposal(e)
    gnu_prop = gnuf(proposal,eta)
    
    choice = exp(log(gnu_prop)-log(gnu_e))##symmetric proposal so q's cancel
    
    alpha = min(1,choice)
    
    if(runif(1)<alpha)
      e = proposal
    else
      e = e
    stat[i] = e[]*e[]
  }
  
  #plot(1:N, stat, type="l")
  #return (stat)
  return(sum(stat)/N)
}

Expected_value_test2 = function(e,eta){
  burn = mh_sampler(e,1000,eta)
  return(sample_path_test2(burn,100,eta))
  
}
etest = rep(1,10)
expected_x_us_ca = function(eta){
  from_data = 0
  from_estimate = 0
  from_data = sum(vote_matrix[,1]*vote_matrix[,2])/nrow(vote_matrix)
  from_estimate = Expected_value_test2(etest,eta)
  return(cat("From estimate of eta, the expected value of X_us*X_ca is", from_estimate, ".", 
             "From the data set, the expected value of X_us*x_ca is", from_data))
  
}
print(expected_x_us_ca(eta))

####numero siete

print(Expected_value(etest,eta))
print(prob_us_ru(etest,eta))
for(i in 3:8){
  print(sum(vote_matrix[,2]*vote_matrix[,i])/nrow(vote_matrix))
}
prob_eco = function(){
  eco_iss = which(vote_matrix[,10]==1)
  prob_sum =0
  for(i in 1:length(eco_iss)){
    if(vote_matrix[eco_iss[i],1]==vote_matrix[eco_iss[i],2]){
      prob_sum=prob_sum+1
  
    }
  }
  return(prob_sum/length(eco_iss))
}
print(prob_eco())
