data {
  //London
  int n_days_London;
  int n_days2_London;
  int daily_death_London[n_days_London];
  int sero_London[n_days2_London];
  int t_London[n_days_London];
  int t2_London[n_days2_London];
  real Epidemic_Stage_London[n_days_London];
  
  //NorthEast
  int n_days_NorthEastYorkshireHumber;
  int n_days2_NorthEastYorkshireHumber;
  int daily_death_NorthEastYorkshireHumber[n_days_NorthEastYorkshireHumber];
  int sero_NorthEastYorkshireHumber[n_days2_NorthEastYorkshireHumber];
  int t_NorthEastYorkshireHumber[n_days_NorthEastYorkshireHumber];
  int t2_NorthEastYorkshireHumber[n_days2_NorthEastYorkshireHumber];
  real Epidemic_Stage_NorthEastYorkshireHumber[n_days_NorthEastYorkshireHumber];
  
   //NorthWest
  int n_days_NorthWest;
  int n_days2_NorthWest;
  int daily_death_NorthWest[n_days_NorthWest];
  int sero_NorthWest[n_days2_NorthWest];
  int t_NorthWest[n_days_NorthWest];
  int t2_NorthWest[n_days2_NorthWest];
  real Epidemic_Stage_NorthWest[n_days_NorthWest];
  
  //SouthEast
  int n_days_SouthEast;
  int n_days2_SouthEast;
  int daily_death_SouthEast[n_days_SouthEast];
  int sero_SouthEast[n_days2_SouthEast];
  int t_SouthEast[n_days_SouthEast];
  int t2_SouthEast[n_days2_SouthEast];
  real Epidemic_Stage_SouthEast[n_days_SouthEast];
  
  //SouthWest
  int n_days_SouthWest;
  int n_days2_SouthWest;
  int daily_death_SouthWest[n_days_SouthWest];
  int sero_SouthWest[n_days2_SouthWest];
  int t_SouthWest[n_days_SouthWest];
  int t2_SouthWest[n_days2_SouthWest];
  real Epidemic_Stage_SouthWest[n_days_SouthWest];
  
  //Midlands
  int n_days_Midlands;
  int n_days2_Midlands;
  int daily_death_Midlands[n_days_Midlands];
  int sero_Midlands[n_days2_Midlands];
  int t_Midlands[n_days_Midlands];
  int t2_Midlands[n_days2_Midlands];
  real Epidemic_Stage_Midlands[n_days_Midlands];
  
  //EastofEngland
  int n_days_EastofEngland;
  int n_days2_EastofEngland;
  int daily_death_EastofEngland[n_days_EastofEngland];
  int sero_EastofEngland[n_days2_EastofEngland];
  int t_EastofEngland[n_days_EastofEngland];
  int t2_EastofEngland[n_days2_EastofEngland];
  real Epidemic_Stage_EastofEngland[n_days_EastofEngland];
}

parameters {
  real <lower= 0,upper=1> beta;  
  
  //London
  real <lower= 0,upper=1> eta_London;
  real <lower= 0,upper=1> gamma_London;
  
  //NorthEast
  real <lower= 0,upper=1> eta_NorthEastYorkshireHumber;
  real <lower= 0,upper=1> gamma_NorthEastYorkshireHumber;
  
  //NorthWest
  real <lower= 0,upper=1> eta_NorthWest;
  real <lower= 0,upper=1> gamma_NorthWest;
  
  //SouthEast
  real <lower= 0,upper=1> eta_SouthEast;
  real <lower= 0,upper=1> gamma_SouthEast;
  
  //SouthWest
  real <lower= 0,upper=1> eta_SouthWest;
  real <lower= 0,upper=1> gamma_SouthWest;
  
  //Midlands
  real <lower= 0,upper=1> eta_Midlands;
  real <lower= 0,upper=1> gamma_Midlands;
  
  //EastofEngland
  real <lower= 0,upper=1> eta_EastofEngland;
  real <lower= 0,upper=1> gamma_EastofEngland;
}

model {  
  real store_value;
  int index;
  //beta~ beta(1,1);
  beta~ weibull(3.67,130.4);
  
  //London
  eta_London~  beta(1,1);  
  gamma_London~ beta(1,1); 
 
  //NorthEast
  eta_NorthEastYorkshireHumber~ beta(1,1);
  gamma_NorthEastYorkshireHumber~ beta(1,1); 

  //NorthWest
  eta_NorthWest~ beta(1,1); 
  gamma_NorthWest~ beta(1,1);
  
  //SouthWest
  eta_SouthWest~ beta(1,1); 
  gamma_SouthWest~ beta(1,1);
  
  //SouthEast
  eta_SouthEast~ beta(1,1); 
  gamma_SouthEast~ beta(1,1);
  
  //Midlands
  eta_Midlands~ beta(1,1); 
  gamma_Midlands~ beta(1,1);
  
  //EastofEngland
  eta_EastofEngland~ beta(1,1); 
  gamma_EastofEngland~ beta(1,1);
  
  //London
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_London) {
    for (n in index:t2_London[i]){
      store_value = store_value + exp(beta*t_London[n])*(1-(gamma_London-eta_London*gamma_London*Epidemic_Stage_London[n]))/(gamma_London-(eta_London*gamma_London)*Epidemic_Stage_London[n])*daily_death_London[n];
    }
    sero_London[i] ~ neg_binomial_2(store_value/exp(beta*t2_London[i]),100);
    index = t2_London[i]+1;
  };
  
  //NorthEast
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_NorthEastYorkshireHumber) {
    for (n in index:t2_NorthEastYorkshireHumber[i]){
      store_value = store_value + exp(beta*t_NorthEastYorkshireHumber[n])*(1-(gamma_NorthEastYorkshireHumber-(eta_NorthEastYorkshireHumber*gamma_NorthEastYorkshireHumber)*Epidemic_Stage_NorthEastYorkshireHumber[n]))/(gamma_NorthEastYorkshireHumber-(eta_NorthEastYorkshireHumber*gamma_NorthEastYorkshireHumber)*Epidemic_Stage_NorthEastYorkshireHumber[n])*daily_death_NorthEastYorkshireHumber[n];
    }
    sero_NorthEastYorkshireHumber[i] ~ neg_binomial_2(store_value/exp(beta*t2_NorthEastYorkshireHumber[i]),100);
    index = t2_NorthEastYorkshireHumber[i]+1;
  };
  
  //NorthWest
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_NorthWest) {
    for (n in index:t2_NorthWest[i]){
      store_value = store_value + exp(beta*t_NorthWest[n])*(1-(gamma_NorthWest-(gamma_NorthWest*eta_NorthWest)*Epidemic_Stage_NorthWest[n]))/(gamma_NorthWest-(gamma_NorthWest*eta_NorthWest)*Epidemic_Stage_NorthWest[n])*daily_death_NorthWest[n];
    }
    sero_NorthWest[i] ~ neg_binomial_2(store_value/exp(beta*t2_NorthWest[i]),100);
    index = t2_NorthWest[i]+1;
  };
  
  //SouthEast
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_SouthEast) {
    for (n in index:t2_SouthEast[i]){
      store_value = store_value + exp(beta*t_SouthEast[n])*(1-(gamma_SouthEast-(gamma_SouthEast*eta_SouthEast)*Epidemic_Stage_SouthEast[n]))/(gamma_SouthEast-(gamma_SouthEast*eta_SouthEast)*Epidemic_Stage_SouthEast[n])*daily_death_SouthEast[n];
    }
    sero_SouthEast[i] ~ neg_binomial_2(store_value/exp(beta*t2_SouthEast[i]),100);
    index = t2_SouthEast[i]+1;
  };
  
  //SouthWest
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_SouthWest) {
    for (n in index:t2_SouthWest[i]){
      store_value = store_value + exp(beta*t_SouthWest[n])*(1-(gamma_SouthWest-(gamma_SouthWest*eta_SouthWest)*Epidemic_Stage_SouthWest[n]))/(gamma_SouthWest-(gamma_SouthWest*eta_SouthWest)*Epidemic_Stage_SouthWest[n])*daily_death_SouthWest[n];
    }
    sero_SouthWest[i] ~ neg_binomial_2(store_value/exp(beta*t2_SouthWest[i]),100);
    index = t2_SouthWest[i]+1;
  };
  
  //Midlands
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_Midlands) {
    for (n in index:t2_Midlands[i]){
      store_value = store_value + exp(beta*t_Midlands[n])*(1-(gamma_Midlands-(gamma_Midlands*eta_Midlands)*Epidemic_Stage_Midlands[n]))/(gamma_Midlands-(gamma_Midlands*eta_Midlands)*Epidemic_Stage_Midlands[n])*daily_death_Midlands[n];
    }
    sero_Midlands[i] ~ neg_binomial_2(store_value/exp(beta*t2_Midlands[i]),100);
    index = t2_Midlands[i]+1;
  };
  
  //EastofEngland
  store_value = 0;
  index = 1;
  for (i in 1:n_days2_EastofEngland) {
    for (n in index:t2_EastofEngland[i]){
      store_value = store_value + exp(beta*t_EastofEngland[n])*(1-(gamma_EastofEngland-(gamma_EastofEngland*eta_EastofEngland)*Epidemic_Stage_EastofEngland[n]))/(gamma_EastofEngland-(gamma_EastofEngland*eta_EastofEngland)*Epidemic_Stage_EastofEngland[n])*daily_death_EastofEngland[n];
    }
    sero_EastofEngland[i] ~ neg_binomial_2(store_value/exp(beta*t2_EastofEngland[i]),100);
    index = t2_EastofEngland[i]+1;
  };
}

 


