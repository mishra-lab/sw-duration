data {
  int Ni; // num integration steps
  int Nm; // num means
  int Nq; // num quantiles
  vector[Ni] z_i; // integral steps
  real m_n[Nm]; // mean sample sizes
  real m_v[Nm]; // mean estimate
  real q_n[Nq]; // quantile sample sizes
  real q_p[Nq]; // quantile CDF estimate
  int  q_i[Nq]; // quantile value index
  real eps;     // a small number
}

parameters {
  real<lower=-1,upper= 3> lEx;  // log true duration mean E[x]
  real<lower=-2,upper= 2> lCVx; // log true duration CV[x]
  real<lower=-4,upper= 2> lti;  // log time to 95% ID
}

transformed parameters {
  real Ex  = exp(lEx);  // true duration mean E[x]
  real CVx = exp(lCVx); // true duration CV[x]
  real ti  = exp(lti);  // time to 95% ID
  real Ez;              // obs duration mean E[z]
  real Vz;              // obs duration variance V[z]
  real CVz;             // obs duration CV[z]
  vector[Ni] fz_i;      // obs duration PDF[z|s]
  vector[Ni] Fz_i;      // obs duration CDF[z|s]
  for (i in 1:Ni){
    fz_i[i] = eps + // stability hack
      (1 - gamma_cdf(z_i[i] | 1/CVx^2, 1/CVx^2/Ex)) .* // p(z|a)
      (1 - exp(z_i[i] * log(.05) / ti)); // p(s|z)
  }
  fz_i = fz_i / sum(fz_i);          // PDF[z|s] normalized
  Ez   = sum(fz_i .* z_i);          // E[z|s]
  Vz   = sum(fz_i .* (Ez - z_i)^2); // V[z|s]
  CVz  = sqrt(Vz)/Ez;               // CV[z|s]
  Fz_i = cumulative_sum(fz_i);      // CDF[z|s]s
}

model {
  for (i in 1:Nm){ // mean estimates
    m_v[i] ~ gamma(m_n[i]/CVz^2, m_n[i]/CVz^2/Ez);
  }
  for (i in 1:Nq){ // quantile CDF estimates
    q_p[i] ~ beta(q_n[i]*   Fz_i[q_i[i]],
                  q_n[i]*(1-Fz_i[q_i[i]]));
  }
}
