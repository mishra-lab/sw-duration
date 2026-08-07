functions {
  vector sum1(vector x){
    return x / sum(x);
  }
  real fE(vector f_i, vector d_i){
    return sum(f_i .* d_i);
  }
  real fCV(vector f_i, vector d_i, real E){
    return sqrt(sum(f_i .* (E - d_i)^2)) / E;
  }
  vector psam(vector z_i, real sp0, real st95){
    return 1 + (sp0 - 1) * exp( z_i * log(.05 / (1 - sp0)) / st95 );
  }
}

data {
  // prior bounds
  real BEu[2];   // full* duration mean E[u]
  real BCVu[2];  // full* duration CV[u]
  real Bsp0[2];  // proportion sampled immediately
  real Bst95[2]; // time to 95% sampled
  // counts
  int Ni; // num integration steps
  int Nm; // num means
  int Nq; // num quantiles
  // data
  real m_n[Nm]; // mean sample sizes
  real m_v[Nm]; // mean estimate
  real q_n[Nq]; // quantile sample sizes
  real q_p[Nq]; // quantile CDF estimate
  int  q_i[Nq]; // quantile value index
  // misc
  vector[Ni] d_i; // integral steps
  real eps; // a small number
  int gps;  // generate posterior samples
  int fam;  // family: {1: gamma, 2: lognormal}
}

parameters {
  real<lower=BEu[1],  upper=BEu[2]>   Eu;  // full* duration mean E[x]
  real<lower=BCVu[1], upper=BCVu[2]>  CVu; // full* duration CV[x]
  real<lower=Bsp0[1], upper=Bsp0[2]>  sp0;  // proportion sampled immediately
  real<lower=Bst95[1],upper=Bst95[2]> st95; // time to 95% sampled
}

transformed parameters {
  real Ex, CVx; // full duration E[x],   CV[x]
  real Ez, CVz; // obs  duration E[z|s], CV[z|s]
  vector[Ni] fx_i, fz_i, Fz_i; // PDF[x],PDF[z|s],CDF[z|s]
  real a, b; // distribution parameters
  // numeric integration
  if (fam==1){ # gamma
    a = 1/CVu^2;
    b = 1/CVu^2/Eu;
    for (i in 1:Ni){
      fx_i[i] = exp(gamma_lpdf(d_i[i] | a, b));
      fz_i[i] = 1 - gamma_cdf (d_i[i] | a, b);
  }}
  if (fam==2){ # lognormal
    a = log(Eu/sqrt(1+CVu^2));
    b = sqrt(log(1+CVu^2));
    for (i in 1:Ni){
      fx_i[i] = exp(lognormal_lpdf(d_i[i] | a, b));
      fz_i[i] = 1 - lognormal_cdf (d_i[i] | a, b);
  }}
  // sampling adjustment & normalization
  fx_i = sum1(fx_i + eps);
  fz_i = sum1(fz_i .* psam(d_i, sp0, st95) + eps);
  // downstream parameters
  Ex   = fE (fx_i, d_i);       // E[x]
  CVx  = fCV(fx_i, d_i, Ex);   // CV[x]
  Ez   = fE (fz_i, d_i);       // E[z|s]
  CVz  = fCV(fz_i, d_i, Ez);   // CV[z|s]
  Fz_i = cumulative_sum(fz_i); // CDF[z|s]
}

model {
  // model vs data
  for (i in 1:Nm){ # mean estimates
    m_v[i] ~ gamma(m_n[i]/CVz^2,
                   m_n[i]/CVz^2/Ez);
  }
  for (i in 1:Nq){ # proportion estimates
    q_p[i] ~ beta(q_n[i]*   Fz_i[q_i[i]],
                  q_n[i]*(1-Fz_i[q_i[i]]));
  }
}

generated quantities {
  // posterior samples
  real m_vs[Nm];
  real q_ps[Nq];
  if (gps){
    for (i in 1:Nm){ # mean estimates
      m_vs[i] = gamma_rng(m_n[i]/CVz^2,
                          m_n[i]/CVz^2/Ez);
    }
    for (i in 1:Nq){ # proportion estimates
      q_ps[i] = beta_rng(q_n[i]*   Fz_i[q_i[i]],
                         q_n[i]*(1-Fz_i[q_i[i]]));
    }
  }
}
