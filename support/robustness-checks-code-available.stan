data {
  int <lower = 1> n_sub;
  int <lower = 0, upper = 6> n_results [n_sub];
  int <lower = 0, upper = 6> x [n_sub];
  int <lower = 0, upper = 1> code [n_sub];
}

parameters{
  real a;
  real b_code;

  real <lower = 0> paper_sigma;
  vector [n_sub] z_a_paper;
}

transformed parameters{
  vector [n_sub] a_paper;
  a_paper = z_a_paper * paper_sigma;
}

model{
  vector [n_sub] q;

  a ~ normal(0, 2.5);

  b_code ~ normal(0, 1);

  paper_sigma ~ exponential(1);

  z_a_paper ~ normal(0, 1);

  for (i in 1:n_sub) {   
      q[i] = inv_logit(a + a_paper[i] + b_code * code[i]);
    }
  
  x ~ binomial(n_results, q);
}
