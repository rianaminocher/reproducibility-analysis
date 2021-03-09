functions{
    matrix cov_GPL3(matrix x, real rho, real alpha, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = 1;
          for (j in (i + 1):N) {
            K[i, j] = rho * exp(-alpha * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = 1;
        return delta * cholesky_decompose(K);
    }
}

data {
  int <lower = 1> n;
  int <lower = 0, upper = 1> d [n];
  int <lower = 0> age [n];
  int <lower = 1, upper = 4> type [n];

  int <lower = 1> n_sub;
  int <lower = 1> index_sub [n_sub];
  int <lower = 0, upper = 6> x [n_sub];
  int <lower = 0, upper = 6> n_results [n_sub];
  int <lower = 0> age_id [n];
  int <lower = 1> n_age;
  matrix[n_age, n_age] age_dist;
}

parameters{
  vector <lower = 0, upper = 1> [4] A;
  vector <lower = 0> [4] b;

  real a;

  real <lower = 0> paper_sigma;

  real <lower = 0, upper = 1> rho;
  real <lower = 0> alpha;
  real <lower = 0> delta;

  vector [n_sub] z_a_paper;
  vector [n_age] z_age;

}

transformed parameters{
  vector [n_age] k_age;
  matrix [n_age, n_age] L;
  vector [n_sub] a_paper;

  L = cov_GPL3(age_dist, rho, alpha, delta);
  k_age = L * z_age;
  
  a_paper = z_a_paper * paper_sigma;
}

model{
  vector [n] p;
  vector [n_sub] q;

  A ~ beta(4, 2);
  b ~ lognormal(-2, 1);

  a ~ normal(0, 2.5);
  paper_sigma ~ exponential(1);

  rho ~ beta(10, 0.5);
  alpha ~ exponential(10);
  delta ~ exponential(1);

  z_a_paper ~ normal(0, 1);
  z_age ~ normal(0, 1);
  
  for (i in 1:n) {
    p[i] = A[type[i]] * exp(-b[type[i]] * age[i]);
  }
  d ~ bernoulli(p);

  
  for (i in 1:n_sub) {   
      q[i] = inv_logit(a + a_paper[i] + k_age[age_id[index_sub[i]]]);
    }
  x ~ binomial(n_results, q);
}


generated quantities{
  vector [n] p;
  vector [n_sub] q;
  
  for (i in 1:n) {
    p[i] = A[type[i]] * exp(-b[type[i]] * age[i]);
  }

  for (i in 1:n_sub) {   
      q[i] = inv_logit(a + a_paper[i] + k_age[age_id[index_sub[i]]]);
  }
}

