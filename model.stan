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
  int <lower = 1, upper = 3> n_stages;
  int <lower = 0, upper = 6> x [4, n_sub];
  int <lower = 0> age_id [n];
  int <lower = 1> n_age;
  matrix[n_age, n_age] age_dist;
}

parameters{
  vector <lower = 0, upper = 1> [4] A;
  vector <lower = 0> [4] b;

  vector [n_stages] a;

  vector <lower = 0> [n_stages] paper_sigma;

  vector <lower = 0, upper = 1> [n_stages] rho;
  vector <lower = 0> [n_stages] alpha;
  vector <lower = 0> [n_stages] delta;

  vector [n_sub] z_a_paper [n_stages];
  vector [n_age] z_age [n_stages];

}

transformed parameters{
  vector [n_age] k_age [n_stages];
  matrix [n_age, n_age] L [n_stages];
  vector [n_sub] a_paper [n_stages];

  for (j in 1:n_stages) {
    L[j] = cov_GPL3(age_dist, rho[j], alpha[j], delta[j]);
    k_age[j] = L[j] * z_age[j];
    a_paper[j] = z_a_paper[j] * paper_sigma[j];
  }
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

  for (j in 1:n_stages) {
    z_a_paper[j] ~ normal(0, 1);
    z_age[j] ~ normal(0, 1);
  }

  for (i in 1:n) {
    p[i] = A[type[i]] * exp(-b[type[i]] * age[i]);
  }
  d ~ bernoulli(p);

  for (j in 1:n_stages) { 
    for (i in 1:n_sub) {   
      q[i] = inv_logit(a[j] + a_paper[j, i] + k_age[j, age_id[index_sub[i]]]);
    }
    x[j+1] ~ binomial(x[j], q);
  }
}

generated quantities{
  vector [n] p;
  vector [n_sub] q[n_stages];
  
  for (i in 1:n) {
    p[i] = A[type[i]] * exp(-b[type[i]] * age[i]);
  }

  for (j in 1:n_stages) { 
    for (i in 1:n_sub) {   
      q[j, i] = inv_logit(a[j] + a_paper[j, i] + k_age[j, age_id[index_sub[i]]]);
    }
  }
}

