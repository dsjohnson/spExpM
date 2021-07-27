arma::mat SPS_exp_spQ_v(const arma::sp_mat Q, const arma::mat v, double prec, bool ren=true, bool tt=true, int *pm=NULL);
arma::mat SPS_v_exp_spQ(const arma::mat v, const arma::sp_mat Q, double prec, bool ren=true, bool tt=true, int *pm=NULL);
arma::mat MUSPS_exp_spQ_v(arma::sp_mat Q, arma::vec v, double prec, arma::vec ts);
arma::mat MUSPS_v_exp_spQ(arma::rowvec v, arma::sp_mat Q, double prec, arma::vec ts);
int get_m(double rho, double prec, int mlo=-1);

