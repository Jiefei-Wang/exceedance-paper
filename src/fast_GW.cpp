#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP C_get_range_by_bound(SEXP R_sx, SEXP R_l, SEXP R_h) {
	double* sx = REAL(R_sx);
	double* l = REAL(R_l);
	double* h = REAL(R_h);
	R_xlen_t m = XLENGTH(R_sx);
	R_xlen_t n = XLENGTH(R_l);
	NumericVector L(n);
	NumericVector H(n);
	int index_x_l = 0;
	int index_x_h = m-1;
	for (auto i = 0; i < n; i++) {
		//Rprintf("i:%d,l:%d\n", i, index_x_l);
		while (sx[index_x_l] <= l[i] || index_x_l < i) {
			index_x_l += 1;
			if (index_x_l == m) {
				return(R_NilValue);
			}
		}
		L[i] = index_x_l + 1;

		auto j = n - i - 1;
		//Rprintf("j:%d,h:%d\n",j , index_x_h);
		while (sx[index_x_h] > h[j] || m - index_x_h - 1 < i) {
			if (index_x_h == 0) {
				return(R_NilValue);
			}
			index_x_h -= 1;
		}
		H[j] = index_x_h + 1;
		index_x_h -= 1;
	}
	//Rf_PrintValue(L);
	//Rf_PrintValue(H);
	for (auto i = 0; i < n; i++) {
		if (L[i] > H[i]) {
			return(R_NilValue);
		}
	}
	return List::create(Named("P") = L, Named("Q") = H);
}

// [[Rcpp::export]]
double C_GW_compute_FDR(SEXP sorted_i, SEXP R_P, SEXP R_Q,int rj_num,int n) {
	double* P = REAL(R_P);
	double* Q = REAL(R_Q);
	//int* sorted_i = INTEGER(R_sorted_i);
	int i = 0;
	int j = 0;
	int FP = 0;
	while (i<rj_num&&j<n) {
		int cur_elt = INTEGER_ELT(sorted_i, i);
		if (cur_elt >= P[j] && cur_elt <= Q[j]) {
			FP += 1;
			i += 1;
			j += 1;
		}
		else {
			if (cur_elt < P[j]) {
				i += 1;
			}
			else {
				if (cur_elt > Q[j]) {
					j += 1;
				}
			}
		}
	}
	return (double)FP / rj_num;
}



// [[Rcpp::export]]
SEXP C_get_range_by_bound2(SEXP R_sx, SEXP R_l, SEXP R_h) {
	double* sx = REAL(R_sx);
	double* l = REAL(R_l);
	double* h = REAL(R_h);
	size_t m = XLENGTH(R_sx);
	size_t n = XLENGTH(R_l);
	NumericVector P(n);
	NumericVector Q(n);
	size_t index_x_l = 0;
	size_t index_x_h = m-1;
	for (size_t i = 0; i < n; i++) {
        // Computing P
		while (sx[index_x_l] < l[i]) {
			if (index_x_l == m - 1) {
				return(R_NilValue);
			}
			index_x_l += 1;
		}
        // The index is 1-based
		P[i] = index_x_l + 1;
        if (index_x_l == m - 1 && i != n-1) {
				return(R_NilValue);
		}
        index_x_l++;

        // Computing Q
		size_t j = n - i - 1;
		while (sx[index_x_h] > h[j]) {
			if (index_x_h == 0) {
				return(R_NilValue);
			}
			index_x_h -= 1;
		}
        // The index is 1-based
		Q[j] = index_x_h + 1;
        if (index_x_h == 0 && i != n-1) {
				return(R_NilValue);
		}
		index_x_h --;
	}
	for (size_t i = 0; i < n; i++) {
		if (P[i] > Q[i]) {
			return(R_NilValue);
		}
	}
	return List::create(Named("P") = P, Named("Q") = Q);
}