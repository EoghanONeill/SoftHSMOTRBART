
# include <RcppArmadillo.h>
# include <cmath>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec get_children_cpp(arma::mat tree_mat, int parent) {

  // arma::vec all_children;
  // termcol = tree_mat.col(0);

  if(tree_mat(parent-1, 0) == 1){
    arma::uvec outputvec = {  parent};

    return(outputvec);

  }else{

    int curr_child_left = tree_mat(parent - 1, 1);
    int curr_child_right = tree_mat(parent - 1, 2);


    arma::uvec lefts_vec = get_children_cpp(tree_mat, curr_child_left);
    arma::uvec rights_vec = get_children_cpp(tree_mat, curr_child_right);

    arma::uvec outputvec = arma::join_cols(lefts_vec, lefts_vec );

    // arma::uvec outputvec = {lefts_vec, rights_vec };

    // arma::vec outputvec = { curr_child_left, curr_child_right};
    return(outputvec);

  }
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat phi_app_soft(arma::mat X_stand, arma::mat anc, double tau) {


  arma::vec ancsplits = anc.col(2);
  arma::vec anclefts = anc.col(3);
  arma::uvec leftinds = arma::find(anclefts);
  arma::vec ancterms = arma::unique(anc.col(0));
  arma::vec anctermcol = anc.col(0);
  arma::uvec ancvars = arma::conv_to<arma::uvec>::from(anc.col(4) - 1) ;

  arma::mat phimat(X_stand.n_rows, ancterms.n_elem);

  // Rcpp::Rcout << "Line 15. leftinds = " << leftinds  << ". \n" ;


  for(unsigned int rowind=0; rowind < X_stand.n_rows; rowind++){


    // Rcpp::Rcout << "Line 15. rowind = " << rowind  << ". \n" ;


    arma::vec xrow = (X_stand.row(rowind)).t();

    // Rcpp::Rcout << "Line 18. rowind = " << rowind  << ". \n" ;

    arma::vec xvars = xrow(ancvars);

    arma::vec input_temp = (xvars - ancsplits)/tau;
    arma::vec psi_temp = 1/(1+ arma::exp(- input_temp));


    // Rcpp::Rcout << "Line 31. rowind = " << rowind  << ". \n" ;


    //test that this gives the right numbers
    arma::vec probvec = 1 - psi_temp;


    probvec.elem(leftinds) = psi_temp.elem(leftinds);


    // Obtain products for each terminal node

    // Rcpp::Rcout << "Line 45. rowind = " << rowind  << ". \n" ;



    arma::vec retvec(ancterms.n_elem);
    // Rcpp::Rcout << "Line 51. rowind = " << rowind  << ". \n" ;

    for(unsigned int termind=0; termind < ancterms.n_elem; termind++){
      double term_temp = ancterms(termind);

      arma::uvec tempinds = arma::find(anctermcol == term_temp);

      // Rcpp::Rcout << "Line 62. tempinds = " << tempinds  << ". \n" ;
      // Rcpp::Rcout << "Line 63. term_temp = " << term_temp  << ". \n" ;
      //
      // Rcpp::Rcout << "Line 65. termind = " << termind  << ". \n" ;

      arma::vec termprobs = probvec.elem(tempinds);

      double temp_prod = arma::prod(termprobs);

      retvec(termind) = temp_prod;
    }

    // Rcpp::Rcout << "Line 65. rowind = " << rowind  << ". \n" ;


    phimat.row(rowind) = retvec.t();

  }

  return( phimat);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat phi_app_softHS( arma::mat treemat,
                      arma::mat internalmat,
                      arma::mat xmat,
                      arma::mat probmat) {

  // int treenrow = treemat.n_rows;

  // arma::field<std::string> treemat(treemat1.nrow(), treemat1.ncol());
  //
  // arma::field<std::string> internalmat();

  // Rcpp::Rcout << "treemat = " << treemat << " .\n " ;
  // Rcpp::Rcout << "internalmat = " << internalmat << " .\n " ;

  //ensure that int is a matrix (maybe use as.matrix)

  int ntemp = xmat.n_rows;

  // Rcpp::Rcout << "internalmat.n_rows = " << internalmat.n_rows << " .\n " ;

  // Rcpp::Rcout << "ntemp = " << ntemp << " .\n " ;

  arma::mat phi_matrix(ntemp, internalmat.n_rows);

  for(unsigned int j=0; j < internalmat.n_rows ; j++){

    // Rcpp::Rcout << "Line 151, j = " <<  j << " .\n " ;

    int child_left = internalmat(j,2);
    int child_right = internalmat(j,3);

    // Rcpp::Rcout << "Line 156, j = " <<  j << " .\n " ;


    arma::uvec children_left = get_children_cpp(treemat, child_left);
    arma::uvec children_right = get_children_cpp(treemat, child_right);

    // Rcpp::Rcout << "Line 162, j = " <<  j << " .\n " ;

    double n_left = internalmat(j, 4);

    // Rcpp::Rcout << "Line 166, j = " <<  j << " .\n " ;

    double n_right = internalmat(j, 5);


    // Rcpp::Rcout << "Line 171, j = " <<  j << " .\n " ;

    // Rcpp::Rcout << "Line 47, j = " <<  j << " .\n " ;

    double tempprod = n_left*n_right;
    double temp_denom = std::sqrt(tempprod );

    for(int i=0; i < ntemp ; i++){

      // Rcpp::Rcout << "Line 174, i = " <<  i << " .\n " ;

      // Rcpp::Rcout << "Line 182, probmat = " <<  probmat << " .\n " ;
      arma::vec probvec = (probmat.row(i)).t();
      // Rcpp::Rcout << "Line 54, probvec = " <<  probvec << " .\n " ;

      // In original R code subtract 1 from indices
      // Subtract 2 here becasue C++ indexes from 1
      arma::vec tempprobs_left = probvec.elem(children_left - 2);
      arma::vec tempprobs_right = probvec.elem(children_right - 2);

      // Rcpp::Rcout << "Line 184, i = " <<  i << " .\n " ;

      double sumleft = arma::sum(tempprobs_left);
      double sumright = arma::sum(tempprobs_right);

      phi_matrix(i,j) = (n_right*sumleft - n_left*sumright )/temp_denom;


    }

  }


  return( phi_matrix);
}



