#include <fstream>
#include <string>

#include "include/ReadCon.hpp"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List readCon(std::string filename) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(filename);
  auto singleCon =
      yodecon::create_single_con<yodecon::types::ConFrameVec>(fconts);
  auto atmNumVector = yodecon::symbols_to_atomic_numbers(singleCon.symbol);
  // Convert atmNumVector to IntegerVector
  Rcpp::IntegerVector atmNumVectorRcpp = Rcpp::wrap(atmNumVector);

  // Create a DataFrame with the vectors
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
      Rcpp::_["symbol"] = singleCon.symbol,
      Rcpp::_["atmNum"] = atmNumVectorRcpp, Rcpp::_["x"] = singleCon.x,
      Rcpp::_["y"] = singleCon.y, Rcpp::_["z"] = singleCon.z,
      Rcpp::_["is_fixed"] = singleCon.is_fixed,
      Rcpp::_["atom_id"] = singleCon.atom_id,
      Rcpp::_["stringsAsFactors"] = false);

  Rcpp::List conFrame = Rcpp::List::create(
      Rcpp::Named("prebox_header") =
          yodecon::helpers::string::to_csv_string(singleCon.prebox_header),
      Rcpp::Named("boxl") = singleCon.boxl,
      Rcpp::Named("angles") = singleCon.angles,
      Rcpp::Named("postbox_header") =
          yodecon::helpers::string::to_csv_string(singleCon.postbox_header),
      Rcpp::Named("natm_types") = singleCon.natm_types,
      Rcpp::Named("natms_per_type") = singleCon.natms_per_type,
      Rcpp::Named("masses_per_type") = singleCon.masses_per_type,
      Rcpp::Named("atom_data") = df);
  conFrame.attr("class") = "conFrame";
  return conFrame;
}
