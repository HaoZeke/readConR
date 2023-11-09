#include <fstream>
#include <string>

#include "include/ReadCon.hpp"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List readCon(std::string filename) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(filename);
  auto tmp = yodecon::create_single_con<yodecon::types::ConFrameVec>(fconts);
  auto atmNumVector = yodecon::symbols_to_atomic_numbers(tmp.symbol);
  // Convert atmNumVector to IntegerVector
  Rcpp::IntegerVector atmNumVectorRcpp = Rcpp::wrap(atmNumVector);

  // Create a DataFrame with the vectors
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
      Rcpp::_["symbol"] = tmp.symbol, Rcpp::_["atmNum"] = atmNumVectorRcpp,
      Rcpp::_["x"] = tmp.x, Rcpp::_["y"] = tmp.y, Rcpp::_["z"] = tmp.z,
      Rcpp::_["is_fixed"] = tmp.is_fixed, Rcpp::_["atom_id"] = tmp.atom_id,
      Rcpp::_["stringsAsFactors"] = false);

  Rcpp::List conFrame = Rcpp::List::create(
      Rcpp::Named("prebox_header") =
          yodecon::helpers::string::to_csv_string(tmp.prebox_header),
      Rcpp::Named("boxl") = tmp.boxl, Rcpp::Named("angles") = tmp.angles,
      Rcpp::Named("postbox_header") =
          yodecon::helpers::string::to_csv_string(tmp.postbox_header),
      Rcpp::Named("natm_types") = tmp.natm_types,
      Rcpp::Named("natms_per_type") = tmp.natms_per_type,
      Rcpp::Named("masses_per_type") = tmp.masses_per_type,
      Rcpp::Named("atom_data") = df);
  conFrame.attr("class") = "conFrame";
  return conFrame;
}
