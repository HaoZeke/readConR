#include <fstream>
#include <string>

#include "include/ReadCon.hpp"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List readCon(std::string filename) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(filename);
  yodecon::types::ConFrame tmp;
  yodecon::process_header(
      (fconts | ranges::views::take(yodecon::constants::HeaderLength)), tmp);
  yodecon::process_coordinates(fconts, tmp);

  // Create vectors for each field
  std::vector<std::string> symbolVector;
  std::vector<double> xVector, yVector, zVector;
  std::vector<bool> isFixedVector;
  std::vector<int> atomIdVector;

  // Iterate through the data and append it to the vectors
  for (const auto &atomDatum : tmp.atom_data) {
    symbolVector.push_back(atomDatum.symbol);
    xVector.push_back(atomDatum.x);
    yVector.push_back(atomDatum.y);
    zVector.push_back(atomDatum.z);
    isFixedVector.push_back(atomDatum.is_fixed);
    atomIdVector.push_back(atomDatum.atom_id);
  }

  auto atmNumVector = yodecon::symbols_to_atomic_numbers(symbolVector);
  // Convert atmNumVector to IntegerVector
  Rcpp::IntegerVector atmNumVectorRcpp = Rcpp::wrap(atmNumVector);

  // Create a DataFrame with the vectors
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
      Rcpp::_["symbol"] = symbolVector, Rcpp::_["atmNum"] = atmNumVectorRcpp,
      Rcpp::_["x"] = xVector, Rcpp::_["y"] = yVector, Rcpp::_["z"] = zVector,
      Rcpp::_["is_fixed"] = isFixedVector, Rcpp::_["atom_id"] = atomIdVector,
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
