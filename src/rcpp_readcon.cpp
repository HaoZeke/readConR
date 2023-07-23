#include <fstream>
#include <string>

#include "include/BaseTypes.hpp"
#include "include/FormatConstants.hpp"
#include "include/Helpers.hpp"
#include "include/ReadCon.hpp"
#include "include/helpers/StringHelpers.hpp"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::DataFrame readCon(std::string filename) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(filename);
  yodecon::types::ConFrame tmp;
  yodecon::process_header(
      (fconts | std::views::take(yodecon::constants::HeaderLength)), tmp);
  yodecon::process_coordinates(fconts, tmp);

  // Create vectors for each field
  std::vector<std::string> symbolVector;
  std::vector<double> xVector, yVector, zVector;
  std::vector<bool> isFixedVector;
  std::vector<uint64_t> atomIdVector;

  // Iterate through the data and append it to the vectors
  for (const auto &atomDatum : tmp.atom_data) {
    symbolVector.push_back(atomDatum.symbol);
    xVector.push_back(atomDatum.x);
    yVector.push_back(atomDatum.y);
    zVector.push_back(atomDatum.z);
    isFixedVector.push_back(atomDatum.is_fixed);
    atomIdVector.push_back(atomDatum.atom_id);
  }

  // Create a DataFrame with the vectors
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
      Rcpp::_["symbol"] = symbolVector, Rcpp::_["x"] = xVector,
      Rcpp::_["y"] = yVector, Rcpp::_["z"] = zVector,
      Rcpp::_["is_fixed"] = isFixedVector, Rcpp::_["atom_id"] = atomIdVector,
      Rcpp::_["stringsAsFactors"] = true);

  return df;
}
