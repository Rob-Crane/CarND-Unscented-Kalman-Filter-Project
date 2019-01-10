#include "tools.h"

using Eigen::VectorXd;
using std::vector;

VectorXd Tools::calculateRMSE(const vector<VectorXd> &estimations,
                       const vector<VectorXd> &ground_truth) {
  assert (estimations.size() == ground_truth.size());
  using v_sz = vector<VectorXd>::size_type;
  v_sz n = estimations.size();
  
  VectorXd ret(4);
  ret << 0, 0, 0, 0;
  for (v_sz i = 0; i < n; ++i) {
    VectorXd residual = ground_truth[i] - estimations[i];
    ret += VectorXd(residual.array() * residual.array());
  }

  ret /= n;
  ret = ret.array().sqrt();
  return ret;
}
