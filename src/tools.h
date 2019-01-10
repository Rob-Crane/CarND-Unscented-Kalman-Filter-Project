#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

namespace Tools {

/**
* A helper method to calculate RMSE.
*/
Eigen::VectorXd calculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                              const std::vector<Eigen::VectorXd> &ground_truth);

} // namespace Tools

#endif  // TOOLS_H_
