//
// Created by gabriel on 13/12/18.
//

#ifndef SYMMETRY_DETECTION_VERIFICATION_H
#define SYMMETRY_DETECTION_VERIFICATION_H

#include <Eigen/Geometry>
#include "patch.h"
#include "transformation.h"

class Verification {
private:
    double threshold;
    std::vector<std::vector<int>> oneRing;
public:
    Verification(Eigen::MatrixXi& F, double threshold);

    void verifyCluster(Eigen::MatrixXd& V, std::vector<Transformation> &cluster, std::vector<PatchPair> &patches);
};


#endif //SYMMETRY_DETECTION_VERIFICATION_H
