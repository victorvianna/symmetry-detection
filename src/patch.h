//
// Created by gabriel on 13/12/18.
//

#ifndef SYMMETRY_DETECTION_PATCH_H
#define SYMMETRY_DETECTION_PATCH_H

#include <Eigen/Geometry>
#include <unordered_set>
#include "transformation.h"

class PatchPair {
private:
    Transformation transform;
    std::vector<int> origin, image, oneRing;
    std::unordered_set<int> imageFrontier;
    double threshold;

    int closestFrontier(Eigen::MatrixXd point, Eigen::MatrixXd &V);
public:
    PatchPair(Transformation _transform, std::vector<std::vector<int>> &oneRing, double _threshold);

    bool insert(int newVertex, Eigen::MatrixXd &V, std::vector<std::vector<int>> &oneRing);

    int size();
};


#endif //SYMMETRY_DETECTION_PATCH_H
