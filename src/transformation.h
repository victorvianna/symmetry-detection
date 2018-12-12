//
// Created by victorvianna on 12/12/18.
//

#ifndef SYMMETRY_DETECTION_TRANSFORMATION_H
#define SYMMETRY_DETECTION_TRANSFORMATION_H

#include "signature.h"

class Transformation {
public:
    Transformation(Signature &a, Signature &b);

    explicit Transformation(std::vector<double> &point);

    static std::vector <std::vector<double>> to_points(std::vector <Transformation> transf_space);

    std::vector<double> to_point();

private:
    int origin_index, image_index;
    double s;
    std::vector<double> R, t;
};

#endif //SYMMETRY_DETECTION_TRANSFORMATION_H