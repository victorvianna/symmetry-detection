//
// Created by victorvianna on 12/12/18.
//

#ifndef SYMMETRY_DETECTION_TRANSFORMATION_H
#define SYMMETRY_DETECTION_TRANSFORMATION_H

#include "signature.h"
#include "io.h"

class Transformation {
public:
    Transformation(Signature &a, Signature &b, bool rigid, bool reflection);
    explicit Transformation(std::vector<double> &point, bool reflection);

    static void to_points(std::vector <Transformation> transf_space,
                                                       std::vector <std::vector<double>> &points);

    std::vector<double> to_point();

    /**
    * @param point Point as a 1x3 matrix
    * @return Image of the transformation as 1x3 matrix
    */
    Eigen::MatrixXd apply(Eigen::MatrixXd &point);

    std::pair<int, int> getVertices();

    int getOrigin();

    int getImage();

private:
  int origin_index, image_index;
  bool has_reflection;
  double s;
  Eigen::Matrix3d R;
    Eigen::MatrixXd t;
};

std::ostream& operator << (std::ostream& os, Transformation &t);

#endif //SYMMETRY_DETECTION_TRANSFORMATION_H
