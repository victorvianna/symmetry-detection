//
// Created by victorvianna on 12/12/18.
//
#include "transformation.h"
#include "signature.h"

Transformation::Transformation(Signature &a, Signature &b) {
    origin_index = a.get_point_index();
    image_index = b.get_point_index();
    // TODO: CHECK OVERFLOW
    s = a.getKMin() / b.getKMin() + a.getKMax() / b.getKMax();

    Eigen::Matrix3d A, B;
    A << a.getMinCurv(), a.getMaxCurv(), a.getNormal();
    B << b.getMinCurv(), b.getMaxCurv(), b.getNormal();

    Eigen::Matrix3d rotation_matrix = B * A.transpose();
    Eigen::Vector3d euler_angles = rotation_matrix.eulerAngles(0, 1 , 2);
    R = std::vector<double>(euler_angles.data(), euler_angles.data() + euler_angles.rows() * euler_angles.cols());

    Eigen::MatrixXd translation = b.getPointCoordinates() - s * rotation_matrix * a.getPointCoordinates();
    t = std::vector<double>(translation.data(), translation.data() + translation.rows() * translation.cols());
}

Transformation::Transformation(std::vector<double> &point) {
    if (point.size() != 7 + 2)
        throw std::length_error("invalid point size");
    s = point[0];
    R = std::vector<double>(point.begin() + 1, point.begin() + 4);
    t = std::vector<double>(point.begin() + 4, point.begin() + 7);
    origin_index = point[7];
    image_index = point[8];
}

void Transformation::to_points(std::vector <Transformation> transf_space,
                                                            std::vector <std::vector<double>> &points) {
    points.clear();
    for (Transformation &t : transf_space)
        points.push_back(t.to_point());
}

std::vector<double> Transformation::to_point() {
    return {s, R[0], R[1], R[2], t[0], t[1], t[2], (double) origin_index, (double) image_index};
}

std::ostream& operator << (std::ostream& os, Transformation &t) {
    auto temp = t.to_point();
    os << std::vector<double>(temp.begin(), temp.end() - 2);
    return os;
}