//
// Created by victorvianna on 12/12/18.
//
#include "transformation.h"
#include "signature.h"

Transformation::Transformation(Signature &a, Signature &b) {
    origin_index = a.get_point_index();
    image_index = b.get_point_index();
    // TODO: CHECK OVERFLOW
    s = (a.getKMin() / b.getKMin() + a.getKMax() / b.getKMax()) / 2;

    Eigen::Matrix3d A, B;
    A << a.getMinCurv(), a.getMaxCurv(), a.getNormal();
    B << b.getMinCurv(), b.getMaxCurv(), b.getNormal();

    R = B * A.transpose();
    t = b.getPointCoordinates() - s * R * a.getPointCoordinates();
}

Transformation::Transformation(std::vector<double> &point) {
    if (point.size() != 7 + 2)
        throw std::length_error("invalid point size");
    s = point[0];
    R = (Eigen::AngleAxisd(point[1], Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(point[2], Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(point[3], Eigen::Vector3d::UnitZ())).toRotationMatrix();
    t.resize(3, 1);
    t << point[4], point[5], point[6];
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
    auto t_vector = std::vector<double>(t.data(), t.data() + t.rows() * t.cols());
    Eigen::Vector3d ea = R.eulerAngles(0, 1 , 2);
    auto ea_vector = std::vector<double>(ea.data(), ea.data() + ea.rows() * ea.cols());
    return {s, ea_vector[0], ea_vector[1], ea_vector[2], t_vector[0], t_vector[1], t_vector[2],
            (double) origin_index, (double) image_index};
}

std::ostream& operator << (std::ostream& os, Transformation &t) {
    auto temp = t.to_point();
    os << std::vector<double>(temp.begin(), temp.end() - 2);
    return os;
}

Eigen::MatrixXd Transformation::apply(Eigen::MatrixXd& point){
    return (t + s * R * point.transpose()).transpose();
}