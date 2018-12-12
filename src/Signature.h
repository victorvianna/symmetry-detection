//
// Created by gabriel on 11/12/18.
//

#ifndef SYMMETRY_DETECTION_SIGNATURE_H
#define SYMMETRY_DETECTION_SIGNATURE_H

#include <Eigen/Geometry>
#include <igl/opengl/glfw/Viewer.h>

class Signature {
private:
    double kMin, kMax;
    Eigen::MatrixXd minCurv, maxCurv, normal;
    int point_index;
    Signature(double _kMin, double _kMax, Eigen::MatrixXd _minCurv, Eigen::MatrixXd _maxCurv, Eigen::MatrixXd _normal, int _point_index);

public:
    static void build_signatures(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<Signature> &signatures);

    int get_point_index();

    static int dimension();

    void plot_directions(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, double length,
                         bool showMin = true, bool showMax = true, bool showNormal = true);

    static void plot_all_directions(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<Signature> &signatures,
                                    bool showMin = true, bool showMax = true, bool showNormal = true);

    std::vector<double> flatten();

    static double* flatten(std::vector<Signature>& signatures);

    bool is_not_umbilic_point();
};


#endif //SYMMETRY_DETECTION_SIGNATURE_H
