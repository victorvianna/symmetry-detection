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
    int point_index;

    Signature(double _kMin, double _kMax, Eigen::MatrixXd _minCurv, Eigen::MatrixXd _maxCurv,
              Eigen::MatrixXd _normal, int _point_index, Eigen::MatrixXd _pointCoordinates);

    Eigen::MatrixXd minCurv, maxCurv, normal, pointCoordinates;

public:

    double getKMin() const;

    double getKMax() const;

    const Eigen::MatrixXd &getMinCurv() const;

    const Eigen::MatrixXd &getMaxCurv() const;

    const Eigen::MatrixXd &getNormal() const;

    const Eigen::MatrixXd &getPointCoordinates() const;

    static void build_signatures(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<Signature> &signatures);

    const int &get_point_index() const;

    static int dimension();

    void plot_directions(igl::opengl::glfw::Viewer &viewer, Eigen::MatrixXd &V, double length,
                         bool showMin = true, bool showMax = true, bool showNormal = true);

    static void plot_all_directions(igl::opengl::glfw::Viewer &viewer, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                                    std::vector<Signature> &signatures,
                                    bool showMin = true, bool showMax = true, bool showNormal = true);

    std::vector<double> flatten(bool rigid);

    static double *flatten(std::vector<Signature> &signatures, bool rigid);

    bool is_not_umbilical_point();
};


#endif //SYMMETRY_DETECTION_SIGNATURE_H
