//
// Created by gabriel on 11/12/18.
//

#include "signature.h"
#include <igl/principal_curvature.h>
#include <igl/avg_edge_length.h>
#include <iostream>

Signature::Signature(double _kMin, double _kMax, Eigen::MatrixXd _minCurv, Eigen::MatrixXd _maxCurv,
                     Eigen::MatrixXd _normal, int _point_index, Eigen::MatrixXd _pointCoordinates) :
        kMin(_kMin), kMax(_kMax), minCurv(_minCurv), maxCurv(_maxCurv), normal(_normal), point_index(_point_index),
        pointCoordinates(_pointCoordinates){}

void Signature::build_signatures(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<Signature> &signatures) {
    Eigen::MatrixXd minCurvDir, maxCurvDir, normalDir;
    Eigen::VectorXd minCurvVal, maxCurvVal;

    // Compute normal direction
    igl::per_vertex_normals(V, F, normalDir);

    // Compute curvature directions via quadric fitting
    igl::principal_curvature(V, F, minCurvDir, maxCurvDir, minCurvVal, maxCurvVal);

    for (int i = 0; i < V.rows(); i++) {
        // Compute the normal as the cross product between the two principal curvature directions
        Eigen::Vector3d d1, d2, n;

        d1(0) = minCurvDir(i, 0);
        d1(1) = minCurvDir(i, 1);
        d1(2) = minCurvDir(i, 2);

        d2(0) = maxCurvDir(i, 0);
        d2(1) = maxCurvDir(i, 1);
        d2(2) = maxCurvDir(i, 2);

        n = d1.cross(d2);
        n.normalize();

        // Fix the orientation if the normal is inwards
        if (normalDir.row(i) * n < 0) {
            n = -n;
            minCurvDir.row(i) = -minCurvDir.row(i);
        }

        Eigen::MatrixXd coords;
        coords = V.row(i);
        coords.transposeInPlace();

        // Create the signature
        signatures.push_back(Signature(minCurvVal(i), maxCurvVal(i), minCurvDir.row(i), maxCurvDir.row(i),
                n.transpose(), i, coords));
    }
}

const int &Signature::get_point_index() const {
    return point_index;
}

int Signature::dimension() {
    return 3 * 3 + 2;
}

void Signature::plot_directions(igl::opengl::glfw::Viewer &viewer, Eigen::MatrixXd &V, double length,
                                bool showMin, bool showMax, bool showNormal) {
    Eigen::RowVector3d red(0.8, 0.2, 0.2), green(0.2, 0.8, 0.2), blue(0.2, 0.2, 0.8);

    // Draw the three directions
    if (showMin)
        viewer.data().add_edges(V.row(point_index) - minCurv * length, V.row(point_index) + minCurv * length, blue);
    if (showMax)
        viewer.data().add_edges(V.row(point_index) - maxCurv * length, V.row(point_index) + maxCurv * length, red);
    if (showNormal)
        viewer.data().add_edges(V.row(point_index), V.row(point_index) + normal * length, green);
}

void Signature::plot_all_directions(igl::opengl::glfw::Viewer &viewer, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                                    std::vector<Signature> &signatures,
                                    bool showMin, bool showMax, bool showNormal) {
    double length = igl::avg_edge_length(V, F);

    for (Signature s : signatures)
        s.plot_directions(viewer, V, length, showMin, showMax, showNormal);
}

std::vector<double> Signature::flatten(bool rigid) {
    std::vector<double> flattened;

    if (rigid) {
        flattened.push_back(kMin);
        flattened.push_back(kMax);
    } else {
        flattened.push_back(kMin / kMax);
    }

    return flattened;
}

double *Signature::flatten(std::vector<Signature> &signatures, bool rigid) {
    auto *all_flattened = new double[Signature::dimension() * signatures.size()];
    std::vector<Signature>::iterator it;
    double *p;
    for (it = signatures.begin(), p = all_flattened; it != signatures.end(); it++) {
        auto flattened = it->flatten(rigid);
        copy(flattened.begin(), flattened.end(), p);
        p += flattened.size();
    }
    return all_flattened;
}

bool Signature::is_not_umbilical_point() {
    return fabs(kMin / kMax) < 0.75;
}

double Signature::getKMin() const {
    return kMin;
}

double Signature::getKMax() const {
    return kMax;
}

const Eigen::MatrixXd &Signature::getMinCurv() const {
    return minCurv;
}

const Eigen::MatrixXd &Signature::getMaxCurv() const {
    return maxCurv;
}

const Eigen::MatrixXd &Signature::getNormal() const {
    return normal;
}

const Eigen::MatrixXd &Signature::getPointCoordinates() const {
    return pointCoordinates;
}
