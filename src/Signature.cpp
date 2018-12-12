//
// Created by gabriel on 11/12/18.
//

#include "Signature.h"
#include <igl/principal_curvature.h>
#include <igl/avg_edge_length.h>
#include <iostream>

Signature::Signature(double _kMin, double _kMax, Eigen::MatrixXd _minCurv, Eigen::MatrixXd _maxCurv, Eigen::MatrixXd _normal, int _point_index) :
    kMin(_kMin), kMax(_kMax), minCurv(_minCurv), maxCurv(_maxCurv), normal(_normal), point_index(_point_index) {}

void Signature::build_signatures(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<Signature> &signatures) {
    Eigen::MatrixXd minCurvDir, maxCurvDir, normalDir;
    Eigen::VectorXd minCurvVal, maxCurvVal;

    // Compute normal direction
    igl::per_vertex_normals(V, F, normalDir);

    // Compute curvature directions via quadric fitting
    igl::principal_curvature(V, F, minCurvDir, maxCurvDir, minCurvVal, maxCurvVal);

    for (int i = 0; i < V.rows(); i++)
        signatures.push_back(Signature(minCurvVal(i), maxCurvVal(i), minCurvDir.row(i), maxCurvDir.row(i), normalDir.row(i), i));
}

int Signature::get_point_index() {
    return point_index;
}

int Signature::dimension() {
    return 3 * 3 + 2;
}

void Signature::plot_directions(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, double length){
    Eigen::RowVector3d red(0.8, 0.2, 0.2), green(0.2, 0.8, 0.2), blue(0.2, 0.2, 0.8);

    // Draw the three directions
    viewer.data().add_edges(V.row(point_index) + minCurv * length, V.row(point_index) - minCurv * length, blue);
    viewer.data().add_edges(V.row(point_index) + maxCurv * length, V.row(point_index) - maxCurv * length, red);
    viewer.data().add_edges(V.row(point_index) + normal * length, V.row(point_index) - normal * length, green);
}

void Signature::plot_all_directions(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::vector<Signature>& signatures){
    double length = igl::avg_edge_length(V, F);

    for(Signature s : signatures)
        s.plot_directions(viewer, V, length);
}

std::vector<double> Signature::flatten() {
    std::vector<double> flattened;
    flattened.push_back(kMin);
    flattened.push_back(kMax);
    flattened.insert(flattened.end(), minCurv.data(), minCurv.data() + minCurv.cols());
    flattened.insert(flattened.end(), maxCurv.data(), maxCurv.data() + maxCurv.cols());
    flattened.insert(flattened.end(), normal.data(), normal.data() + normal.cols());
    return flattened;
}

double* Signature::flatten(std::vector<Signature>& signatures) {
    double *all_flattened = new double[Signature::dimension() * signatures.size()];
    std::vector<Signature>::iterator it;
    double *p;
    for (it = signatures.begin(), p = all_flattened; it != signatures.end(); it++) {
        auto flattened = it->flatten();
        copy(flattened.begin(), flattened.end(), p);
        p += flattened.size();
    }
    return all_flattened;
}

bool Signature::not_umbilic_point() {
    return fabs(kMin / kMax) < 0.75;
}