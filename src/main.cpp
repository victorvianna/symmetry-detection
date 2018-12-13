#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/principal_curvature.h>
#include "flann/flann.hpp"
#include "mean_shift.h"
#include "kernel.h"
#include "transformation.h"

#include "signature.h"

using namespace std;

vector<Signature> prune_points(vector<Signature> &signatures) {
    vector<Signature> pruned;

    for (Signature &s : signatures) {
        if (s.is_not_umbilical_point())
            pruned.push_back(s);
    }

    return pruned;
}


template<typename T>
vector<T> random_sample(vector<T> &v, int num_samples) {
    std::vector<int> indices(v.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    vector<T> sampled;
    for (int idx : indices)
        sampled.push_back(v[idx]);
    return sampled;
};

void build_pairing(vector<Signature> &signatures, vector<Transformation> &transf, bool rigid,
    const int NUM_SAMPLES = std::min<int>((int) signatures.size(), 100);
    vector<Signature> samples = random_sample(signatures, NUM_SAMPLES);
    flann::Matrix<double> datapoints(Signature::flatten(signatures, rigid), signatures.size(), Signature::dimension(rigid));
    flann::Matrix<double> query(Signature::flatten(samples, rigid), samples.size(), Signature::dimension(rigid));
    flann::Index<flann::L2<double>> index(datapoints, flann::KDTreeIndexParams(4));
    index.buildIndex();
    double radius = 0.6;
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<double>> dists;
    index.radiusSearch(query, indices, dists, radius, flann::SearchParams(128));
    for (int i = 0; i < samples.size(); i++) {
        vector<int> &neighbors = indices[i];
        Signature &p_a = samples[i];
        for (int j : neighbors) {
            Signature &p_b = signatures[j];
            transf.push_back(Transformation(p_a, p_b));
            transf.push_back(Transformation(p_b, p_a));
        }
    }
}


void run_clustering(vector<Transformation>& transf_space, vector<vector<Transformation>>& clusters_transf) {
    double kernel_bandwidth = 3;
    double beta_1 = 1, beta_2 = 1, beta_3 = 1;
    vector<double> weights = {beta_1, beta_2, beta_2, beta_2, beta_3, beta_3, beta_3, 0, 0};
    MeanShift *msp = new MeanShift(epanechnikov_kernel, weights);

    vector<vector<double> > points = Transformation::to_points(transf_space);
    vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);
    for (Cluster &cluster : clusters) {
        vector<Transformation> cluster_transf;
        for (vector<double> &point : cluster.original_points)
            cluster_transf.push_back(Transformation(point));
        clusters_transf.push_back(cluster_transf);
    }

}

int main() {

    bool rigid = false;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("mesh/bunny.off", V, F);

    // Plot the mesh
    //igl::opengl::glfw::Viewer viewer;
    //viewer.data().set_mesh(V, F);

    // signature computation
    vector<Signature> signatures;
    Signature::build_signatures(V, F, signatures);

    // pruning
    signatures = prune_points(signatures);

    // Plot after pruning
    //Signature::plot_all_directions(viewer, V, F, signatures);

    // pairing
    vector<Transformation> transf_space;
    build_pairing(signatures, transf_space, rigid);

    // clustering
    vector<vector<Transformation>> clusters;
    run_clustering(transf_space, clusters);

    // run verification and build patches
    /// TODO

    // display patches or write to file
    /// TODO

    //viewer.launch();

    return 0;
}
