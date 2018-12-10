#include <iostream>
#include <igl/readOFF.h>
#include "flann/flann.hpp"
#include "mean_shift.h"
#include "kernel.h"

using namespace std;

class Signature {
private:
    double k1, k2;
    vector<double> d1, d2, d3;
    int point_index;
public:
    Signature(double _k1, double _k2, vector<double> _d1, vector<double> _d2, vector<double> _d3, int _point_index) :
            k1(_k1), k2(_k2), d1(_d1), d2(_d2), d3(_d3), point_index(_point_index) {}

    static int dimension() {
        return 3 * 3 + 2;
    }

    int get_point_index() {
        return point_index;
    }

    static void build_signatures(Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<Signature> &signatures) {
        ///TODO
        for (int i = 0; i < V.rows(); i++)
            signatures.push_back(Signature(1, 1, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, i));
    }

    vector<double> flatten() {
        vector<double> flattened;
        flattened.push_back(k1);
        flattened.push_back(k2);
        flattened.insert(flattened.end(), d1.begin(), d1.end());
        flattened.insert(flattened.end(), d2.begin(), d2.end());
        flattened.insert(flattened.end(), d3.begin(), d3.end());
        return flattened;
    }

    static double *flatten(vector<Signature> &signatures) {
        double *all_flattened = new double[Signature::dimension() * signatures.size()];
        vector<Signature>::iterator it;
        double *p;
        for (it = signatures.begin(), p = all_flattened; it != signatures.end(); it++) {
            auto flattened = it->flatten();
            copy(flattened.begin(), flattened.end(), p);
            p += flattened.size();
        }
        return all_flattened;
    }
};

vector<Signature> prune_points(vector<Signature> &signatures) {
    /// TODO
    vector<Signature> pruned(signatures.begin(), signatures.end());
    return pruned;
}

class Transformation {
public:
    Transformation(Signature &a, Signature &b) {
        // TODO
        origin_index = a.get_point_index();
        image_index = b.get_point_index();
        s = 1;
        R = {0, 0, 0};
        t = {0, 0, 0};
    }

    Transformation(vector<double> &point) {
        if (point.size() != 7 + 2)
            throw string("Invalid point size");
        s = point[0];
        R = vector<double>(point.begin() + 1, point.begin() + 4);
        t = vector<double>(point.begin() + 4, point.begin() + 7);
        origin_index = point[7];
        image_index = point[8];
    }

    static vector<vector<double>> to_points(vector<Transformation> transf_space) {
        vector<vector<double>> points;
        for (Transformation &t : transf_space) {
            points.push_back(t.to_point());
        }
        return points;
    }

    vector<double> to_point() {
        return {s, R[0], R[1], R[2], t[0], t[1], t[2], (double) origin_index, (double) image_index};
    }

private:
    int origin_index, image_index;
    double s;
    vector<double> R, t;
};

template<typename T>
vector<T> random_sample(vector<T> &v) {
    /// TODO
    vector<T> sampled(v.begin(), v.end());
    return sampled;
};

void build_pairing(vector<Signature> &signatures, vector<Transformation> &transf) {
    vector<Signature> samples = random_sample(signatures);
    flann::Matrix<double> datapoints(Signature::flatten(signatures), signatures.size(), Signature::dimension());
    flann::Matrix<double> query(Signature::flatten(samples), samples.size(), Signature::dimension());
    flann::Index<flann::L2<double>> index(datapoints, flann::KDTreeIndexParams(4));
    index.buildIndex();
    double radius = 0.6;
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<double>> dists;
    index.radiusSearch(query, indices, dists, radius, flann::SearchParams(128));
    for (int i = 0; i < samples.size(); i++) {
        vector<int> &neighbors = indices[i];
        Signature &p_i = samples[i];
        for (int j : neighbors) {
            Signature &p_j = samples[j];
            transf.push_back(Transformation(p_i, p_j));
            transf.push_back(Transformation(p_j, p_i));
        }
    }
}


void run_clustering(vector<Transformation> &transf_space, vector<vector<Transformation>> &clusters_transf) {
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

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOFF("mesh/tetrahedron.off", V, F);

    // signature computation
    vector<Signature> signatures;
    Signature::build_signatures(V, F, signatures);

    // pruning
    signatures = prune_points(signatures);

    // pairing
    vector<Transformation> transf_space;
    build_pairing(signatures, transf_space);

    // clustering
    vector<vector<Transformation>> clusters;
    run_clustering(transf_space, clusters);


    // run verification and build patches
    /// TODO

    // display patches or write to file
    /// TODO
    return 0;
}