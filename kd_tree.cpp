#include "kd_tree.h"
#include <algorithm>
#include <cmath>
#include <limits>

KDTree::KDTree(statistics& stat) : stat(stat) {
    std::vector<int> indices(stat.n);
    for (int i = 0; i < stat.n; ++i) indices[i] = i;
    root = buildTree(indices, 0);
}

std::unique_ptr<KDTree::KDNode> KDTree::buildTree(std::vector<int>& indices, int depth) {
    if (indices.empty()) return nullptr;

    const int axis = depth % stat.dim;
    const size_t medianPos = indices.size() / 2;


    splitByMedian(indices, axis, medianPos);

    auto node = std::make_unique<KDNode>(indices[medianPos], axis);


    std::vector<int> leftIndices(indices.begin(), indices.begin() + medianPos);
    std::vector<int> rightIndices(indices.begin() + medianPos + 1, indices.end());

    node->left = buildTree(leftIndices, depth + 1);
    node->right = buildTree(rightIndices, depth + 1);

    return node;
}

void KDTree::splitByMedian(std::vector<int>& indices, int axis, size_t k) {

    std::nth_element(indices.begin(), indices.begin() + k, indices.end(),
                     [this, axis](int a, int b) {
                         return stat.featureVector[a][axis] < stat.featureVector[b][axis];
                     });
}

std::vector<int> KDTree::rangeQuery(double qx_LB, double qx_UB,
                                    double qy_LB, double qy_UB) const {
    std::vector<int> result;
    rangeSearch(root.get(), qx_LB, qx_UB, qy_LB, qy_UB, result);
    return result;
}

void KDTree::rangeSearch(const KDNode* node,
                         double qx_LB, double qx_UB,
                         double qy_LB, double qy_UB,
                         std::vector<int>& result) const {
    if (!node) return;

    const double x = stat.featureVector[node->pointIdx][0];
    const double y = stat.featureVector[node->pointIdx][1];


    if (x >= qx_LB && x <= qx_UB && y >= qy_LB && y <= qy_UB) {
//        result.push_back(node->pointIdx);
        stat.density = stat.density + compute_density(stat.q[0],stat.featureVector[node->pointIdx][0],stat.k_type_x,stat.b_x)*compute_density(stat.q[1],stat.featureVector[node->pointIdx][1],stat.k_type_y,stat.b_y);
    }


    if (node->axis == 0) { 
        if (qx_LB <= x) rangeSearch(node->left.get(), qx_LB, qx_UB, qy_LB, qy_UB, result);
        if (qx_UB >= x) rangeSearch(node->right.get(), qx_LB, qx_UB, qy_LB, qy_UB, result);
    } else { 
        if (qy_LB <= y) rangeSearch(node->left.get(), qx_LB, qx_UB, qy_LB, qy_UB, result);
        if (qy_UB >= y) rangeSearch(node->right.get(), qx_LB, qx_UB, qy_LB, qy_UB, result);
    }
}

int KDTree::nearestNeighbor(const double* query) const {
    int bestIdx = -1;
    double minDistSq = std::numeric_limits<double>::max();
    nearestSearch(root.get(), query, bestIdx, minDistSq);
    return bestIdx;
}

void KDTree::nearestSearch(const KDNode* node, const double* query,
                           int& bestIdx, double& minDistSq) const {
    if (!node) return;

    const double dx = stat.featureVector[node->pointIdx][0] - query[0];
    const double dy = stat.featureVector[node->pointIdx][1] - query[1];
    const double distSq = dx * dx + dy * dy;

    if (distSq < minDistSq) {
        minDistSq = distSq;
        bestIdx = node->pointIdx;
    }

    const int axis = node->axis;
    const double nodeValue = stat.featureVector[node->pointIdx][axis];
    const double queryValue = query[axis];


    KDNode* first = (queryValue < nodeValue) ? node->left.get() : node->right.get();
    KDNode* second = (queryValue < nodeValue) ? node->right.get() : node->left.get();

    if (first) nearestSearch(first, query, bestIdx, minDistSq);


    if (second && std::pow(queryValue - nodeValue, 2) < minDistSq) {
        nearestSearch(second, query, bestIdx, minDistSq);
    }

}
