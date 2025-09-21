#ifndef KDTREE_H
#define KDTREE_H

#include "alg_visual.h"
#include <vector>
#include <memory>

class KDTree {
public:
    explicit KDTree(statistics& stat);

    std::vector<int> rangeQuery(double qx_LB, double qx_UB,
                                double qy_LB, double qy_UB) const;


    int nearestNeighbor(const double* query) const;

private:
    struct KDNode {
        int pointIdx;      
        int axis;       
        std::unique_ptr<KDNode> left, right;
        KDNode(int idx, int a) : pointIdx(idx), axis(a) {}
    };

    statistics& stat;     
    std::unique_ptr<KDNode> root;


    std::unique_ptr<KDNode> buildTree(std::vector<int>& indices, int depth);
    void splitByMedian(std::vector<int>& indices, int axis, size_t medianPos);


    void rangeSearch(const KDNode* node,
                     double qx_LB, double qx_UB,
                     double qy_LB, double qy_UB,
                     std::vector<int>& result) const;

    void nearestSearch(const KDNode* node, const double* query,
                       int& bestIdx, double& minDistSq) const;
};


#endif // KDTREE_H
