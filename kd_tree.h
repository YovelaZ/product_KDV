#ifndef KDTREE_H
#define KDTREE_H

#include "alg_visual.h"
#include <vector>
#include <memory>

class KDTree {
public:
    explicit KDTree(statistics& stat);

    // 范围查询 [qx_LB,qx_UB]×[qy_LB,qy_UB]
    std::vector<int> rangeQuery(double qx_LB, double qx_UB,
                                double qy_LB, double qy_UB) const;

    // 最近邻查询（返回数据索引）
    int nearestNeighbor(const double* query) const;

private:
    struct KDNode {
        int pointIdx;      // 数据索引
        int axis;          // 分割轴 (0=x, 1=y)
        std::unique_ptr<KDNode> left, right;
        KDNode(int idx, int a) : pointIdx(idx), axis(a) {}
    };

    statistics& stat;      // 引用外部数据
    std::unique_ptr<KDNode> root;

    // 构建相关（改进分割策略）
    std::unique_ptr<KDNode> buildTree(std::vector<int>& indices, int depth);
    void splitByMedian(std::vector<int>& indices, int axis, size_t medianPos);

    // 查询相关
    void rangeSearch(const KDNode* node,
                     double qx_LB, double qx_UB,
                     double qy_LB, double qy_UB,
                     std::vector<int>& result) const;

    void nearestSearch(const KDNode* node, const double* query,
                       int& bestIdx, double& minDistSq) const;
};

#endif // KDTREE_H