#pragma once
#ifndef RANGE_TREE_H
#define RANGE_TREE_H

#include "alg_visual.h"
#include <vector>
#include <algorithm>

class RangeTree {
public:
    RangeTree(statistics& stat);
    ~RangeTree();

    void query(double x_low, double x_high, double y_low, double y_high);

private:
    struct Node {
        int key_index;          // 关键点索引
        Node* left;
        Node* right;
        Node* y_tree;           // y轴子树
        bool is_leaf;           // 标记是否为叶子节点
    };

    Node* root;
    statistics& stat;


    Node* build_x_tree(std::vector<int>& indices);
    Node* build_y_tree(std::vector<int>& indices);
    Node* FINDSPLITNODE_X(Node* node, double low, double high);
    Node* FINDSPLITNODE_Y(Node* node, double low, double high);
    void query_tree(Node* node, double x_low, double x_high, double y_low, double y_high);

//    void query_x_tree(Node* node, double x_low, double x_high, double y_low, double y_high);
    void traverse(Node* node);
    void query_y_tree(Node* node, double y_low, double y_high);
    void delete_tree(Node* node);
};

#endif // RANGE_TREE_H