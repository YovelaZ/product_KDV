#include "range_tree.h"

using namespace std;
RangeTree::RangeTree(statistics& stat_ref) : stat(stat_ref), root(nullptr) {
    vector<int> indices(stat.n);
    for (int i = 0; i < stat.n; ++i) indices[i] = i;

    sort(indices.begin(), indices.end(), [&](int a, int b) {
        if (stat.featureVector[a][0] != stat.featureVector[b][0])
            return stat.featureVector[a][0] < stat.featureVector[b][0];
        if (stat.featureVector[a][1] != stat.featureVector[b][1])
            return stat.featureVector[a][1] < stat.featureVector[b][1];
        return a < b; // 完全相同时，按原始索引排序
    });

    root = build_x_tree(indices);
}

RangeTree::~RangeTree() {
    delete_tree(root);
}
void RangeTree::delete_tree(Node* node) {
    if (!node) return; // 空节点直接返回

    // 递归删除左子树
    delete_tree(node->left);

    // 递归删除右子树
    delete_tree(node->right);

    // 递归删除Y子树（如果存在）
    delete_tree(node->y_tree);

    // 释放当前节点内存
    delete node;
}
RangeTree::Node* RangeTree::build_x_tree(vector<int>& indices) {
    if (indices.empty()) return nullptr;

    Node* node = new Node();
    node->left = node->right = node->y_tree = nullptr;
    node->is_leaf = false;

    if (indices.size() == 1) {
        node->key_index = indices[0];
        node->is_leaf = true;
        vector<int> y_indices = indices;
        sort(y_indices.begin(), y_indices.end(), [&](int a, int b) {
            if (stat.featureVector[a][1] != stat.featureVector[b][1])
                return stat.featureVector[a][1] < stat.featureVector[b][1];
            return a < b;
        });

        node->y_tree = build_y_tree(y_indices);

        // 叶子节点不需要y树，直接存储值
        return node;
    }

    size_t mid = indices.size() / 2;
    node->key_index = indices[mid];

    vector<int> left_indices(indices.begin(), indices.begin() + mid);
    vector<int> right_indices(indices.begin() + mid + 1, indices.end());

    node->left = build_x_tree(left_indices);
    node->right = build_x_tree(right_indices);

    // 非叶子节点构建y树
    vector<int> y_indices = indices;
    sort(y_indices.begin(), y_indices.end(), [&](int a, int b) {
        if (stat.featureVector[a][1] != stat.featureVector[b][1])
            return stat.featureVector[a][1] < stat.featureVector[b][1];
        return a < b;
    });

    node->y_tree = build_y_tree(y_indices);

    return node;
}
RangeTree::Node* RangeTree::build_y_tree(vector<int>& indices) {
    if (indices.empty()) return nullptr;

    Node* node = new Node();
    node->left = node->right = nullptr;
    node->y_tree = nullptr;
    node->is_leaf = false;

    if (indices.size() == 1) {
        node->key_index = indices[0];
        node->is_leaf = true;
        return node;
    }

    size_t mid = indices.size() / 2;
    node->key_index = indices[mid];

    vector<int> left_indices(indices.begin(), indices.begin() + mid);
    vector<int> right_indices(indices.begin() + mid + 1, indices.end());

    node->left = build_y_tree(left_indices);
    node->right = build_y_tree(right_indices);

    return node;
}




void RangeTree::query(double x_low, double x_high, double y_low, double y_high) {
    stat.density = 0;
    query_tree(root, x_low, x_high, y_low, y_high);
}
RangeTree::Node* RangeTree::FINDSPLITNODE_X(Node* node, double low, double high) {
    if (!node) {return nullptr;}
    double x = stat.featureVector[node->key_index][0];
    while(!node->is_leaf && (high <= x ||low > x)){

        if(high<=x){
            if(node->left){
                node = node->left;
                x = stat.featureVector[node->key_index][0];
            }else{
                break;
            }
        } else {
            if(node->right){
                node = node->right;
                x = stat.featureVector[node->key_index][0];
            }else{
                break;
            }
        }
    }
    return node;
}
RangeTree::Node* RangeTree::FINDSPLITNODE_Y(Node* node, double low, double high) {
    if (!node) {return nullptr;}
    double y = stat.featureVector[node->key_index][1];
    while(!node->is_leaf && (high <= y ||low > y)){
        if(high<=y){
            if(node->left){
                node = node->left;
                y = stat.featureVector[node->key_index][1];
            }else{
                break;
            }

        } else {
            if(node->right){
                node = node->right;
                y = stat.featureVector[node->key_index][1];
            }else{
                break;
            }

        }
    }
    return node;
}


void RangeTree::query_tree(Node* node, double x_low, double x_high, double y_low, double y_high) {
    if (!node) return;
    Node* v_split = FINDSPLITNODE_X(node,x_low,x_high);
    double x = stat.featureVector[v_split->key_index][0];

    if (v_split->is_leaf) {
        // 处理叶子节点
        double y = stat.featureVector[v_split->key_index][1];
        if (x >= x_low && x <= x_high && y >= y_low && y <= y_high) {
            stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                            compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
        }
        return;
    }else{
        double y = stat.featureVector[v_split->key_index][1];
        if (x >= x_low && x <= x_high && y >= y_low && y <= y_high) {
            stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                            compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
        }
        if(v_split->left){
            Node* v_left = v_split->left;

            while(!v_left->is_leaf){
                x = stat.featureVector[v_left->key_index][0];
                y = stat.featureVector[v_left->key_index][1];
                if(x_low<=x){
                    if ( x>=x_low && x<=x_high &&y >= y_low && y <= y_high) {
                        stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                                        compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
                    }
                    if(v_left->right){
                        query_y_tree(v_left->right->y_tree, y_low,y_high);
                    }
                    if(v_left->left){
                        v_left = v_left->left;
                    } else{
                        break;
                    }

                } else {
                    if(v_left->right){
                        v_left = v_left->right;
                    }
                    else{
                        break;
                    }

                }
            }
            if(v_left->is_leaf){
                x = stat.featureVector[v_left->key_index][0];
                y = stat.featureVector[v_left->key_index][1];
                if ( x>=x_low && x<=x_high && y >= y_low && y <= y_high) {
                    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                                    compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
                }
            }

        }

        if(v_split->right){
            Node* v_right = v_split->right;

            while(!v_right->is_leaf){
                double x = stat.featureVector[v_right->key_index][0];
                double y = stat.featureVector[v_right->key_index][1];
                if(x_high>=x){
                    if ( x>=x_low && x<=x_high &&y >= y_low && y <= y_high) {
                        stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                                        compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
                    }
                    if(v_right->left){
                        query_y_tree(v_right->left->y_tree, y_low,y_high);
                    }

                    if(v_right->right){
                        v_right = v_right->right;
                    }else{
                        break;
                    }
                } else {
                    if(v_right->left){
                        v_right = v_right->left;
                    }else{
                        break;
                    }

                }
            }
            if(v_right->is_leaf){
                double x = stat.featureVector[v_right->key_index][0];
                double y = stat.featureVector[v_right->key_index][1];
                if ( x>=x_low && x<=x_high && y >= y_low && y <= y_high) {
                    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                                    compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
                }
            }
        }


    }


}
void RangeTree::traverse(Node* node) {
    if (!node) return; // 终止条件
    double x = stat.featureVector[node->key_index][0];
    double y = stat.featureVector[node->key_index][1];
    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                    compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
    // 递归遍历左子树
    traverse(node->left);
    // 递归遍历右子树
    traverse(node->right);
}
void RangeTree::query_y_tree(Node* node, double y_low, double y_high) {
    if (!node) return;
    Node* v_split = FINDSPLITNODE_Y(node,y_low,y_high);
    double y = stat.featureVector[v_split->key_index][1];
    if (v_split->is_leaf) {
        // 处理叶子节点
        double x = stat.featureVector[v_split->key_index][0];
        if ( y >= y_low && y <= y_high) {
            stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                            compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
        }
    }else{
        double x = stat.featureVector[v_split->key_index][0];
        if ( y >= y_low && y <= y_high) {
            stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                            compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
        }
        if( v_split->left){
            Node* v_left = v_split->left;

            while(!v_left->is_leaf){
                double x = stat.featureVector[v_left->key_index][0];
                double y = stat.featureVector[v_left->key_index][1];
                if(y_low<=y){

                    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                            compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);


                    if(v_left->right){
                        traverse(v_left->right);
                    }

                    if(v_left->left){
                        v_left = v_left->left;
                    } else{
                        break;
                    }

                }else{
                    if(v_left->right){
                        v_left = v_left->right;
                    } else{
                        break;
                    }

                }
            }
            if(v_left->is_leaf){
                double x = stat.featureVector[v_left->key_index][0];
                double y = stat.featureVector[v_left->key_index][1];
                if ( y >= y_low && y <= y_high) {
                    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                                    compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
                }
            }

        }

        if(v_split->right){
            Node* v_right = v_split->right;

            while(!v_right->is_leaf){
                double x = stat.featureVector[v_right->key_index][0];
                double y = stat.featureVector[v_right->key_index][1];
                if(y_high>=y){
                    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                            compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);

                    if(v_right->left){
                        traverse(v_right->left);
                    }

                    if(v_right->right){
                        v_right = v_right->right;
                    }else{
                        break;
                    }

                }else{
                    if(v_right->left){
                        v_right = v_right->left;
                    }else{
                        break;
                    }

                }
            }
            if(v_right->is_leaf){
                double x = stat.featureVector[v_right->key_index][0];
                double y = stat.featureVector[v_right->key_index][1];
                if ( y >= y_low && y <= y_high) {
                    stat.density += compute_density(stat.q[0], x, stat.k_type_x, stat.b_x) *
                                    compute_density(stat.q[1], y, stat.k_type_y, stat.b_y);
                }
            }

        }


    }
    return;
}