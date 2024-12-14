#include "mesh_simplification.h"

void buildAdjacencyList(std::vector<std::set<int>>& adjacency_list, const std::vector<veci3>& faces, int numVertices) {
    adjacency_list.resize(numVertices);
    for (const auto& face : faces) {
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                adjacency_list[face[i]].insert(face[j]);
                adjacency_list[face[j]].insert(face[i]);
            }
        }
    }
}
struct Compare {
    bool operator()(struct Edge A, struct Edge B) {
        return B < A;
    }
};
Model* simplify_mesh(
    const std::vector<vecf3>& _vertices,    // positions of vertices in the mesh
    const std::vector<veci3>& _faces,       // indices of vertices in each face
    float ratio                             // the ratio of the number of vertices after simplification to the original number of vertices
) {

    // avoid modifying the original mesh
    std::vector<vecf3> vertices = _vertices;
    std::vector<veci3> faces = _faces;
    std::vector<std::set<int>> adjacency_list;
    // record whether the vertex is deleted
    std::deque<bool> vertices_deleted(vertices.size(), false);

    // record the face index of each vertex,
    std::vector<std::vector<int>> faces_of_vertices(vertices.size());

    for (int i = 0; i < faces.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            faces_of_vertices[faces[i][j]].push_back(i);
        }
    }
    for (int i = 0; i < vertices.size(); ++i) {
        if (faces_of_vertices[i].empty()) {
            vertices_deleted[i] = true;
        }
    }
    //std::cout << faces.size() << std::endl;;

    // TODO 3.1:
    // compute the Q matrices for all the initial vertices
    //
    std::vector < Eigen::Matrix4f> Q_faces(faces.size());
    for (int i = 0;i < faces.size();++i) {
        vecf3 x = vertices[faces[i][0]] - vertices[faces[i][1]];
        vecf3 y = vertices[faces[i][0]] - vertices[faces[i][2]];
        vecf3 n = x.cross(y);
        n.normalize();
        double a = n[0];
        double b = n[1];
        double c = n[2];
        double d = -(a * vertices[faces[i][0]][0] + b * vertices[faces[i][0]][1] + c * vertices[faces[i][0]][2]);
        Q_faces[i] << a * a, a* b, a* c, a* d,
            a* b, b* b, b* c, b* d,
            a* c, b* c, c* c, d* c,
            a* d, b* d, c* d, d* d;
    }
    //
    std::vector < Eigen::Matrix4f> Q_vertices(vertices.size(), Eigen::Matrix4f::Zero());
    for (int i = 0;i < vertices.size();++i) {
        for (int j = 0;j < faces_of_vertices[i].size();++j) {
            Q_vertices[i] += Q_faces[faces_of_vertices[i][j]];
        }
    }


    // TODO 3.2:
    // select all valid pairs(edges) and compute the cost of each edge
    //std::cout << vertices.size();
    buildAdjacencyList(adjacency_list, faces, vertices.size());
    std::priority_queue<Edge, std::vector<struct Edge>, Compare> edges;
    // 
    //std::cout << adjacency_list.size();
    for (int vertex = 0; vertex < adjacency_list.size(); ++vertex) {
        
        for (int adjacentVertex : adjacency_list[vertex]) {
            //
            if (vertex < adjacentVertex) {
                //       (vertex, adjacentVertex)
                struct Edge E;
                E.first = vertex;
                E.second = adjacentVertex;

                Eigen::Matrix4f Q_merge;
                Q_merge = Q_vertices[vertex] + Q_vertices[adjacentVertex];
                Eigen::Matrix4f Q_trans;
                Q_trans.block<3, 4>(0, 0) = Q_merge.block<3, 4>(0, 0);
                Q_trans(3,0) = 0;
                Q_trans(3, 1) = 0;
                Q_trans(3, 2) = 0;
                Q_trans(3, 3) = 1;
                vecf4 X;
                X[0] = 0;
                X[1] = 0;
                X[2] = 0;
                X[3] = 1;
                if (Q_trans.determinant() != 0) {
                    E.new_V = Q_trans.inverse() * X;
                    E.inver = 1;
                }
                else {
                    E.inver = 0;
                    E.new_V[0] = vertices[adjacentVertex][0];
                    E.new_V[1] = vertices[adjacentVertex][1];
                    E.new_V[2] = vertices[adjacentVertex][2];
                    E.new_V[3] = 1;

                }
                    E.cost = E.new_V.transpose() * Q_merge * E.new_V;
                edges.push(E);
                //std::cout << edges.empty() << std::endl;

            }
        }
    }

    //std::cout << "complete 3.2" << std::endl;
    // TODO 3.3:
    // iteratively remove the pair of the least cost from the heap
    uint32_t face_cnt = faces.size();
    uint32_t target_face_cnt = face_cnt * ratio;
    int count = 0;
    while (face_cnt > target_face_cnt) {

        // remove the min edge from the heap
        if (edges.empty())break;
        //count++;
        //std::cout << count << std::endl;
        struct Edge E = edges.top();
        edges.pop();
        //std::cout << E.first << " " << E.second << std::endl;
        if (vertices_deleted[E.first] || vertices_deleted[E.second])continue;

        vertices_deleted[E.first] = true;

        //更新second坐标
        if (E.inver == 1) {
            vertices[E.second][0] = E.new_V[0];
            vertices[E.second][1] = E.new_V[1];
            vertices[E.second][2] = E.new_V[2];

        }

        for (int i = 0;i < faces_of_vertices[E.first].size();++i) {
            int one_face = faces_of_vertices[E.first][i];
            //                    
            for (int j = 0;j < 3;j++) {
                if (faces[one_face][j] == E.first) {
                    faces[one_face][j] = E.second;
                    break;
                }
            }
            // 
            if (faces[one_face][0] == faces[one_face][1] || faces[one_face][1] == faces[one_face][2] || faces[one_face][2] == faces[one_face][0]) {
                faces[one_face][0] = -1;
                face_cnt--;
            }
            //      ڽӱ  cost
            for (int x : adjacency_list[E.first]) {
                if (x != E.second) {
                    adjacency_list[x].insert(E.second);
                    adjacency_list[E.second].insert(x);
                    int vertex;
                    int adjacentVertex;
                    if (x < E.second) {
                        vertex = x;
                        adjacentVertex = E.second;
                    }
                    else {
                        vertex = E.second;
                        adjacentVertex = x;
                    }
                    struct Edge Ed;
                    Ed.first = vertex;
                    Ed.second = adjacentVertex;

                    Eigen::Matrix4f Q_merge = Q_vertices[vertex] + Q_vertices[adjacentVertex];
                    Q_vertices[adjacentVertex] = Q_merge;
                    Eigen::Matrix4f Q_trans;
                    Q_trans.block<3, 4>(0, 0) = Q_merge.block<3, 4>(0, 0);
                    Q_trans(3, 0) = 0;
                    Q_trans(3, 1) = 0;
                    Q_trans(3, 2) = 0;
                    Q_trans(3, 3) = 1;
                    vecf4 X;
                    X[0] = 0;
                    X[1] = 0;
                    X[2] = 0;
                    X[3] = 1;
                    if (Q_trans.determinant() != 0) {
                        Ed.new_V = Q_trans.inverse() * X;
                        Ed.inver = 1;
                    }
                    else {
                        Ed.inver = 0;
                        Ed.new_V[0] = vertices[adjacentVertex][0];
                        Ed.new_V[1] = vertices[adjacentVertex][1];
                        Ed.new_V[2] = vertices[adjacentVertex][2];
                        Ed.new_V[3] = 1;
                    }
                    Ed.cost = Ed.new_V.transpose() * Q_merge * Ed.new_V;
                    edges.push(Ed);
                }
            }
            //    faces_of_vertices    
            for (auto face : faces_of_vertices[E.first]) {
                assert(face != -1);
                faces_of_vertices[E.second].push_back(face);
            }

        }
        // update the costs of all valid pairs


        // maintain the faces
        // set face invalid (with -1, -1, -1)


    }
    //std::cout << "complete while";

    // create the new mesh
    int new_vert_cnt = 0;
    int new_face_cnt = 0;
    for (auto i = 0; i < vertices.size(); ++i) {
        if (!vertices_deleted[i]) {
            vertices[new_vert_cnt] = vertices[i];
            for (auto face : faces_of_vertices[i]) {
                //assert(face != -1);
                for (int j = 0; j < 3; ++j) {
                    if (faces[face][j] == i) {
                        faces[face][j] = new_vert_cnt;
                    }
                }
            }
            new_vert_cnt += 1;
        }
    }
    //std::cout << faces.size() << std::endl;;
    for (int i = 0; i < faces.size(); ++i) {
        if (faces[i][0] != -1) {
            faces[new_face_cnt] = faces[i];
            //std::cout << faces[i][0] <<" "<< faces[i][1] <<" "<< faces[i][2] << std::endl;
            new_face_cnt += 1;
        }
    }
    //std::cout << new_vert_cnt<<std::endl;
    //std::cout << new_face_cnt<<std::endl;
    vertices.resize(new_vert_cnt);
    faces.resize(new_face_cnt);
    return Model::load(vertices, faces);
}
