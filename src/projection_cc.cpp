#include <gismo.h>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <queue>

#include "knncpp.h"

#include <gsGSplines/gsGSplines.h>

#include <limits>

#include <cstdlib>  // for rand
#include <ctime>    // for rand

const float infinity = std::numeric_limits<float>::max();

#define PI 3.14159265358979323846


using namespace gismo;



typedef gsSurfMesh::Point Point;
typedef gsSurfMesh::Point Vector3d;
typedef gsSurfMesh::Vertex Vertex;
typedef gsSurfMesh::Edge Edge;


typedef gsEigen::MatrixXd Matrix;
typedef knncpp::Matrixi Matrixi;



float randFloat(float X){
    // Source - https://stackoverflow.com/a/686373
    // Posted by John Dibling, modified by community. See post 'Timeline' for change history
    // Retrieved 2026-05-29, License - CC BY-SA 3.0

    return static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/X));

}



double find_mean_edge_length (gsSurfMesh & mesh){
    double mean = 0;
    double length;
    for (auto e : mesh.edges()){
        length = (mesh.position(mesh.vertex(e,0)) - mesh.position(mesh.vertex(e, 1))).norm();
        mean += length / mesh.n_edges();
    }
    return mean;
}



void store_neighbours (gsSurfMesh & mesh){

    gsSurfMesh::Vertex_property <std::vector <Vertex> > neighbours = mesh.get_vertex_property<std::vector <Vertex> >("v:neighbours");
    if (neighbours) return;


    neighbours = mesh.add_vertex_property<std::vector <Vertex> >("v:neighbours");
    // initialize as empty
    for (auto v : mesh.vertices()){
        neighbours[v] = {};
    }

    for (auto v : mesh.vertices()){
        auto he_start = mesh.halfedge(v);
        neighbours[v].push_back(mesh.to_vertex(he_start));
        auto he = mesh.ccw_rotated_halfedge(he_start);
        while (he != he_start){
            neighbours[v].push_back(mesh.to_vertex(he));
            he = mesh.ccw_rotated_halfedge(he);
        }
    }
    
}


inline bool is_ev (gsSurfMesh & mesh, Vertex vert){
    return mesh.valence(vert) != 4;
}



bool has_connected_evs(gsSurfMesh & mesh){

    // load the neighbours of the verties
    gsSurfMesh::Vertex_property <std::vector <Vertex> > neighbours = mesh.get_vertex_property<std::vector <Vertex> >("v:neighbours");
    if (! neighbours){
        store_neighbours(mesh);
        neighbours = mesh.get_vertex_property<std::vector <Vertex> >("v:neighbours");
    }

    for (auto v : mesh.vertices()){
        if (mesh.valence(v) != 4){
            // it's ev
            for (auto n : neighbours[v]){
                if (mesh.valence(n) != 4){
                    // connected evs
                    return true;
                }
            }
        }
    }
    return false;
}




void mark_sharp_edges(gsSurfMesh& mesh, int angle_thr){

    bool  has_sharp = false;
    gsSurfMesh::Halfedge_property<bool> sharp = mesh.get_halfedge_property<bool>("h:sharp");
    if (sharp) has_sharp = true;

    // if sharp information exists, keep those 
    if (!has_sharp){ 

        gsSurfMesh::Halfedge_property<bool> sharp = mesh.add_halfedge_property<bool>("h:sharp");
        gsSurfMesh::Halfedge he;

        // iterate all edges
        for (gsSurfMesh::Edge_iterator it = mesh.edges_begin();
            it != mesh.edges_end();
            ++it)
        {
            auto e = *it;
            if (!mesh.is_boundary(e)){
                auto f1 = mesh.face(e, 0);
                auto f2 = mesh.face(e, 1);
                auto n1 = mesh.compute_face_normal(f1);
                auto n2 = mesh.compute_face_normal(f2);
                //auto v1 = mesh.vertex(e, 0);
                //auto v2 = mesh.vertex(e, 1);
                // dihedral angle
                float angle = std::acos(n1.dot(n2)) * 180 / PI;

                if (angle > angle_thr){
                    he = mesh.halfedge(e, 0);
                    sharp[he] = true;
                    he = mesh.halfedge(e, 1);
                    sharp[he] = true;
                }
            }
        }
    }
}







struct GridPoint{
    Point coords;
    int quad_idx = -1;
    double u;
    double v;
    bool is_sharp = false;
    gsSurfMesh::Edge sharp_edge;
};

struct FitPoint{
    Point coords;
    double u;
    double v;
};


inline float norm (Vector3d a){
    return std::sqrt(a.dot(a));
}


inline float angle (Vector3d a, Vector3d b){
    return std::acos(a.dot(b) / (norm(a) * norm(b))) * 180 / PI;
}



inline float calculate_distance(Point a, Point b){
    return std::sqrt( std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + std::pow(a[2] - b[2], 2));
}



bool project_point_on_triangle(Point p, Point t1, Point t2, Point t3, Point & result){
    Vector3d normal = (t3 - t1).cross(t2 - t1);
    Vector3d t1_to_point = p - t1;
    Vector3d triangle_to_point_orth = ( t1_to_point.dot(normal) ) / ( normal.dot(normal) ) * normal;
    result = p - triangle_to_point_orth;

    // Barycentric coordinates of the projection p' of p onto the triangle:
    // γ=(u×w)⋅n/‖n‖²
    float gamma = (t3-t1).cross(t1_to_point).dot(normal) / normal.dot(normal);
    // β=(w×v)⋅n/‖n‖²
    float beta = t1_to_point.cross(t2-t1).dot(normal) / normal.dot(normal);
    float alpha = 1 - gamma - beta;
    // The point p' lies inside T if:
    return ((0 <= alpha) && (alpha <= 1) &&
            (0 <= beta)  && (beta  <= 1) &&
            (0 <= gamma) && (gamma <= 1));
}


bool project_point_on_edge(Point p, Point e1, Point e2, Point & result){
    Vector3d pointer = p - e1;
    Vector3d e = e2 - e1;
    Vector3d pointer_proj = (pointer.dot(e)) / (e.dot(e)) * e;
    result = pointer_proj + e1;
    Vector3d diff = e - pointer_proj;
    // the point is on the edge if
    return (diff.squaredNorm() < e.squaredNorm() && diff.dot(e) > 0);
}




std::vector <gsSurfMesh::Face> triangle_1neighbourhood (gsSurfMesh & mesh, Vertex cv){
    // returns the triangles that vertex cv is part of

    if (mesh.is_isolated(cv)) return {};

    std::vector <gsSurfMesh::Face> neighbourhood = {};

    if (mesh.is_boundary(cv)){
        
        // find the boundary halfedge that points to cv
        auto he = mesh.halfedge(cv);
        he = mesh.next_halfedge(mesh.next_halfedge(he));
        while (! mesh.touches_boundary(he)){
            he = mesh.opposite_halfedge(mesh.next_halfedge(mesh.next_halfedge(he)));
        }
        // found ..

        neighbourhood.push_back(mesh.face(he));
        he = mesh.next_halfedge(he);
        while (! mesh.touches_boundary(he)){
            he = mesh.next_halfedge(mesh.opposite_halfedge(he));
            neighbourhood.push_back(mesh.face(he));
        }

        return neighbourhood;
    }

    
    auto he_start = mesh.halfedge(cv);
    neighbourhood.push_back(mesh.face(he_start));

    auto he = mesh.opposite_halfedge(mesh.next_halfedge(mesh.next_halfedge(he_start)));
    while (he != he_start){
        neighbourhood.push_back(mesh.face(he));
        he = mesh.opposite_halfedge(mesh.next_halfedge(mesh.next_halfedge(he)));
    }

    return neighbourhood;
}






void store_1neighbourhood (gsSurfMesh & mesh){
    // store for every vertex the triangles that it is part of

    if (! mesh.is_triangle_mesh())  GISMO_ERROR("\nT1-neighbourhood is defined only for triangle meshes!");

    gsSurfMesh::Vertex_property <std::vector <gsSurfMesh::Face> > T1neighbours = mesh.get_vertex_property<std::vector <gsSurfMesh::Face> >("v:T1neighbours");
    if (T1neighbours) return;


    T1neighbours = mesh.add_vertex_property<std::vector <gsSurfMesh::Face> >("v:T1neighbours");
    // initialize as empty
    for (auto v : mesh.vertices()){
        T1neighbours[v] = {};
    }

    // loop over all triangles and append them on their vertices' neighbour list
    for (auto triangle : mesh.faces()){
        auto he = mesh.halfedge(triangle);
        Vertex v1 = mesh.from_vertex(he);
        Vertex v2 = mesh.to_vertex(he);
        Vertex v3 = mesh.to_vertex(mesh.next_halfedge(he));
        T1neighbours[v1].push_back(triangle);
        T1neighbours[v2].push_back(triangle);
        T1neighbours[v3].push_back(triangle);
    }

}

void store_2neighbourhood (gsSurfMesh &mesh){
    // store on every vertex the triangles with vertices at most two edges away

    if (! mesh.is_triangle_mesh())  GISMO_ERROR("\nT2-neighbourhood is defined only for triangle meshes!");

    gsSurfMesh::Vertex_property <std::vector <gsSurfMesh::Face> > T2neighbours = mesh.get_vertex_property<std::vector <gsSurfMesh::Face> >("v:T2neighbours");
    if (T2neighbours) return;

    // load sharp tri vertices
    gsSurfMesh::Vertex_property <bool> sharp_tri_verts = mesh.get_vertex_property <bool> ("v:sharp_tri_verts");
    if (! sharp_tri_verts)
        GISMO_ERROR("\ninsert sharp vertex info before calling T2-neighbourhood\n");


    T2neighbours = mesh.add_vertex_property<std::vector <gsSurfMesh::Face> >("v:T2neighbours");
    // initialize as empty
    for (auto v : mesh.vertices()){
        T2neighbours[v] = {};
    }    

    // add property to acoount for already visited vertices
    gsSurfMesh::Vertex_property <bool> is_visited = mesh.add_vertex_property <bool> ("v:is_visited");
    for (auto v : mesh.vertices()){
        is_visited[v] = false;
    } 

    store_neighbours(mesh);
    // load the neighbours property
    gsSurfMesh::Vertex_property <std::vector <Vertex> > neighbours = mesh.get_vertex_property<std::vector <Vertex> >("v:neighbours");

    for (auto triangle : mesh.faces()){

        auto he = mesh.halfedge(triangle);
        Vertex v1 = mesh.from_vertex(he);
        Vertex v2 = mesh.to_vertex(he);
        Vertex v3 = mesh.to_vertex(mesh.next_halfedge(he));
        std::vector <Vertex> vs = {};

        // include only the non sharp verts in the vector vs
        if (!sharp_tri_verts[v1])
            vs.push_back(v1);
        if (!sharp_tri_verts[v2])
            vs.push_back(v2);
        if (!sharp_tri_verts[v3])
            vs.push_back(v3);

        for (auto v : vs){
            T2neighbours[v].push_back(triangle);

            for (auto n : neighbours[v]){

                if (is_visited[n]) continue;

                if (sharp_tri_verts[n]) continue;

                T2neighbours[n].push_back(triangle);
                is_visited[n] = true;
            }
        }

        // clear the is_visited property
        for (auto v : vs){
            for (auto n :neighbours[v]){
                is_visited[n] = false;
            }
        }
    } 
}








Point project_on_mesh(GridPoint gp, gsSurfMesh & mesh_quad_original, gsSurfMesh & mesh_quad, gsSurfMesh & mesh_tri, std::vector<Vertex> closest_verts, bool fix_degeneracies = true){


    Point p = gp.coords;
    Vertex closest_vert = closest_verts[0];

    std::vector <Point> possible_projections = {};
    Point projected;


    // if it is marked as sharp, project on its sharp edge....
    if (gp.is_sharp){

        // read map to sharp triangle edges
        gsSurfMesh::Edge_property <std::vector <gsSurfMesh::Edge> > sharps_map = mesh_quad_original.get_edge_property <std::vector <gsSurfMesh::Edge> > ("e:sharps_map");
        if (! sharps_map)
            GISMO_ERROR("\n No map to read\n");

        for (auto tri_edge : sharps_map[gp.sharp_edge]){
            Point e1 = mesh_tri.position( mesh_tri.vertex( tri_edge, 0 ) );
            Point e2 = mesh_tri.position( mesh_tri.vertex( tri_edge, 1 ) );

            bool is_inside = project_point_on_edge(p, e1, e2, projected);

            if (is_inside){
                possible_projections.push_back(projected);
            }
        }
    

        if (possible_projections.size() == 1){
            return possible_projections[0];
        }

        if (possible_projections.size() > 1){
            // return the closest one

            Point closest_proj;
            double min_distance = std::numeric_limits<double>::max();
            double distance;
            for (Point proj : possible_projections){
                distance = (proj - p).squaredNorm();
                if (distance < min_distance){
                    min_distance = distance;
                    closest_proj = proj;
                }
            }

            return closest_proj;
        }


        // else no projection on edge is possible
        // return the closest vertex that is on the sharp edge

        for (int current_cv_id = 0; current_cv_id < closest_verts.size(); current_cv_id++){
            closest_vert = closest_verts[current_cv_id];

            for (Edge sharp_tri : sharps_map[gp.sharp_edge]){
                Vertex v1 = mesh_tri.vertex( sharp_tri, 0 );
                Vertex v2 = mesh_tri.vertex( sharp_tri, 1 );

                if (closest_vert == v1 || closest_vert == v2)
                    return mesh_tri.position(closest_vert);
            }
            
        }
        
        // Not enough knn points to project on the actual sharp edge (though it should work with enough points)
        return mesh_tri.position(closest_verts[0]);

        //GISMO_ERROR("\nNOT ENOUGH KNN POINTS, REQUEST MORE!\n");
        

    }


    // it is not marked sharp...


    // if it is part of a sharp quad, check the direction of the closest vertex 
    // if it is on the wrong quad, fix it
    // we know it's in the wrong quad when the vector pointing to the projection of it in the sharp edge has a big angle
    // with the diretion to the sharp edge on the quad 
    
    // load the properties
    gsSurfMesh::Face_property<bool>     has_one_sharp      = mesh_quad.get_face_property<bool>    ("f:has_one_sharp");  
    gsSurfMesh::Face_property<Vector3d> direction_to_sharp = mesh_quad.get_face_property<Vector3d>("f:direction_to_sharp");
    gsSurfMesh::Face_property<gsSurfMesh::Halfedge> f_sharp_he = mesh_quad.get_face_property<gsSurfMesh::Halfedge> ("f:sharp_he");
    if (! has_one_sharp)
        GISMO_ERROR("\ndidn't find property...\n");
    if (! direction_to_sharp)
        GISMO_ERROR("\ndidn't find property...\n");
    if (! f_sharp_he)
        GISMO_ERROR("\ndidn't find property...\n");

    if (gp.quad_idx == -1)
        GISMO_ERROR("\nNon valid quad index.\n");
    auto gp_face = gsSurfMesh::Face(gp.quad_idx);

    int current_cv_id = 0;
    int vector_size = closest_verts.size();

    
    if (has_one_sharp[gp_face] && fix_degeneracies){
        // gsInfo << "found one...\n";

        // load sharp tri vertices
        gsSurfMesh::Vertex_property <bool> sharp_tri_verts = mesh_tri.get_vertex_property <bool> ("v:sharp_tri_verts");
        if (! sharp_tri_verts)
            GISMO_ERROR("\ninsert sharp vertex info before calling project_on_mesh\n");
        
        Vector3d dir = direction_to_sharp[gp_face];
        auto sharp_he = f_sharp_he[gp_face];
        Point from = mesh_quad.position( mesh_quad.from_vertex(sharp_he) );
        Point to   = mesh_quad.position( mesh_quad.to_vertex(sharp_he) );

        Point cvproj;

        while (true){
            Point cvpoint = mesh_tri.position(closest_vert);
            project_point_on_edge(cvpoint, from, to, cvproj);
            Vector3d pointer = cvproj - cvpoint;

            float ang = angle(pointer, dir);
            if ( (ang > 30 && ang < 330)  ||  sharp_tri_verts[closest_vert]){
                // discard this closest vert, get the next one
                // gsInfo << "skipped one...\n";
                // gsInfo << "------\n";
                // gsInfo << pointer << "\n" << dir  << "\n" << angle(pointer, dir) << "\n";
                if (current_cv_id < vector_size - 2){
                    current_cv_id++;
                    closest_vert = closest_verts[current_cv_id];
                }
                else{
                    GISMO_ERROR("\nNOT ENOUGH KNN POINTS, REQUEST MORE!\n");
                }
            }
            else{
                break;
            }
        }
    }




    // now closest vert is on the correct quad...



    // get adgjacent triangles of this vertex of distance 2
    gsSurfMesh::Vertex_property <std::vector <gsSurfMesh::Face> > T2neighbours = mesh_tri.get_vertex_property<std::vector <gsSurfMesh::Face> >("v:T2neighbours");
    if (! T2neighbours)
        GISMO_ERROR ("\nStore the T2-neighbours before calling project_on_mesh !\n");


    // project on each one triangle checking if it is possible
    

    for (auto triangle : T2neighbours[closest_vert]){   // project only on the neighbouring triangles 

        auto he  = mesh_tri.halfedge(triangle);
        Point t1 = mesh_tri.position(mesh_tri.from_vertex(he));
        Point t2 = mesh_tri.position(mesh_tri.to_vertex(he));
        Point t3 = mesh_tri.position(mesh_tri.to_vertex(mesh_tri.next_halfedge(he)));

        
        bool is_inside = project_point_on_triangle(p, t1, t2, t3, projected);

        if (is_inside){
            possible_projections.push_back(projected);
        }
    }


    // if only one possible projection, return it
    if (possible_projections.size() == 1){
        return possible_projections[0];
    }


    // if multiple projections are possible keep the closest one
    if (possible_projections.size() > 1){
        Point closest_proj;
        double min_distance = std::numeric_limits<double>::max();
        double distance;
        for (Point proj : possible_projections){
            distance = (proj - p).squaredNorm();
            if (distance < min_distance){
                min_distance = distance;
                closest_proj = proj;
            }
        }
        
        return closest_proj;
    }


    // if no projection is possible project on edges 


    // project on every edge of the t2-neighbourhood
    for (auto triangle : T2neighbours[closest_vert]){

        auto he    = mesh_tri.halfedge(triangle);
        for (int i = 0; i < 3; i++){
            Point e1 = mesh_tri.position( mesh_tri.from_vertex( he ) );
            Point e2 = mesh_tri.position( mesh_tri.to_vertex( he ) );

            bool is_inside = project_point_on_edge(p, e1, e2, projected);

            if (is_inside){
                possible_projections.push_back(projected);
            }

            he = mesh_tri.next_halfedge(he);
        }
    }

    if (possible_projections.size() == 1){
        return possible_projections[0];
    }

    if (possible_projections.size() > 1){
        // return the closest one

        Point closest_proj;
        double min_distance = std::numeric_limits<double>::max();
        double distance;
        for (Point proj : possible_projections){
            distance = (proj - p).squaredNorm();
            if (distance < min_distance){
                min_distance = distance;
                closest_proj = proj;
            }
        }
        
        return closest_proj;
    }
    

    // if it is not possible also on edges, return the closest vertex 

    return mesh_tri.position(closest_vert);
}







gsSurfMesh::Halfedge find_halfedge(gsSurfMesh & mesh, Vertex start, Vertex end){
    auto he = mesh.halfedge(start);
    for (unsigned int i = 0; i < mesh.valence(start); i++){
        if (mesh.to_vertex(he) == end)
            return he;
        else 
            he = mesh.ccw_rotated_halfedge(he);
    }
    GISMO_ERROR ("\nthere is no halfedge connecting the provided inputs.\n");
}



gsSurfMesh::Halfedge find_halfedge(gsSurfMesh & mesh, Vertex start, gsSurfMesh::Face f){
    auto he = mesh.halfedge(start);
    for (unsigned int i = 0; i < mesh.valence(start); i++){
        if (mesh.face(he) == f)
            return he;
        else 
            he = mesh.ccw_rotated_halfedge(he);
    }
    GISMO_ERROR ("\nthere is no halfedge from this vert on this face.\n");
}







void inherit_og_edge_info (gsSurfMesh & mesh, gsSurfMesh & smesh){
    // smesh (subdivided mesh) inherits the original edge information from mesh

    // read original edge property
    gsSurfMesh::Edge_property <gsSurfMesh::Edge> og_edge = mesh.get_edge_property <gsSurfMesh::Edge> ("e:og_edge");
    if (! og_edge)
        GISMO_ERROR("\nThere is no original edge property to read.\n");

    // create original edge property on the subdivided mesh
    gsSurfMesh::Edge_property <gsSurfMesh::Edge> og_edge_sub = smesh.get_edge_property <gsSurfMesh::Edge> ("e:og_edge");
    if (! og_edge_sub)
        og_edge_sub = smesh.add_edge_property <gsSurfMesh::Edge> ("e:og_edge");


    // read sharp info 
    gsSurfMesh::Halfedge_property <bool> sharp = mesh.get_halfedge_property <bool> ("h:sharp");
    if (!sharp)
        GISMO_ERROR("\n there is no sharp info to read\n");
    
    
    for (auto edge : mesh.edges()){
        // we only care on the sharps
        if (sharp[mesh.halfedge(edge,0)] || sharp[mesh.halfedge(edge,1)]){
            og_edge_sub[ gsSurfMesh::Edge( edge.idx() ) ] = og_edge[edge];
            og_edge_sub[ gsSurfMesh::Edge( edge.idx() + mesh.n_edges() ) ] = og_edge[edge];
        }

    }
}









void update_GridPoints_on_face (std::vector <GridPoint> & gps, 
                                gsSurfMesh & mesh, 
                                gsSurfMesh::Face & f, 
                                gsSurfMesh & final_mesh, 
                                double u_start, double u_end, 
                                double v_start, double v_end,
                                Vertex uv_starting_vert,
                                int quad_index){

    /*
        v4-----v3         v           x---e3---x      new points:  ----p3----      new faces:  -----------
        |      |          |           |        |                   |   |    |                  | f3 | f2 |
        |      |          |           e4       e2                  p4--p5---p2                 |----|----|
        v1-----v2         ---> u      |        |                   |   |    |                  | f4 | f1 |
                                      x---e1---x                   ----p1----                  -----------

        for new vertices on edges (p1,p2,p3,p4): new vert idx = old.n_vertices + edge.idx
        for new vertices on faces (p5)         : new vert idx = old.n_vertices + old.n_edges + face.idx

        v1 is the start of u,v coordinates

    */ 

    // auto he = mesh.halfedge(f);


    // read the property of which grid points have been added already
    gsSurfMesh::Vertex_property <bool> gp_added = final_mesh.get_vertex_property<bool>("v:gp_added");
    if (! gp_added)
        GISMO_ERROR("\ncreate a grid point has been added (v:gp_added) property before calling update_GridPoints_on_face.\n");



    // read the sharp info 
    gsSurfMesh::Halfedge_property<bool> sharp = mesh.get_halfedge_property<bool>("h:sharp");
    if (! sharp)   
        GISMO_ERROR("\nNo sharp information to read!\n");


    // read the original edge info 
    gsSurfMesh::Edge_property <gsSurfMesh::Edge> og_edge = mesh.get_edge_property <gsSurfMesh::Edge> ("e:og_edge");
    if (! og_edge)
        GISMO_ERROR("\nThere is no original edge info to read.\n");


    auto he = find_halfedge(mesh, uv_starting_vert, f);


    auto e1 = mesh.edge(he);
    he = mesh.next_halfedge(he); 
    auto e2 = mesh.edge(he);
    he = mesh.next_halfedge(he);
    auto e3 = mesh.edge(he);
    auto e4 = mesh.edge(mesh.next_halfedge(he));

    GridPoint gp;
    gp.quad_idx = quad_index;


    Vertex vp1 = Vertex(mesh.n_vertices() + e1.idx());
    if (! gp_added[vp1]){
        Point p1 = final_mesh.position(vp1);
        gp.coords = p1;
        gp.u = u_start + (u_end - u_start) / 2; 
        gp.v = v_start;
        

        // check if the edge is sharp
        if (sharp[mesh.halfedge(e1,0)] || sharp[mesh.halfedge(e1,1)]){
            gp.is_sharp = true;
            gp.sharp_edge = og_edge[e1];
        }
        else{
            gp.is_sharp = false;
        }

        gps.push_back(gp);

        gp_added[vp1] = true;
    }

    Vertex vp2 = Vertex(mesh.n_vertices() + e2.idx());
    if (! gp_added[vp2]){
        Point p2 = final_mesh.position(vp2);
        gp.coords = p2;
        gp.u = u_end; 
        gp.v = v_start + (v_end - v_start) / 2;

        // check if the edge is sharp
        if (sharp[mesh.halfedge(e2,0)] || sharp[mesh.halfedge(e2,1)]){
            gp.is_sharp = true;
            gp.sharp_edge = og_edge[e2];
        }
        else{
            gp.is_sharp = false;
        }

        gps.push_back(gp);

        gp_added[vp2] = true;
    }

    Vertex vp3 = Vertex(mesh.n_vertices() + e3.idx());
    if (! gp_added[vp3]){
        Point p3 = final_mesh.position(vp3);
        gp.coords = p3;
        gp.u = u_start + (u_end - u_start) / 2; 
        gp.v = v_end;

        // check if the edge is sharp
        if (sharp[mesh.halfedge(e3,0)] || sharp[mesh.halfedge(e3,1)]){
            gp.is_sharp = true;
            gp.sharp_edge = og_edge[e3];
        }
        else{
            gp.is_sharp = false;
        }

        gps.push_back(gp);

        gp_added[vp3] = true;
    }

    Vertex vp4 = Vertex(mesh.n_vertices() + e4.idx());
    if (! gp_added[vp4]){
        Point p4 = final_mesh.position(vp4);
        gp.coords = p4;
        gp.u = u_start; 
        gp.v = v_start + (v_end - v_start) / 2;

        // check if the edge is sharp
        if (sharp[mesh.halfedge(e4,0)] || sharp[mesh.halfedge(e4,1)]){
            gp.is_sharp = true;
            gp.sharp_edge = og_edge[e4];
        }
        else{
            gp.is_sharp = false;
        }

        gps.push_back(gp);

        gp_added[vp4] = true;
    }

    Vertex vp5 = Vertex(mesh.n_vertices() + mesh.n_edges() + f.idx());

    if (! gp_added[vp5]){
        Point p5 = final_mesh.position(vp5);
        gp.coords = p5;
        gp.u = u_start + (u_end - u_start) / 2; 
        gp.v = v_start + (v_end - v_start) / 2;
        gp.is_sharp = false;
        gps.push_back(gp);

        gp_added[vp5] = true;
    }

}






int main(int argc, char** argv){



    std::string fn;
    bool smooth = false;
    bool fix_degeneracies = false;
    bool use_c0 = false;
    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("smooth", "Use the smoothed quad mesh or the original.", smooth);
    cmd.addSwitch("fd", "Fix edge degeneracies.", fix_degeneracies);
    cmd.addSwitch("c0", "Use C0 basis on sharps", use_c0);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    
    std::string fn_tri = fn + "_p0.obj";
    std::string fn_quad;
    if (smooth)
        fn_quad = fn + "_quadrangulation_smooth.obj";
    else
        fn_quad = fn + "_quadrangulation.obj";
    std::string sharps_fn = fn + "_sharps.txt";
    
    gsSurfMesh mesh_tri;
    mesh_tri.read(fn_tri);
    gsSurfMesh mesh_quad;
    mesh_quad.read(fn_quad);




    // for noise on grid points 
    // double mean_edge_length = find_mean_edge_length(mesh_quad);
    // double noise_range = mean_edge_length / 8;


    // read sharp vert pairs
    std::ifstream sharps_in (sharps_fn);


    size_t quad_sharps_count;
    sharps_in >> quad_sharps_count;
    std::vector<std::pair <size_t, size_t> > sharp_vert_pairs (quad_sharps_count);
    for (size_t i = 0; i < quad_sharps_count; i++){
        sharps_in >> sharp_vert_pairs[i].first;
        sharps_in >> sharp_vert_pairs[i].second;
    }

    // mark sharp features
    gsSurfMesh::Halfedge_property<bool> sharp = mesh_quad.get_halfedge_property<bool>("h:sharp");
    if (!sharp){

        sharp = mesh_quad.add_halfedge_property<bool>("h:sharp");

        for (auto vert_pair : sharp_vert_pairs){
            auto he = mesh_quad.find_halfedge(gsSurfMesh::Vertex(vert_pair.first), gsSurfMesh::Vertex(vert_pair.second));
            sharp[he] = true;
            he = mesh_quad.opposite_halfedge(he);
            sharp[he] = true;
        }
    }



    // map sharp quad edges to sharp tri edges
    gsSurfMesh::Edge_property <std::vector <gsSurfMesh::Edge> > sharps_map = mesh_quad.add_edge_property <std::vector <gsSurfMesh::Edge> > ("e:sharps_map");

    std::ifstream in(fn + "_quad_to_tri_sharp_correspondence.txt");
    // gpt
    if (!in.is_open())
    {
        std::cerr << "Cannot open sharp correspondence file\n";
    }
    else
    {
        std::string line;

        // skip header
        std::getline(in, line);

        while (std::getline(in, line))
        {
            if (line.empty()) continue;

            std::stringstream ss(line);

            int qv0, qv1;
            ss >> qv0 >> qv1;

            std::string dummy;
            ss >> dummy; // '|'

            int feature_group;
            ss >> feature_group;

            ss >> dummy; // '|'

            //---------------------------------------------------
            // FIND EDGE INDEX IN MESH
            //---------------------------------------------------
            int edge_index = -1;
            int tri_edge_index = -1;

            for (auto e : mesh_quad.edges())
            {
                int v0 = mesh_quad.vertex(e, 0).idx();
                int v1 = mesh_quad.vertex(e, 1).idx();

                if ((v0 == qv0 && v1 == qv1) ||
                    (v0 == qv1 && v1 == qv0))
                {
                    edge_index = e.idx();
                    break;
                }
            }

            if (edge_index == -1)
            {
                GISMO_ERROR("\n did not find quad edge\n");
                std::cerr << "Warning: quad edge not found: "
                        << qv0 << " " << qv1 << "\n";
                continue;
                
            }

            gsSurfMesh::Edge qEdge(edge_index);

            //---------------------------------------------------
            // READ ALL TRIANGLE EDGES AND STORE THEM
            //---------------------------------------------------
            int tv0, tv1;

            while (ss >> tv0 >> tv1)
            {
                // WARNING: interpreting (tv0,tv1) as triangle vertex indices
                // and converting them into an EDGE in tri mesh is NOT possible
                // unless you also search the tri mesh.

                // If you actually stored TRI EDGE INDICES before:
                // -> then treat tv0 as edge index directly

                // searching the tri mesh instead: slow!!!

                for (auto e : mesh_tri.edges()){
                    int v0 = mesh_tri.vertex(e,0).idx();
                    int v1 = mesh_tri.vertex(e,1).idx();


                    if ((v0 == tv0 && v1 == tv1) ||
                        (v0 == tv1 && v1 == tv0)){
                            tri_edge_index = e.idx();
                            
                            break;
                        }
                }

                if (tri_edge_index == -1){
                    GISMO_ERROR("\ndid not find tri edge.\n");
                }

                gsSurfMesh::Edge triEdge(tri_edge_index);
                sharps_map[qEdge].push_back(triEdge);

                // NOTE: tv1 is ignored in this format version
            }
        }

        in.close();
    }


    // store sharp halfedge info on the triangle mesh...
    mark_sharp_edges(mesh_tri, 35);
    gsSurfMesh::Halfedge_property <bool> sharp_tri_he = mesh_tri.get_halfedge_property <bool> ("h:sharp");



    // mark vertices as sharp in the triangle mesh
    gsSurfMesh::Vertex_property <bool> sharp_tri_verts = mesh_tri.add_vertex_property <bool> ("v:sharp_tri_verts");
    // initialize
    for (auto v : mesh_tri.vertices()){
        sharp_tri_verts[v] = false;
    }
    for (auto e : mesh_tri.edges()){
        auto he = mesh_tri.halfedge(e,0);
        if (sharp_tri_he[he]){
            auto v1 = mesh_tri.from_vertex(he);
            auto v2 = mesh_tri.to_vertex(he);
            sharp_tri_verts[v1] = true;
            sharp_tri_verts[v2] = true;
        }
    }


    // 2 - neighbourhood calculate once
    store_2neighbourhood(mesh_tri);




    // create the original edges property on the first quad mesh
    gsSurfMesh::Edge_property <gsSurfMesh::Edge> og_edge = mesh_quad.add_edge_property <gsSurfMesh::Edge> ("e:og_edge");
    // mark every edge as the original
    for (auto e : mesh_quad.edges()){
        og_edge[e] = e;
    }

    gsSurfMesh mesh_quad_original = mesh_quad;


    // connected evs are not allowed
    if (has_connected_evs(mesh_quad)){
        // do one cc
        gsSubdivScheme smesh(mesh_quad);
        smesh.options().setInt("scheme", 0);
        smesh.subdivide(1);
        inherit_og_edge_info(mesh_quad_original, mesh_quad);
        og_edge = mesh_quad.get_edge_property <gsSurfMesh::Edge> ("e:og_edge");

        gsInfo << "Subdivided once to avoid connected evs.\n";
    }




    // subdivide one on a copy
    gsSurfMesh meshs1 = mesh_quad;
    gsSubdivScheme smesh1(meshs1);
    smesh1.options().setInt("scheme", 0);//Catmull-Clark
    smesh1.subdivide(1);
    inherit_og_edge_info(mesh_quad, meshs1);
    

    // subdivide a second time on a copy
    gsSurfMesh meshs2 = meshs1;
    gsSubdivScheme smesh2(meshs2);
    smesh2.options().setInt("scheme", 0);
    smesh2.subdivide(1);
    inherit_og_edge_info(meshs1, meshs2);
    

    // third time
    gsSurfMesh meshs3 = meshs2;
    gsSubdivScheme smesh3(meshs3);
    smesh3.options().setInt("scheme", 0);
    smesh3.subdivide(1);
    inherit_og_edge_info(meshs2, meshs3);





    // add properties for every face of the mesh_quad
    // if it has exactly one sharp edge 
    // the vector pointing from the sharp edge to its opposite edge
    gsSurfMesh::Face_property<bool>     has_one_sharp      = mesh_quad.add_face_property<bool>    ("f:has_one_sharp");  
    gsSurfMesh::Face_property<Vector3d> direction_to_sharp = mesh_quad.add_face_property<Vector3d>("f:direction_to_sharp");
    gsSurfMesh::Face_property<gsSurfMesh::Halfedge> f_sharp_he = mesh_quad.add_face_property<gsSurfMesh::Halfedge> ("f:sharp_he");
    
    for (auto f : mesh_quad.faces()){
        auto he = mesh_quad.halfedge(f);
        int sharp_count = 0;
        gsSurfMesh::Halfedge sharp_he;
        for (int i = 0; i < 4; i++){
            // for all 4 halfedges
            if (sharp[he]){
                sharp_count++;
                sharp_he = he;
            }
            he = mesh_quad.next_halfedge(he);
        }

        if (sharp_count == 1){
            has_one_sharp[f] = true;
            gsSurfMesh::Halfedge opposite_sh_he = mesh_quad.next_halfedge(mesh_quad.next_halfedge(sharp_he));
            // to find the vector pointing from the opposite perpendicular to the sharp
            // we project a point of the opposite - say the from_vertex - to the sharp
            // and calculate their difference
            Point random_point_on_opposite = mesh_quad.position( mesh_quad.from_vertex( opposite_sh_he ) );
            Point from = mesh_quad.position( mesh_quad.from_vertex(sharp_he) );
            Point to   = mesh_quad.position( mesh_quad.to_vertex  (sharp_he) );
            Point projected_on_sh;
            project_point_on_edge(random_point_on_opposite, from, to, projected_on_sh);
            direction_to_sharp[f] = projected_on_sh - random_point_on_opposite;
            // gsInfo << "HERE: " << projected_on_sh << "\n" << from << "\n" << projected_on_sh - from << "\n";

            f_sharp_he[f] = sharp_he;
        }
        else{
            has_one_sharp[f] = false;
            direction_to_sharp[f] = {0,0,0};
        }
    }





    // initialize a kdtree with the triangle mesh vertices 
    Matrix dataPoints(3, mesh_tri.n_vertices());
    int counter = 0;
    for (auto tri_v : mesh_tri.vertices()){
        Point pos = mesh_tri.position(tri_v);
        dataPoints(0, counter) = pos[0];
        dataPoints(1, counter) = pos[1];
        dataPoints(2, counter) = pos[2];
        counter++;
    }
    knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(dataPoints);
    kdtree.build();



    // mark if a tri vert has been used already
    gsSurfMesh::Vertex_property<bool> is_used = mesh_tri.add_vertex_property<bool>("v:is_used");
    for (auto v : mesh_tri.vertices()){
        is_used[v] = false;
    }


    // face property for grid points inside it
    //gsSurfMesh::Face_property<std::vector<GridPoint>> points_in_face = mesh_quad.add_face_property<std::vector<GridPoint>>("f:points");
    std::vector <std::vector <FitPoint> > points_in_face (mesh_quad.n_faces());

    // for debug and plotting
    gsSurfMesh::Vertex_property<int> corrFace = mesh_tri.add_vertex_property<int>("v:corrFace");


    // matrixes for fitting
    int amount_of_points_per_face = 81; // for 2 catmull clarks
    gsVector <index_t> offset (mesh_quad.n_faces() + 1);
    gsMatrix <real_t> uvs;
    gsMatrix <real_t> points;
    uvs.resize(2, amount_of_points_per_face * mesh_quad.n_faces());
    points.resize(3, amount_of_points_per_face * mesh_quad.n_faces());
    int current_offset = 0;
    int point_counter = 0;





    std::vector <GridPoint> gps;
    GridPoint gp;


    // property for visited grid points 
    gsSurfMesh::Vertex_property <bool> gp_added = meshs3.add_vertex_property<bool>("v:gp_added");
    



    for (auto f : mesh_quad.faces()){

        // update offset matrix
        offset(f.idx()) = current_offset;



        gps.clear();
        // add the four vertices as grid points
        // the coords should be after the subdivisions

        // first find the four vertices and the four edges
        auto he = mesh_quad.halfedge(f);
        auto e1 = mesh_quad.edge(he);
        Vertex v1 = mesh_quad.from_vertex(he);

        Vertex v2 = mesh_quad.to_vertex(he);

        he = mesh_quad.next_halfedge(he);
        auto e2 = mesh_quad.edge(he);
        Vertex v3 = mesh_quad.to_vertex(he);

        he = mesh_quad.next_halfedge(he);
        auto e3 = mesh_quad.edge(he);
        Vertex v4 = mesh_quad.to_vertex(he);

        auto e4 = mesh_quad.edge(mesh_quad.next_halfedge(he));


        // create 4 grid points for the 4 vertices
        GridPoint gp1;
        GridPoint gp2;
        GridPoint gp3;
        GridPoint gp4;

        // initialise
        gp1.is_sharp = false;
        gp2.is_sharp = false;
        gp3.is_sharp = false;
        gp4.is_sharp = false;

        // gsInfo << e1.idx() << ", " << e2.idx() << ", " << e3.idx() << ", " << e4.idx() << "\n";
        // int a;
        // std::cin >> a; 


        // identify from the edges which vertices should be sharp 
        if (sharp[mesh_quad.halfedge(e1,0)] || sharp[mesh_quad.halfedge(e1,1)]){
            gp1.is_sharp = true;
            gp1.sharp_edge = og_edge[e1];
            gp2.is_sharp = true;
            gp2.sharp_edge = og_edge[e1];
        }

        if (sharp[mesh_quad.halfedge(e2,0)] || sharp[mesh_quad.halfedge(e2,1)]){
            gp2.is_sharp = true;
            gp2.sharp_edge = og_edge[e2];
            gp3.is_sharp = true;
            gp3.sharp_edge = og_edge[e2];
        }

        if (sharp[mesh_quad.halfedge(e3,0)] || sharp[mesh_quad.halfedge(e3,1)]){
            gp3.is_sharp = true;
            gp3.sharp_edge = og_edge[e3];
            gp4.is_sharp = true;
            gp4.sharp_edge = og_edge[e3];
        }

        if (sharp[mesh_quad.halfedge(e4,0)] || sharp[mesh_quad.halfedge(e4,1)]){
            gp4.is_sharp = true;
            gp4.sharp_edge = og_edge[e4];
            gp1.is_sharp = true;
            gp1.sharp_edge = og_edge[e4];
        }

        // add the coordinate, u, v and quad index info on the grid points
        gp1.coords = meshs3.position(v1);
        gp1.u = 0; gp1.v = 0;
        gp1.quad_idx = f.idx();
        gps.push_back(gp1);

        
        gp2.coords = meshs3.position(v2);
        gp2.u = 1; gp2.v = 0;
        gp2.quad_idx = f.idx();
        gps.push_back(gp2);

        
        gp3.coords = meshs3.position(v3);
        gp3.u = 1; gp3.v = 1;
        gp3.quad_idx = f.idx();
        gps.push_back(gp3);

        
        gp4.coords = meshs3.position(v4);
        gp4. u = 0; gp4.v = 1;
        gp4.quad_idx = f.idx();
        gps.push_back(gp4);

        

        Vertex vp1 = Vertex(mesh_quad.n_vertices() + e1.idx());  // vertex of point p1
        Vertex vp2 = Vertex(mesh_quad.n_vertices() + e2.idx());
        Vertex vp3 = Vertex(mesh_quad.n_vertices() + e3.idx());
        Vertex vp4 = Vertex(mesh_quad.n_vertices() + e4.idx());
        Vertex vp5 = Vertex(mesh_quad.n_vertices() + mesh_quad.n_edges() + f.idx());


        update_GridPoints_on_face(gps, mesh_quad, f, meshs3, 0, 1, 0, 1, v1, f.idx());
        

        // find sub faces of f
        // the creation of the new quads start from v2, so f1 is the quad v2->p2->p5->p1
        // thus the one that corresponds to the halfedge v2->p2
        gsSurfMesh::Halfedge he2 = find_halfedge(meshs1, v2, vp2);
        gsSurfMesh::Face     f1  = meshs1.face(he2);


        // 
        //     v4-------p3--he3--v3             v               (1,1)
        //     |        |        |              |
        //     he4  f3  |   f2   |              | 
        //     |        |        |              |
        //     p4-------p5-------p2             |
        //     |        |        |              |
        //     |   f4   |   f1   he2            |
        //     |        |        |              |
        //     v1--he1--p1-------v2             -----------------> u
        // 

        auto he3 = meshs1.next_halfedge(meshs1.next_halfedge(meshs1.opposite_halfedge(meshs1.next_halfedge(he2))));
        auto he4 = meshs1.next_halfedge(meshs1.next_halfedge(meshs1.opposite_halfedge(meshs1.next_halfedge(he3))));
        auto he1 = meshs1.next_halfedge(meshs1.next_halfedge(meshs1.opposite_halfedge(meshs1.next_halfedge(he4))));
        auto f2  = meshs1.face(he3);
        auto f3  = meshs1.face(he4);
        auto f4  = meshs1.face(he1);

        // update grid points with appropriate uvs
        update_GridPoints_on_face(gps, meshs1, f1, meshs3, 0.5,  1,    0,    0.5, vp1, f.idx());
        update_GridPoints_on_face(gps, meshs1, f2, meshs3, 0.5,  1,    0.5,  1,   vp5, f.idx());
        update_GridPoints_on_face(gps, meshs1, f3, meshs3, 0,    0.5,  0.5,  1,   vp4, f.idx());
        update_GridPoints_on_face(gps, meshs1, f4, meshs3, 0,    0.5,  0,    0.5, v1,  f.idx());



        // add the 4 subfaces on a vector
        std::vector <gsSurfMesh::Face> sfs = {f1, f2, f3, f4};
        std::vector <Vertex> starting_verts = {vp1, vp5, vp4, v1};
        int i = 0;
        for (auto sf : sfs){

            i++;
            double sub_u_start;
            double sub_u_end;
            double sub_v_start;
            double sub_v_end;
            switch(i){
                case 1:
                    sub_u_start = 0.5;
                    sub_u_end = 1;
                    sub_v_start = 0;
                    sub_v_end = 0.5;
                    break;
                case 2:
                    sub_u_start = 0.5;
                    sub_u_end = 1;
                    sub_v_start = 0.5;
                    sub_v_end = 1;
                    break;
                case 3:
                    sub_u_start = 0;
                    sub_u_end = 0.5;
                    sub_v_start = 0.5;
                    sub_v_end = 1;
                    break;
                case 4:
                    sub_u_start = 0;
                    sub_u_end = 0.5;
                    sub_v_start = 0;
                    sub_v_end = 0.5;
                    break;
                default:
                    GISMO_ERROR("\nMore than four faces!?\n");
            }



            // find the sub-verts again:

            // v1 should be the start of the uv directions!
            auto she = find_halfedge(meshs1, starting_verts[i-1], sf);

            auto se1 = meshs1.edge(she);
            Vertex sv1 = meshs1.from_vertex(she);

            Vertex sv2 = meshs1.to_vertex(she);

            she = meshs1.next_halfedge(she);
            auto se2 = meshs1.edge(she);
            Vertex sv3 = meshs1.to_vertex(she);

            she = meshs1.next_halfedge(she);
            auto se3 = meshs1.edge(she);
            Vertex sv4 = meshs1.to_vertex(she);

            auto se4 = meshs1.edge(meshs1.next_halfedge(she));

            Vertex svp1 = Vertex(meshs1.n_vertices() + se1.idx());  // vertex of point p1
            Vertex svp2 = Vertex(meshs1.n_vertices() + se2.idx());
            Vertex svp3 = Vertex(meshs1.n_vertices() + se3.idx());
            Vertex svp4 = Vertex(meshs1.n_vertices() + se4.idx());
            Vertex svp5 = Vertex(meshs1.n_vertices() + meshs1.n_edges() + sf.idx());

            

            // find sub faces of sf
            // the creation of the new quads start from v2, so f1 is the quad v2->p2->p5->p1
            // thus the one that corresponds to the halfedge v2->p2
            gsSurfMesh::Halfedge she2 = find_halfedge(meshs2, sv2, svp2);
            gsSurfMesh::Face     sf1  = meshs2.face(she2);

            auto she3 = meshs2.next_halfedge(meshs2.next_halfedge(meshs2.opposite_halfedge(meshs2.next_halfedge(she2))));
            auto she4 = meshs2.next_halfedge(meshs2.next_halfedge(meshs2.opposite_halfedge(meshs2.next_halfedge(she3))));
            auto she1 = meshs2.next_halfedge(meshs2.next_halfedge(meshs2.opposite_halfedge(meshs2.next_halfedge(she4))));
            auto sf2  = meshs2.face(she3);
            auto sf3  = meshs2.face(she4);
            auto sf4  = meshs2.face(she1);

            // update grid points with appropriate uvs
            double sub_u_mid = sub_u_start + (sub_u_end - sub_u_start) / 2;
            double sub_v_mid = sub_v_start + (sub_v_end - sub_v_start) / 2;
            update_GridPoints_on_face(gps, meshs2, sf1, meshs3, sub_u_mid,    sub_u_end,  sub_v_start,  sub_v_mid,   svp1, f.idx());
            update_GridPoints_on_face(gps, meshs2, sf2, meshs3, sub_u_mid,    sub_u_end,  sub_v_mid,    sub_v_end,   svp5, f.idx());
            update_GridPoints_on_face(gps, meshs2, sf3, meshs3, sub_u_start,  sub_u_mid,  sub_v_mid,    sub_v_end,   svp4, f.idx());
            update_GridPoints_on_face(gps, meshs2, sf4, meshs3, sub_u_start,  sub_u_mid,  sub_v_start,  sub_v_mid,   sv1,  f.idx());
        }







        // grid points are set...

        for (auto gp : gps){

            Point x = gp.coords;

            // for debug
            // gsInfo << x << ": (" << gp.u << ", " << gp.v << ")\n";
            

            // find for this grid point the closest tri vertex
            // knn:
            Matrix queryPoints(3, 1);
            queryPoints << x[0], x[1], x[2];
            Matrixi indices;
            Matrix distances;
            int amount = 40;
            kdtree.query(queryPoints, amount, indices, distances);   // searching for 40 nearest neighbours
            // int closestVert_idx = indices(0);  // this index corresponds to the closest vert in the tri mesh

            // Vertex closestVert = Vertex(closestVert_idx);
            std::vector<Vertex> closestVerts;
            closestVerts.reserve(amount);
            for (int i = 0; i < indices.size(); ++i){
                closestVerts.emplace_back(indices(i));
            }


            // find actual closest point on the tri mesh 
            Point closestPoint = project_on_mesh(gp, mesh_quad_original, mesh_quad, mesh_tri, closestVerts, fix_degeneracies);

            // // override for debug
            // Point closestPoint = mesh_tri.position(closestVert);


            FitPoint fp;
            fp.coords = closestPoint;
            fp.u = gp.u;
            fp.v = gp.v;
            points_in_face[f.idx()].push_back(fp);

            // corrFace[closestVert] = f.idx();

            // is_used[closestVert] = true;


            // update fitting matrixes immediately
            current_offset++; // a grid point has been added
            points(0, point_counter) = fp.coords[0];
            points(1, point_counter) = fp.coords[1];
            points(2, point_counter) = fp.coords[2];

            uvs(0, point_counter) = fp.u;
            uvs(1, point_counter) = fp.v;
            point_counter++;
            
        }
    }
    // add last offset value
    offset(mesh_quad.n_faces()) = current_offset;




   

    // shrink to actual used size
    uvs.conservativeResize(gsEigen::NoChange, point_counter);
    points.conservativeResize(gsEigen::NoChange, point_counter);

    
    gsWriteParaview(mesh_tri, "mesh_in", { "v:corrFace"});
    gsWriteParaview(points, "points");


    gsFileData<> out;
    out.add(uvs, 0);// << Mpar;
    out.add(points, 1);// << fval;
    out.add(gsMatrix<index_t>(offset), 2); // << ; //conversion
    //out.add(mesh_quad);
    out.save("fitting_information");


    // write the multi basis
    gsFileData<> out_mesh;
    gsKnotVector<> kv(0,1,0,6);
    gsTensorBSplineBasis<2> tb(kv,kv);
    gsMultiBasis<> mb;
    for (unsigned int i = 0; i < mesh_quad.n_faces(); i++){
        mb.addBasis(tb.clone());
    }
    out_mesh << mb;

    
    
    gsG1SplineBasis mm(mesh_quad, use_c0);
    gsInfo << "initializng...\n";
    mm.initialize(5);
    gsInfo << "MultiBasis created!\n";
    
    gsSparseMatrix<> aa;
    mm.matrix().toSparseMatrix_into(aa);
    out_mesh << aa;
    out_mesh.save("quad_mesh_to_fit");


    // write the fully subdivided mesh
    gsWriteParaview(meshs3, "fully_subdivided_mesh");



    
    return 0;
}