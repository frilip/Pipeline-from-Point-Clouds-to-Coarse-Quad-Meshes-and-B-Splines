#include <iostream>
#include <clocale>

#include "functions.h"

#include <bits/stdc++.h>



std::vector <size_t> path_on_field(VertexFieldGraph<TraceMesh> &VFGraph, 
                                   const typename TraceMesh::ScalarType &Drift,
                                   const typename TraceMesh::ScalarType &MaxWeight,
                                   size_t start,
                                   size_t end)
{
    // start, finish are vertices -> find nodes:
    int min_len = std::numeric_limits<int>::max();
    std::vector <size_t> min_path;
    std::vector <size_t> path;
    for (size_t n2 = end * 4; n2 < end * 4 + 4; n2++){
            
        VFGraph.Select(n2);
        
        //VFGraph.Deselect(n2);
        // if (path.size() < min_len){
        //     min_len = path.size();
        //     min_path = path;
        // }
    }
    for (size_t n1 = start * 4; n1 < start * 4 + 4; n1++){
        TraceDijkstraPath(VFGraph, n1, Drift, MaxWeight, path);
        std::cout << ": " << path.size() << "\n";
        if (path.size() < min_len){
            min_len = path.size();
            min_path = path;
        }

        for (auto p : path){
            std::cout << p / 4 << ", ";
        }
        std::cout << "\n";
    }

    for (size_t n2 = end * 4; n2 < end * 4 + 4; n2++){
        VFGraph.DeSelect(n2);
    }


    return min_path;
}           






// int main(int argc, char *argv[])
// {

//     TraceMesh traceTrimesh;


//     PolyMesh quadmesh;
//     std::vector<std::vector<size_t>> quadmeshPartitions;
//     std::vector<std::vector<size_t>> quadmeshCorners;
//     std::vector<int> ilpResult;

//     if (argc==1)
//     {
//         std::cout<<"Please specify a mesh (OBJ or PLY) as argument"<<std::endl;
//         exit(0);
//     }

//     std::string name=std::string(argv[1]);
//     std::string tri_filename = name + "_rem.obj";
//     std::string quad_filename = name + "_quadrangulation.obj";
//     std::string field_fn = name + "_rem.rosy";


//     //Use "." as decimal separator
//     std::setlocale(LC_NUMERIC, "en_US.UTF-8");



//     bool loaded=traceTrimesh.LoadMesh(tri_filename);
//     traceTrimesh.UpdateAttributes();
//     bool loaded_f = traceTrimesh.LoadField(field_fn);
//     //trimesh.UpdateDataStructures();

//     if (!loaded)
//     {
//         std::cout<<"Wrong mesh filename"<<std::endl;
//         exit(0);
//     }

//     int mask = 0;
//     int err = vcg::tri::io::ImporterOBJ<PolyMesh>::Open(quadmesh, quad_filename.c_str(), mask);

//     if(err) {
//         printf("Error loading quad OBJ\n");
//     }

//     std::cout<<"Tri loaded "<<traceTrimesh.fn<<" faces and "<<traceTrimesh.vn<<" vertices"<<std::endl;
//     std::cout<<"Quad loaded "<<quadmesh.fn<<" faces and "<<quadmesh.vn<<" vertices"<<std::endl;

//     std::vector<size_t> PathN;
//     TraceMesh::ScalarType Drift=100;
//     TraceMesh::ScalarType MaxWeight=std::numeric_limits<float>::max();
//     //initializing graph
//     VertexFieldGraph<TraceMesh> VFGraph(traceTrimesh);
//     VFGraph.InitGraph(false);

//     std::cout << "ALl good here\n";

//     // tri vert pairs that correspond to quad edges
//     std::ifstream in ("quad_to_tri_edges_map.txt");
//     int quad_edge_number;
//     in >> quad_edge_number;
//     std::cout << "number of quad edges: " << quad_edge_number << "\n";
//     size_t start, end;
//     std::vector <size_t> path;
//     for (int i = 0; i < quad_edge_number; i++){
//         in >> start;
//         in >> end;
//         std::cout << "starting path..\nstart: " << start << " end: " << end << "\n";
//         path = path_on_field(VFGraph, Drift, MaxWeight, start, end);
//         std::cout << "path found!\n";

//         std::cout << "size = " << path.size();
//         for (auto p : path){
//             std::cout << p << ", ";
//         }            
//     }
// }


int main(int argc, char *argv[])
{
    Parameters parameters;

    FieldTriMesh trimesh;
    TraceMesh traceTrimesh;

    TriangleMesh trimeshToQuadrangulate;
    std::vector<std::vector<size_t>> trimeshPartitions;
    std::vector<std::vector<size_t>> trimeshCorners;
    std::vector<std::pair<size_t,size_t> > trimeshFeatures;
    std::vector<size_t> trimeshFeaturesC;

    PolyMesh quadmesh;
    std::vector<std::vector<size_t>> quadmeshPartitions;
    std::vector<std::vector<size_t>> quadmeshCorners;
    std::vector<int> ilpResult;

    if (argc==1)
    {
        std::cout<<"Please specify a mesh (OBJ or PLY) as argument"<<std::endl;
        exit(0);
    }

    std::string name = std::string(argv[1]);
    std::string meshFilename = name + ".obj";
    std::string sharpFilename;
    std::string fieldFilename;

    //Use "." as decimal separator
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");

    std::cout<<"Reading input..."<<std::endl;
    loadConfigFile("basic_setup.txt", parameters);

    for (int i=2;i<argc;i++)
    {
        int position;

        std::string pathTest=std::string(argv[i]);

        position=pathTest.find(".sharp");
        if (position!=-1)
        {
            sharpFilename=pathTest;
            parameters.hasFeature=true;
            continue;
        }

        position=pathTest.find(".txt");
        if (position!=-1)
        {
           loadConfigFile(pathTest.c_str(), parameters);
           continue;
        }

        position=pathTest.find(".rosy");
        if (position!=-1)
        {
           fieldFilename=pathTest;
           parameters.hasField=true;
           continue;
        }
    }


    std::cout<<"Loading:"<<meshFilename.c_str()<<std::endl;

    bool allQuad;
    bool loaded=trimesh.LoadTriMesh(meshFilename,allQuad);
    trimesh.UpdateDataStructures();

    if (!loaded)
    {
        std::cout<<"Wrong mesh filename"<<std::endl;
        exit(0);
    }

    std::cout<<"Loaded "<<trimesh.fn<<" faces and "<<trimesh.vn<<" vertices"<<std::endl;

    std::cout<<std::endl<<"--------------------- 1 - Remesh and field ---------------------"<<std::endl;
    remeshAndField(trimesh, parameters, meshFilename, sharpFilename, fieldFilename);

    std::cout<<std::endl<<"--------------------- 2 - Tracing ---------------------"<<std::endl;
    trace(meshFilename, traceTrimesh);

    std::cout<<std::endl<<"--------------------- 3 - Quadrangulation ---------------------"<<std::endl;
    bool DoSmoothing = false;
    quadrangulate(meshFilename, trimeshToQuadrangulate, quadmesh, trimeshPartitions, trimeshCorners, trimeshFeatures, trimeshFeaturesC, quadmeshPartitions, quadmeshCorners, ilpResult, parameters, DoSmoothing);




    std::vector<size_t> PathN;
    TraceMesh::ScalarType Drift=100;
    TraceMesh::ScalarType MaxWeight=std::numeric_limits<float>::max();
    //initializing graph
    VertexFieldGraph<TraceMesh> VFGraph(traceTrimesh);
    VFGraph.InitGraph(false);


    // // tri vert pairs that correspond to quad edges
    // std::ifstream in ("quad_to_tri_edges_map.txt");
    // int quad_edge_number;
    // in >> quad_edge_number;
    // size_t start, end;
    // std::vector <size_t> path;
    // for (int i = 0; i < quad_edge_number; i++){
    //     in >> start;
    //     in >> end;
    //     path = path_on_field(VFGraph, Drift, MaxWeight, start, end);

    //     for (auto p : path){
    //         std::cout << p / 4 << ", ";
    //     }            
    //     std::cout << "\n";
    //     int a;
    //     std::cin >> a;
    // }



    // path from vertices: start -> end
    size_t start = 1212; 
    size_t end = 4055;
    PathN = path_on_field(VFGraph, Drift, MaxWeight, start, end);
    std::cout << start << ", " << end << "\n";
    std::cout << PathN.size() << "\n";
    for (auto p : PathN){
        std::cout << p / 4 << ", ";
    }
    std::cout << "\n";
    std::cout << "AAAAAAAA\n";
    int a;
    std::cin >> a;
}