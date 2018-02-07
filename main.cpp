#include <iostream>
#include "ICP.h"
#include "io_obj.h"
#include <chrono>
#include <vtkVersion.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

int main (int argc, char const ** argv)
{   
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> Vertices;
    std::string file_source;
    std::string file_target;
    std::string file_source_reg;

    ///--- Option parsing    
    if(argc==1){
        std::cout << "examples: " << std::endl;
        std::cout << "./sparseicp data/hippo_0.vtk data/hippo_1.vtk data/hippo_0.vtk__REG__hippo_1.vtk" << std::endl;
        // std::cout << "./sparseicp data/bunny.obj data/bunny_cut2.obj data/bunny.obj__REG__bunny_cut2.obj" << std::endl;
        exit(0);
    }
    else if(argc==4)
    {
        file_source = argv[1];
        file_target = argv[2];
        file_source_reg = argv[3];        
    } 
    else 
    {
        std::cout << "argument error!" << std::endl;
        exit(0);
    }
    
    vtkSmartPointer<vtkPolyData> input;
    vtkSmartPointer<vtkPolyDataReader> reader1 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader1->SetFileName(argv[1]);
    reader1->Update();
    input = reader1->GetOutput();

    vtkSmartPointer<vtkPolyData> target;
    vtkSmartPointer<vtkPolyDataReader> reader2 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(argv[2]);
    reader2->Update();
    target = reader2->GetOutput();

    double p[3];
    Vertices vertices_source;
    vtkIdType n_vertices = input->GetPoints()->GetNumberOfPoints();
    vertices_source.resize(3, n_vertices);
    for(unsigned int i=0; i<n_vertices; ++i)
    {
      input->GetPoint(i,p);
      vertices_source(0, i) = p[0];
      vertices_source(1, i) = p[1];
      vertices_source(2, i) = p[2];
    }
    std::cout << "source: " << vertices_source.rows() << "x" << vertices_source.cols() << std::endl;

    Vertices vertices_target;
    n_vertices = target->GetPoints()->GetNumberOfPoints();
    vertices_target.resize(3, n_vertices);
    for(unsigned int i=0; i<n_vertices; ++i)
    {
      target->GetPoint(i,p);
      vertices_target(0, i) = p[0];
      vertices_target(1, i) = p[1];
      vertices_target(2, i) = p[2];
    }
    std::cout << "target: " << vertices_target.rows() << "x" << vertices_target.cols() << std::endl;

/*
    ///--- Model that will be rigidly transformed
    Vertices vertices_source;
    read_obj(vertices_source, file_source);
    std::cout << "source: " << vertices_source.rows() << "x" << vertices_source.cols() << std::endl;
    
    ///--- Model that source will be aligned to
    Vertices vertices_target;
    read_obj(vertices_target, file_target);
    std::cout << "target: " << vertices_target.rows() << "x" << vertices_target.cols() << std::endl;
*/
    ///--- Execute registration
    auto tic = std::chrono::steady_clock::now(); 
        SICP::Parameters pars;
        pars.p = 1.5;
        pars.max_icp = 2500;
        pars.max_inner = 1;
        pars.mu = 1.;
        //pars.use_penalty = true;
        pars.print_icpn = true;
        SICP::point_to_point(vertices_source, vertices_target, pars);
    auto toc = std::chrono::steady_clock::now();

    ///--- Write result to file
    // write_obj_replaceverts(file_source, vertices_source, file_source_reg);
    n_vertices = input->GetPoints()->GetNumberOfPoints();
    for(unsigned int i=0; i<n_vertices; ++i)
    {
      input->GetPoints()->SetPoint(i,
                                   vertices_source(0, i),
                                   vertices_source(1, i),
                                   vertices_source(2, i));
    }
    vtkSmartPointer<vtkPolyDataWriter> writer
      = vtkSmartPointer<vtkPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(input);
#else
    writer->SetInputData(input);
#endif
    writer->SetFileName(argv[3]);
    writer->Update();
    
    ///--- Print execution time
    double time_ms = std::chrono::duration <double, std::milli> (toc-tic).count();
    std::cout << "sparseicp registered source to target in: " << time_ms << "ms" << std::endl;
    return 0;
}
