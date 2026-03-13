#include <cinolib/gl/glcanvas.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/gl/surface_mesh_controls.h>

using namespace cinolib;

int main()
{
    DrawableTrimesh<> m("/Users/cino/Desktop/tmp_data/bunny.obj");

    for(uint nbr : m.adj_v2v(12045))
    {
        std::cout << nbr << std::endl;
    }

    GLcanvas gui;
    gui.push(&m);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>>(&m,&gui,"mesh"));
    gui.show_side_bar = true;
    return gui.launch();
}
