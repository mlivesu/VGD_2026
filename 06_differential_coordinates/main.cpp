#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/linear_solvers.h>
#include <cinolib/laplacian.h>
#include <cinolib/profiler.h>

using namespace cinolib;

Profiler p;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<vec3d> diff_coords(DrawableTrimesh<> & m)
{
    std::vector<vec3d> coords(m.num_verts(),vec3d(0,0,0));
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        for(uint nbr : m.adj_v2v(vid))
        {
            coords.at(vid) += (m.vert(vid) - m.vert(nbr));
        }
    }
    return coords;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void surface_reconstruction(DrawableTrimesh<>        & m,
                            const std::vector<vec3d> & diff_coords,
                            const std::vector<uint>  & anchors)
{
    Eigen::VectorXd rhs_x(m.num_verts()), rhs_y(m.num_verts()), rhs_z(m.num_verts());
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        rhs_x[vid] = diff_coords.at(vid).x();
        rhs_y[vid] = diff_coords.at(vid).y();
        rhs_z[vid] = diff_coords.at(vid).z();
    }

    std::map<uint,double> bcx,bcy,bcz;
    for(uint vid : anchors)
    {
        bcx[vid] = m.vert(vid).x();
        bcy[vid] = m.vert(vid).y();
        bcz[vid] = m.vert(vid).z();
    }

    Eigen::SparseMatrix<double> A = -laplacian(m, UNIFORM, 1);
    Eigen::VectorXd x,y,z;
    solve_square_system_with_bc(A, rhs_x, x, bcx);
    solve_square_system_with_bc(A, rhs_y, y, bcy);
    solve_square_system_with_bc(A, rhs_z, z, bcz);

    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        m.vert(vid) = vec3d(x[vid], y[vid], z[vid]);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    DrawableTrimesh<> m(argv[1]);
    m.normalize_bbox();
    m.updateGL();

    std::vector<uint>  anchors;
    std::vector<vec3d> coords = diff_coords(m);

    DrawableSegmentSoup ss;
    ss.thickness = 5;
    ss.use_gl_lines = true;
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        ss.push_seg(m.vert(vid),m.vert(vid)+coords.at(vid));
    }

    GLcanvas gui;
    gui.push(&m);
    gui.push(&ss);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>>(&m,&gui));

    gui.callback_mouse_left_click = [&](int mod) -> bool
    {
        if(mod & GLFW_MOD_SHIFT)
        {
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click, p))
            {
                uint vid = m.pick_vert(p);
                anchors.push_back(vid);
                m.vert_data(vid).color = Color::BLUE();
                m.show_vert_color();
                m.updateGL();
            }
            return true;
        }
        return false;
    };

    gui.callback_key_pressed = [&](int key, int mod) -> bool
    {
        if(key==GLFW_KEY_SPACE)
        {
            p.push("Reconstruction with hard BCs");
            surface_reconstruction(m,coords,anchors);
            p.pop();
            m.updateGL();
            return true;
        }
        return false;
    };

    return gui.launch();
}
