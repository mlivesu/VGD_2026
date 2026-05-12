#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/laplacian.h>
#include <cinolib/linear_solvers.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

using namespace cinolib;

std::map<uint,double> bc_x;
std::map<uint,double> bc_y;
std::map<uint,double> bc_z;
std::vector<uint>     handles;
std::vector<vec3d>    diff_coords;
enum
{
    ADD_HANDLES = 0,
    MOV_HANDLES = 1
};
int mode = ADD_HANDLES;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

uint pick_handle(const Trimesh<> & m,
                 const vec3d     & click)
{
    return *std::min_element(handles.begin(),
                             handles.end(),
                             [&](const uint hi, const uint hj)
                             {
                                 return click.dist(m.vert(hi))<click.dist(m.vert(hj));
                             });
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<vec3d> differential_coordinates(const DrawableTrimesh<> & m)
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

void surface_reconstruction(DrawableTrimesh<> & m)
{
    Eigen::SparseMatrix<double> A = -laplacian(m, UNIFORM, 1);

    Eigen::VectorXd rhs_x(m.num_verts());
    Eigen::VectorXd rhs_y(m.num_verts());
    Eigen::VectorXd rhs_z(m.num_verts());
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        rhs_x[vid] = diff_coords.at(vid).x();
        rhs_y[vid] = diff_coords.at(vid).y();
        rhs_z[vid] = diff_coords.at(vid).z();
    }

    Eigen::VectorXd x,y,z;
    solve_square_system_with_bc(A, rhs_x, x, bc_x);
    solve_square_system_with_bc(A, rhs_y, y, bc_y);
    solve_square_system_with_bc(A, rhs_z, z, bc_z);

    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        m.vert(vid).x() = x[vid];
        m.vert(vid).y() = y[vid];
        m.vert(vid).z() = z[vid];
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main()
{
    DrawableTrimesh<> m("/Users/cino/Desktop/tmp_data/maxFace.obj");
    m.normalize_bbox();
    m.updateGL();    

    GLcanvas gui;
    gui.push(&m);

    diff_coords = differential_coordinates(m);
    for(uint vid : m.get_boundary_vertices())
    {
        handles.push_back(vid);
        bc_x[vid] = m.vert(vid).x();
        bc_y[vid] = m.vert(vid).y();
        bc_z[vid] = m.vert(vid).z();
    }
    gui.pop_all_markers();
    for(uint vid : handles) gui.push_marker(m.vert(vid), "", Color::BLUE(), 5);

    gui.callback_app_controls = [&]()
    {
        if(ImGui::RadioButton("Add  handles",mode==ADD_HANDLES)) { mode = 0; }
        if(ImGui::RadioButton("Move handles",mode==MOV_HANDLES)) { mode = 1; }
    };

    GLdouble zbuf = 0;
    int curr_handle;
    gui.callback_mouse_left_click = [&](int mod) -> bool
    {
        if(mod & GLFW_MOD_SHIFT)
        {
            vec3d click_3d;
            vec2d click_2d = gui.cursor_pos();
            if(gui.unproject(click_2d, click_3d))
            {
                if(mode==ADD_HANDLES)
                {
                    uint vid  = m.pick_vert(click_3d);
                    handles.push_back(vid);
                    bc_x[vid] = m.vert(vid).x();
                    bc_y[vid] = m.vert(vid).y();
                    bc_z[vid] = m.vert(vid).z();
                    gui.push_marker(m.vert(vid), "", Color::BLUE(), 5);
                }
                else
                {
                    curr_handle = pick_handle(m, click_3d);
                    zbuf = gui.query_Z_buffer(click_2d);
                }
            }
            m.updateGL();
            return true;
        }
        return false;
    };

    gui.callback_mouse_moved = [&](double x_pos, double y_pos) -> bool
    {
        if(mode==ADD_HANDLES) return false;

        bool left  = glfwGetMouseButton(gui.window,GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS;
        bool shift = glfwGetKey(gui.window,GLFW_KEY_LEFT_SHIFT )==GLFW_PRESS ||
                     glfwGetKey(gui.window,GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS;

        if(left && shift)
        {
            vec3d p;
            gui.unproject(vec2d(x_pos,y_pos),zbuf, p);
            vec3d delta = p - m.vert(curr_handle);
            bc_x[curr_handle] += delta.x();
            bc_y[curr_handle] += delta.y();
            bc_z[curr_handle] += delta.z();
            surface_reconstruction(m);
            m.updateGL();
            gui.pop_all_markers();
            for(uint vid : handles) gui.push_marker(m.vert(vid), "", Color::BLUE(), 5);
            gui.draw();
            return true;
        }
        return false;
    };

    return gui.launch();
}
