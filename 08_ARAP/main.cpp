#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/linear_solvers.h>
#include <cinolib/laplacian.h>
#include <cinolib/parallel_for.h>
#include <cinolib/profiler.h>

using namespace cinolib;

Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
Eigen::SparseMatrix<double> A;
std::vector<double>         w;
std::vector<mat3d>          R;
std::map<uint,double>       bc_x;
std::map<uint,double>       bc_y;
std::map<uint,double>       bc_z;
DrawableTrimesh<>           m_local;

std::vector<uint> handles;
enum
{
    ADD_HANDLES = 0,
    MOV_HANDLES = 1
};
int  ARAP_mode        = ADD_HANDLES;
bool show_local_mesh  = false;
bool hard_constraints = true;
Profiler p;

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

void init(DrawableTrimesh<> & m)
{
    p.push("Init");

    R.resize(m.num_verts());

    w.resize(m.num_edges());
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        w.at(eid) = m.edge_weight_cotangent(eid);
    }

    if(hard_constraints)
    {
        A = -laplacian(m,COTANGENT,1);
    }
    else
    {
        uint nv = m.num_verts();
        uint nh = handles.size();
        std::vector<Eigen::Triplet<double>> entries;
        for(uint vid=0; vid<m.num_verts(); ++vid)
        {
            double diag = 0;
            for(uint nbr : m.adj_v2v(vid))
            {
                int eid = m.edge_id(vid,nbr);
                entries.emplace_back(vid,nbr,-w.at(eid));
                diag += w.at(eid);
            }
            entries.emplace_back(vid,vid,diag);
        }
        for(uint i=0; i<handles.size(); ++i)
        {
            entries.emplace_back(nv+i,handles.at(i),1.0);
        }
        A.resize(nv+nh,nv);
        A.setFromTriplets(entries.begin(), entries.end());
        solver.derived().compute(A.transpose()*A);
    }
    p.pop();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void local_step(const DrawableTrimesh<> & m_cur,
                const DrawableTrimesh<> & m_ref)
{
    PARALLEL_FOR(0,m_cur.num_verts(),1000,[&](uint vid)
    {
        mat3d cov = mat3d::ZERO();
        for(uint nbr : m_cur.adj_v2v(vid))
        {
            int   eid   = m_cur.edge_id(vid,nbr);
            vec3d e_ref = (m_ref.vert(vid) - m_ref.vert(nbr));
            vec3d e_cur = (m_cur.vert(vid) - m_cur.vert(nbr));

            cov += w.at(eid) * (e_cur * e_ref.transpose());
        }

        mat3d U,V;
        vec3d S;
        cov.SVD(U,S,V);
        mat3d I = mat3d::DIAG(1);
        I(2,2) = (U*V.transpose()).det();
        R.at(vid) = U*I*V.transpose();
    });
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void build_local_mesh(const DrawableTrimesh<> & m,
                      const DrawableTrimesh<> & m_ref)
{
    m_local = DrawableTrimesh<>();
    for(uint vid=0; vid<m_ref.num_verts(); ++vid)
    {
        std::map<uint,uint> v_map;
        v_map[vid] = m_local.vert_add(m.vert(vid));
        for(uint nbr : m_ref.adj_v2v(vid))
        {
            vec3d P = m_ref.vert(nbr);
            P -= m_ref.vert(vid);
            P  = R.at(vid)*P;
            P += m_ref.vert(vid) + (m.vert(vid) - m_ref.vert(vid));
            v_map[nbr] = m_local.vert_add(P);
        }
        for(uint pid : m_ref.adj_v2p(vid))
        {
            m_local.poly_add({v_map.at(m_ref.poly_vert_id(pid,0)),
                              v_map.at(m_ref.poly_vert_id(pid,1)),
                              v_map.at(m_ref.poly_vert_id(pid,2))});
        }
    }
    m_local.poly_set_color(Color::PASTEL_GREEN());
    m_local.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void global_step(      DrawableTrimesh<> & m_cur,
                 const DrawableTrimesh<> & m_ref)
{
    uint nv = m_cur.num_verts();
    uint nh = (hard_constraints) ? 0 : handles.size();

    Eigen::VectorXd rhs_x(nv+nh);
    Eigen::VectorXd rhs_y(nv+nh);
    Eigen::VectorXd rhs_z(nv+nh);

    for(uint vid=0; vid<m_cur.num_verts(); ++vid)
    {
        vec3d rhs(0,0,0);
        for(uint nbr : m_cur.adj_v2v(vid))
        {
            mat3d Ravg = (R.at(vid)+R.at(nbr))/2.0;
            vec3d e    = (m_ref.vert(vid) - m_ref.vert(nbr));

            int eid = m_cur.edge_id(vid,nbr);
            rhs += w.at(eid) * Ravg * e;
        }
        rhs_x[vid] = rhs.x();
        rhs_y[vid] = rhs.y();
        rhs_z[vid] = rhs.z();
    }

    Eigen::VectorXd x,y,z;
    if(hard_constraints)
    {
        solve_square_system_with_bc(A, rhs_x, x, bc_x);
        solve_square_system_with_bc(A, rhs_y, y, bc_y);
        solve_square_system_with_bc(A, rhs_z, z, bc_z);
    }
    else
    {
        for(uint i=0; i<handles.size(); ++i)
        {
            uint vid = handles.at(i);
            rhs_x[nv+i] = bc_x.at(vid);
            rhs_y[nv+i] = bc_y.at(vid);
            rhs_z[nv+i] = bc_z.at(vid);
        }
        x = solver.solve(A.transpose()*rhs_x).eval();
        y = solver.solve(A.transpose()*rhs_y).eval();
        z = solver.solve(A.transpose()*rhs_z).eval();
    }

    for(uint vid=0; vid<m_cur.num_verts(); ++vid)
    {
        m_cur.vert(vid) = vec3d(x[vid],y[vid],z[vid]);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void ARAP(      DrawableTrimesh<> & m_cur,
          const DrawableTrimesh<> & m_ref,
          const uint n_iters)
{
    p.push("ARAP");
    for(uint i=0; i<n_iters; ++i)
    {
        local_step(m_cur,m_ref);
        global_step(m_cur,m_ref);
    }
    p.pop();
    m_cur.update_normals();
    m_cur.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    DrawableTrimesh<> m(argv[1]);
    DrawableTrimesh<> m_ref = m;
    m_ref.poly_set_color(Color::PASTEL_GREEN());

    init(m_ref);

    GLcanvas gui;
    gui.push(&m);
    gui.push(&m_local);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>>(&m,&gui,"Deformed mesh"));
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>>(&m_local,&gui,"Local mesh"));
    gui.show_side_bar      = true;
    gui.depth_cull_markers = false;

    gui.callback_key_pressed = [&](int key, int mod) -> bool
    {
        if(key==GLFW_KEY_SPACE)
        {
            ARAP(m,m_ref,1);
            if(show_local_mesh) build_local_mesh(m,m_ref);
            return true;
        }
        return false;
    };

    gui.callback_app_controls = [&]()
    {
        if(ImGui::RadioButton("Add  handles",ARAP_mode==ADD_HANDLES)) { ARAP_mode = 0;              }
        if(ImGui::RadioButton("Move handles",ARAP_mode==MOV_HANDLES)) { ARAP_mode = 1; init(m_ref); }

        if(ImGui::Checkbox("Show Local Mesh",&show_local_mesh))
        {
            if(show_local_mesh) build_local_mesh(m,m_ref);
            m_local.show_mesh(show_local_mesh);
        }
        if(ImGui::Checkbox("Hard Constraints",&hard_constraints))
        {
            init(m_ref);
        }
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
                if(ARAP_mode==ADD_HANDLES)
                {
                    uint vid  = m.pick_vert(click_3d);
                    handles.push_back(vid);
                    bc_x[vid] = m.vert(vid).x();
                    bc_y[vid] = m.vert(vid).y();
                    bc_z[vid] = m.vert(vid).z();
                    gui.push_marker(m.vert(vid), "", Color::BLUE(), 10);
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
        if(ARAP_mode==ADD_HANDLES) return false;

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
            ARAP(m,m_ref,4);
            if(show_local_mesh) build_local_mesh(m,m_ref);
            gui.pop_all_markers();
            for(uint vid : handles) gui.push_marker(m.vert(vid), "", Color::BLUE(), 10);
            gui.draw();
            return true;
        }
        return false;
    };
    return gui.launch();
}
