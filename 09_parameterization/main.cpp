#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/geometry/n_sided_poygon.h>
#include <cinolib/serialize_index.h>
#include <cinolib/grid_mesh.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void distortion(const DrawableTrimesh<>   & m_xyz,
                const DrawableTrimesh<>   & m_uv,
                      std::vector<double> & dist)
{
    dist.resize(m_xyz.num_polys());

    for(uint pid=0; pid<m_xyz.num_polys(); ++pid)
    {
        vec3d A0 = m_xyz.poly_vert(pid,0);
        vec3d A1 = m_xyz.poly_vert(pid,1);
        vec3d A2 = m_xyz.poly_vert(pid,2);

        vec3d Tu =  A1-A0;
        vec3d N  = (A1-A0).cross(A2-A0);
        vec3d Tv = Tu.cross(N);
        Tu.normalize();
        Tv.normalize();

        vec2d a0(0,0);
        vec2d a1(A1.dist(A0),0);
        vec2d a2((A2-A0).dot(Tu),(A2-A0).dot(Tv));

        vec2d b0 = m_uv.poly_vert(pid,0).rem_coord();
        vec2d b1 = m_uv.poly_vert(pid,1).rem_coord();
        vec2d b2 = m_uv.poly_vert(pid,2).rem_coord();

        mat2d uv0({a1-a0, a2-a0});
        mat2d uv1({b1-b0, b2-b0});
        mat2d J = uv1 * uv0.inverse();

        vec2d S;
        mat2d U,V;
        J.SVD(U,S,V);

        dist.at(pid) = (S[0]-S[1])*(S[0]-S[1]);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void laplacian(DrawableTrimesh<> & m, const uint vid)
{
    vec3d p(0,0,0);
    for(uint nbr : m.adj_v2v(vid))
    {
        p += m.vert(nbr);
    }
    m.vert(vid) = p / m.vert_valence(vid);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void laplacian_all(DrawableTrimesh<> & m, const uint n_iter)
{
    for(uint i=0; i<n_iter; ++i)
    {
        for(uint vid=0; vid<m.num_verts(); ++vid)
        {
            if(m.vert_is_boundary(vid)) continue;
            laplacian(m,vid);
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void quadmesh(const DrawableTrimesh<>  & m_xyz,
              const DrawableTrimesh<>  & m_uv,
                    DrawableQuadmesh<> & q)
{
    Octree o;
    o.build_from_mesh_polys(m_uv);

    uint n_quads = 100;

    grid_mesh(n_quads,n_quads,q);
    q.center_bbox();
    q.normalize_bbox();
    q.scale(2);

    for(uint i=0; i<=n_quads; ++i)
    for(uint j=0; j<=n_quads; ++j)
    {
        uint vid = serialize_2D_index(i,j,n_quads+1);

        uint pid;
        double d;
        vec3d pos;
        o.closest_point(q.vert(vid),pid,pos,d);

        double w[3];
        triangle_barycentric_coords(m_uv.poly_vert(pid,0),
                                    m_uv.poly_vert(pid,1),
                                    m_uv.poly_vert(pid,2),q.vert(vid),w);

        q.vert(vid) = w[0] * m_xyz.poly_vert(pid,0) +
                      w[1] * m_xyz.poly_vert(pid,1) +
                      w[2] * m_xyz.poly_vert(pid,2);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    DrawableTrimesh<> m(argv[1]);
    m.center_bbox();
    m.normalize_bbox();
    m.updateGL();

    DrawableTrimesh<>  m_ref = m;
    DrawableQuadmesh<> q;

    std::vector<uint>  boundary = m.get_ordered_boundary_vertices();
    std::vector<vec3d> circle   = n_sided_polygon(boundary.size(),CIRCLE);

    GLcanvas gui;
    gui.push(&m);
    gui.push(&q);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>>(&m,&gui));

    std::vector<double> dist;

    int  tot_iters = 0;
    int  iters = 1;
    gui.callback_key_pressed = [&](int key, int modifiers) -> bool
    {
        if(key==GLFW_KEY_I)
        {
            for(uint i=0; i<boundary.size(); ++i)
            {
                m.vert(boundary.at(i)) = circle.at(i);
            }
            m.updateGL();
        }
        if(key==GLFW_KEY_SPACE)
        {
            laplacian_all(m,iters++);
            m.updateGL();
            tot_iters += iters;
            std::cout << "#iters " << tot_iters << std::endl;
            return true;
        }
        if(key==GLFW_KEY_T)
        {
            m.copy_xyz_to_uvw(UVW_param);
            m.show_texture2D(TEXTURE_2D_ISOLINES, 5.0);
            quadmesh(m_ref,m,q);
            q.update_normals();
            q.updateGL();
            return true;
        }
        if(key==GLFW_KEY_U)
        {
            m.vector_verts() = m_ref.vector_verts();
            m.update_normals();
            m.updateGL();
            return true;
        }
        return false;
    };

    gui.callback_mouse_left_click = [&](int modifiers) -> bool
    {
        if(modifiers & GLFW_MOD_SHIFT)
        {
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click,p))
            {
                uint vid = m.pick_vert(p);
                laplacian(m,vid);
                m.updateGL();
            }
        }
        return false;
    };

    gui.callback_app_controls = [&]()
    {
        if(ImGui::Button("Distortion"))
        {
            if(dist.empty()) distortion(m_ref,m,dist);

            double E = 0;
            for(uint pid=0; pid<m.num_polys(); ++pid)
            {
                E += dist.at(pid) * m.poly_area(pid);;
                m.poly_data(pid).color = Color::red_ramp_01(std::min(dist.at(pid),50.0));
                m_ref.poly_data(pid).color = m.poly_data(pid).color;
            }
            m_ref.show_poly_color();
            m.show_poly_color();
            std::cout << "E_conf : " << E << std::endl;
        }
    };

    return gui.launch();
}
