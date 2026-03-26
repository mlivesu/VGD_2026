#undef NDEBUG
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

mat2d change_of_basis_matrix(const DrawableTrimesh<> & m)
{
    vec2d a0 = m.vert(0).rem_coord();
    vec2d a1 = m.vert(1).rem_coord();
    vec2d a2 = m.vert(2).rem_coord();

    vec2d b0 = m.vert(3).rem_coord();
    vec2d b1 = m.vert(4).rem_coord();
    vec2d b2 = m.vert(5).rem_coord();

    mat2d f0({a1-a0, a2-a0});
    mat2d f1({b1-b0, b2-b0});

    return f1*f0.inverse();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

vec2d singular_values(const mat2d & M)
{
    mat2d U,V;
    vec2d S;
    M.SVD(U,S,V);
    return S;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double Dirichlet   (const vec2d & sigma) { return std::pow(sigma[0],2) + std::pow(sigma[1],2); }
double SymDirichlet(const vec2d & sigma) { return std::pow(sigma[0],2) + std::pow(sigma[1],2) + std::pow(sigma[0],-2) + std::pow(sigma[1],-2); }
double MIPS        (const vec2d & sigma) { return sigma[0]/sigma[1] + sigma[1]/sigma[0]; }

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void plot_distortion(DrawableTrimesh<> & m, const vec2d & sigma, const uint dist_energy)
{
    switch(dist_energy)
    {
        case 0 : m.poly_set_color(Color::red_ramp_01((Dirichlet(sigma)-2.)/20.)); break;
        case 1 : m.poly_set_color(Color::red_ramp_01((SymDirichlet(sigma)-4.)/20.)); break;
        case 2 : m.poly_set_color(Color::red_ramp_01((MIPS(sigma)-2.)/20.)); break;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main()
{
    DrawableTrimesh<> m({vec3d(0,0,0),
                         vec3d(1,0,0),
                         vec3d(0,1,0),
                         vec3d(3,0,0),
                         vec3d(4,0,0),
                         vec3d(3,1,0)},{0,1,2,3,4,5});
    m.show_mesh_with_lighting(false);

    mat2d  M     = change_of_basis_matrix(m);
    vec2d  sigma = singular_values(M);

    GLcanvas gui;
    gui.push(&m);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>>(&m,&gui,"mesh"));
    gui.show_side_bar = true;

    // 0 : Dirichelt
    // 1 : Symmetric Dirichlet
    // 2 : MIPS
    uint dist_energy = 0;

    GLdouble zbuf = 0;
    uint vid;
    gui.callback_mouse_left_click = [&](int mod) -> bool
    {
        if(mod & GLFW_MOD_SHIFT)
        {
            vec3d click_3d;
            vec2d click_2d = gui.cursor_pos();
            if(gui.unproject(click_2d, click_3d))
            {
                vid  = m.pick_vert(click_3d);
                zbuf = gui.query_Z_buffer(click_2d);
            }
            return true;
        }
        return false;
    };

    gui.callback_mouse_moved = [&](double x_pos, double y_pos) -> bool
    {
        bool left  = glfwGetMouseButton(gui.window,GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS;
        bool shift = glfwGetKey(gui.window,GLFW_KEY_LEFT_SHIFT )==GLFW_PRESS ||
                     glfwGetKey(gui.window,GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS;

        if(left && shift)
        {
            vec3d p;
            gui.unproject(vec2d(x_pos,y_pos),zbuf, p);
            m.vert(vid) = p;
            m.updateGL();
            M     = change_of_basis_matrix(m);
            sigma = singular_values(M);
            plot_distortion(m,sigma,dist_energy);
            return true;
        }
        return false;
    };

    gui.callback_app_controls = [&]()
    {
        std::string s = "Change of Basis Matrix:\n" +
                        std::to_string(M(0,0)) + "\t" + std::to_string(M(0,1)) + "\n" +
                        std::to_string(M(1,0)) + "\t" + std::to_string(M(1,1)) + "\n\n" +
                        + "Singular Values:\n" +
                        "s0: " + std::to_string(sigma[0]) + "\n" +
                        "s0: " + std::to_string(sigma[1]) + "\n";
        ImGui::Text(s.c_str());

        if(M.det()<=0) ImGui::Text("FLIP!\n");

        s = "\nDirichlet :"  + std::to_string(Dirichlet(sigma)) +
            "\nSym Diric :"  + std::to_string(SymDirichlet(sigma)) +
            "\nMIPS :"  + std::to_string(MIPS(sigma));
        ImGui::Text(s.c_str());

        ImGui::Text("\nColor by\n");
        if(ImGui::RadioButton("Dir ",(dist_energy==0))) { dist_energy = 0; plot_distortion(m,sigma,dist_energy);}
        if(ImGui::RadioButton("SymD",(dist_energy==1))) { dist_energy = 1; plot_distortion(m,sigma,dist_energy);}
        if(ImGui::RadioButton("MIPS",(dist_energy==2))) { dist_energy = 2; plot_distortion(m,sigma,dist_energy);}
    };

    return gui.launch();
}
