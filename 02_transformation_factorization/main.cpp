#undef NDEBUG
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void matrix_to_axis_angle(const mat3d & R, vec3d & axis, float & angle)
{
    // WARNING: degenerate for 0° and 180°
    angle = acos((R.trace()-1.0)/2.0);
    axis  = vec3d(R(2,1)-R(1,2),
                  R(0,2)-R(2,0),
                  R(1,0)-R(0,1)) * 1.0/(2.0*sin(angle));
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main()
{
    DrawableTrimesh<> beg("/Users/cino/Desktop/tmp_data/bunny.obj");
    beg.show_wireframe(false);
    beg.center_bbox();
    DrawableTrimesh<> cur = beg;
    DrawableTrimesh<> end = beg;

    float Euler[3]     = { 0, 0, 0 };
    float Euler_old[3] = { 0, 0, 0 };
    float scale[3]     = { 1, 1, 1 };
    float scale_old[3] = { 1, 1, 1 };
    float t_matrix     = 0.f;
    float t_QR        = 0.f;
    mat3d M            = mat3d::DIAG(1.0);
    vec3d T            = vec3d(cur.bbox().delta_x()*2,0,0);
    vec3d c_end        = cur.bbox().center()+T;

    end.translate(T);
    end.poly_set_color(Color::RED());

    GLcanvas gui;
    gui.show_side_bar = true;
    gui.push(&cur);
    gui.push(&end);

    gui.callback_app_controls = [&]()
    {
        if(ImGui::SliderFloat("Rx",&Euler[0],0,360))
        {
            float delta   = Euler[0] - Euler_old[0];
            mat3d R_delta = mat3d::ROT_3D(vec3d(1,0,0),to_rad(delta));
            Euler_old[0]  = Euler[0];
            M = R_delta*M;
            end.translate(-c_end);
            end.transform(R_delta);
            end.translate(c_end);
            end.updateGL();
        }
        if(ImGui::SliderFloat("Ry",&Euler[1],0,360))
        {
            float delta   = Euler[1] - Euler_old[1];
            mat3d R_delta = mat3d::ROT_3D(vec3d(0,1,0),to_rad(delta));
            Euler_old[1]  = Euler[1];
            M = R_delta*M;
            end.translate(-c_end);
            end.transform(R_delta);
            end.translate(c_end);
            end.updateGL();
        }
        if(ImGui::SliderFloat("Rz",&Euler[2],0,360))
        {
            float delta   = Euler[2] - Euler_old[2];
            mat3d R_delta = mat3d::ROT_3D(vec3d(0,0,1),to_rad(delta));
            Euler_old[2]  = Euler[2];
            M = R_delta*M;
            end.translate(-c_end);
            end.transform(R_delta);
            end.translate(c_end);
            end.updateGL();
        }
        if(ImGui::SliderFloat("Sx",&scale[0],0.1,2))
        {
            float delta   = scale[0]/scale_old[0];
            mat3d S_delta = mat3d::DIAG(vec3d(delta,1,1));
            scale_old[0]  = scale[0];
            M = S_delta*M;
            end.translate(-c_end);
            end.transform(S_delta);
            end.translate(c_end);
            end.updateGL();
        }
        if(ImGui::SliderFloat("Sy",&scale[1],0.1,2))
        {
            float delta   = scale[1]/scale_old[1];
            mat3d S_delta = mat3d::DIAG(vec3d(1,delta,1));
            scale_old[1]  = scale[1];
            M = S_delta*M;
            end.translate(-c_end);
            end.transform(S_delta);
            end.translate(c_end);
            end.updateGL();
        }
        if(ImGui::SliderFloat("Sz",&scale[2],0.1,2))
        {
            float delta   = scale[2]/scale_old[2];
            mat3d S_delta = mat3d::DIAG(vec3d(1,1,delta));
            scale_old[2]  = scale[2];
            M = S_delta*M;
            end.translate(-c_end);
            end.transform(S_delta);
            end.translate(c_end);
            end.updateGL();
        }
        ImGui::NewLine();
        ImGui::Text("Interpolate Matrix");
        if(ImGui::SliderFloat("##Interpolate Matrix",&t_matrix,0,1))
        {
            mat3d M_int = mat3d::DIAG(1.f)*(1-t_matrix) + M*t_matrix;
            cur = beg;
            cur.transform(M_int);
            vec3d T_int = T*t_matrix;
            cur.translate(T_int);
            cur.poly_set_color(Color::red_ramp_01(t_matrix));
            cur.updateGL();
        }
        ImGui::Text("Interpolate QR");
        if(ImGui::SliderFloat("##Interpolate QR",&t_QR,0,1))
        {
            mat3d Q,R;
            M.QR(Q,R);
            vec3d axis;
            float angle;
            matrix_to_axis_angle(Q,axis,angle); // interpolate Q as a pure rotation
            mat3d Q_int = mat3d::ROT_3D(axis,angle*t_QR);
            mat3d R_int = mat3d::DIAG(1.f)*(1-t_QR) + R*t_QR;
            cur = beg;
            cur.transform(Q_int*R_int);
            vec3d T_int = T*t_QR;
            cur.translate(T_int);
            cur.poly_set_color(Color::red_ramp_01(t_QR));
            cur.updateGL();
        }
    };

    return gui.launch();
}
