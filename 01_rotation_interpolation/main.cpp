#undef NDEBUG
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// https://learnopencv.com/rotation-matrix-to-euler-angles/#:~:text=def%20isRotationMatrix(R)%20:%20Rt,(%5Bx%2C%20y%2C%20z%5D)
//
void matrix_to_Euler_angles(const mat3d & R, float Euler_angles[3])
{
    float sy       = sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0));
    bool  singular = sy<1e-6;

    if(!singular)
    {
        Euler_angles[0] = atan2( R(2,1), R(2,2));
        Euler_angles[1] = atan2(-R(2,0), sy    );
        Euler_angles[2] = atan2( R(1,0), R(0,0));
    }
    else
    {
        Euler_angles[0] = atan2(-R(1,2), R(1,1));
        Euler_angles[1] = atan2(-R(2,0), sy    );
        Euler_angles[2] = 0;
    }
}

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
    float t_Euler      = 0.0;
    float t_matrix     = 0.0;
    float t_matrix_fix = 0.0;
    float t_axis_angle = 0.0;
    mat3d R            = mat3d::DIAG(1.0);
    vec3d T            = vec3d(beg.bbox().delta_x()*2,0,0);
    vec3d c_end        = beg.bbox().center()+T;

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
            R = R_delta*R;
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
            R = R_delta*R;
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
            R = R_delta*R;
            end.translate(-c_end);
            end.transform(R_delta);
            end.translate(c_end);
            end.updateGL();
        }
        ImGui::NewLine();
        ImGui::Text("Interpolate Euler");
        if(ImGui::SliderFloat("##euler",&t_Euler,0.f,1.f))
        {
            float R_Euler_rad[3];
            matrix_to_Euler_angles(R, R_Euler_rad);
            mat3d Rx = mat3d::ROT_3D(vec3d(1,0,0),t_Euler*R_Euler_rad[0]);
            mat3d Ry = mat3d::ROT_3D(vec3d(0,1,0),t_Euler*R_Euler_rad[1]);
            mat3d Rz = mat3d::ROT_3D(vec3d(0,0,1),t_Euler*R_Euler_rad[2]);
            mat3d R_int = Rz*Ry*Rx;
            cur = beg;
            cur.transform(R_int);
            vec3d T_int = T*t_Euler;
            cur.translate(T_int);
            cur.poly_set_color(Color::red_ramp_01(t_Euler));
            cur.updateGL();
        }
        ImGui::Text("Interpolate Matrix");
        if(ImGui::SliderFloat("##matrix",&t_matrix,0.f,1.f))
        {
            mat3d R_int = mat3d::DIAG(1.f)*(1-t_matrix) + R*t_matrix;
            cur = beg;
            cur.transform(R_int);
            vec3d T_int = T*t_matrix;
            cur.translate(T_int);
            cur.poly_set_color(Color::red_ramp_01(t_matrix));
            cur.updateGL();
        }
        ImGui::Text("Interpolate Matrix (corrected)");
        if(ImGui::SliderFloat("##matrix2",&t_matrix_fix,0.f,1.f))
        {
            mat3d R_int = mat3d::DIAG(1.f)*(1-t_matrix_fix) + R*t_matrix_fix;
            mat3d U,V;
            vec3d S;
            R_int.SVD(U,S,V);
            mat3d I = mat3d::DIAG(1);
            I(2,2) = (U*V.transpose()).det();
            R_int = U*I*V.transpose();
            cur = beg;
            cur.transform(R_int);
            vec3d T_int = T*t_matrix_fix;
            cur.translate(T_int);
            cur.poly_set_color(Color::red_ramp_01(t_matrix_fix));
            cur.updateGL();
        }
        ImGui::Text("Interpolate Axid/Angle");
        if(ImGui::SliderFloat("##axisangle",&t_axis_angle,0.f,1.f))
        {
            vec3d axis;
            float angle;
            matrix_to_axis_angle(R,axis,angle);
            mat3d R_int = mat3d::ROT_3D(axis,angle*t_axis_angle);
            cur = beg;
            cur.transform(R_int);
            vec3d T_int = T*t_axis_angle;
            cur.translate(T_int);
            cur.poly_set_color(Color::red_ramp_01(t_axis_angle));
            cur.updateGL();
        }
    };

    return gui.launch();
}
