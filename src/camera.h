class Camera
{
    vec3 cam_o, cam_d, cam_x, cam_y;
    float fov;
    int width, height;

public:
    Camera(const vec3& origin, const vec3& lookat, float fovy, int width, int height)
        : cam_o(origin), fov(fovy), width(width), height(height)
    {
        vec3 up(0, 1, 0);
        cam_d = vec3(normalize(lookat - origin));
        cam_x = vec3(normalize(cross(cam_d, up)));
        cam_y = vec3(cross(cam_x, cam_d));
    }
    Ray pixelRay(int i, int j)
    {
        Ray r(cam_o, cam_d
            + cam_x * (signed_map(i, width) * tan(fov * 0.5f / 180.0f * M_PI))
            + cam_y * (signed_map(j, height) * tan(fov * 0.5f / 180.0f * M_PI)));
        return r;
    }
};
