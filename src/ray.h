#ifndef ray_h__
#define ray_h__

struct Ray
{
    vec3 o, d;
    vec3 invdir;
    int sign[3];
    Ray(){}
    Ray(const vec3& o_, const vec3& d_) :o(o_), d(normalize(d_))
    {
        invdir.x = 1.0f / d.x;
        invdir.y = 1.0f / d.y;
        invdir.z = 1.0f / d.z;
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
    }
    vec3 at(float t) const { return o + d * t; }
};

#endif // ray_h__
