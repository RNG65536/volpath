#ifndef scene_h__
#define scene_h__

bool intersect_vol(const Ray& r, const vec3& aabb_min, const vec3& aabb_max, float& t_near, float& t_far)
{
    vec3 bounds[2] = { aabb_min, aabb_max };
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    tmin = (bounds[r.sign[0]].x - r.o.x) * r.invdir.x;
    tmax = (bounds[1 - r.sign[0]].x - r.o.x) * r.invdir.x;
    tymin = (bounds[r.sign[1]].y - r.o.y) * r.invdir.y;
    tymax = (bounds[1 - r.sign[1]].y - r.o.y) * r.invdir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
    tzmin = (bounds[r.sign[2]].z - r.o.z) * r.invdir.z;
    tzmax = (bounds[1 - r.sign[2]].z - r.o.z) * r.invdir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    if (tmin > 0) t_near = tmin; else t_near = 0;
    if (tmax < inf) t_far = tmax;
    if (tmax >= inf){
        return false;
    }
    return true;
}

#endif // scene_h__
