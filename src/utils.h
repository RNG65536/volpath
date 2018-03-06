#ifndef util_h__
#define util_h__

#include <string>

class vec3;

template<typename T>
T *rawPtr(std::vector<T>& x)
{
    return (&x[0]);
}

template<typename T>
const T *rawPtr(const std::vector<T>& x)
{
    return (&x[0]);
}

float f_min(float a, float b);
float f_max(float a, float b);
int clampi(int x, int a, int b);

float signed_map(int x, int n);
float sq(float x);

float rand01();

std::string str(const vec3& v);

class Frame
{
    vec3 n, t, b; // normal, tangent, bitangent

public:
    Frame(const vec3& normal);
    vec3 toWorld(const vec3& c) const;
    const vec3& normal() const;
    const vec3& tangent() const;
    const vec3& bitangent() const;
};

float misWeightPowerHeuristic(float pdf0, float pdf1);
float misWeightBalanceHeuristic(float pdf0, float pdf1);

#endif // util_h__
