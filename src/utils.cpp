#include <random>
#include <sstream>
#include "vector.h"
#include "utils.h"

static std::default_random_engine rng;
static std::uniform_real_distribution<float> dist01(0.0f, 1.0f);

float f_min(float a, float b)
{
    return a < b ? a : b;
}

float f_max(float a, float b)
{
    return a > b ? a : b;
}

float signed_map(int x, int n)
{
    return 2 * (x / (float)n) - 1;
}

float sq(float x)
{
    return x*x;
}

float rand01()
{
    return dist01(rng);
}

std::string str(const vec3& v)
{
    std::stringstream ss;
    ss << "( " << v.x << ", " << v.y << ", " << v.z << " )";
    return ss.str();
}

int clampi(int x, int a, int b)
{
    return x < a ? a : x > b ? b : x;
}

float misWeightPowerHeuristic(float pdf0, float pdf1)
{
    float a = pdf0 * pdf0;
    return a / (a + pdf1 * pdf1);
}

float misWeightBalanceHeuristic(float pdf0, float pdf1)
{
    return pdf0 / (pdf0 + pdf1);
}

vec3 Frame::toWorld(const vec3& c) const
{
    return t * c.x + b * c.y + n * c.z;
}

Frame::Frame(const vec3& normal)
{
    n = normalize(normal);
    vec3 a = fabs(n.x) > 0.1 ? vec3(0, 1, 0) : vec3(1, 0, 0);
    t = normalize(cross(a, n));
    b = cross(n, t);
}

const vec3& Frame::bitangent() const
{
    return b;
}

const vec3& Frame::tangent() const
{
    return t;
}

const vec3& Frame::normal() const
{
    return n;
}
