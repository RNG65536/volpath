#include <cmath>
#include "vector.hpp"

vec3 normalize(const vec3& a)
{
    float len = sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
    return a * (1.0f / len);
}

float dot(const vec3& a, const vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 cross(const vec3& a, const vec3& b)
{
    return vec3(a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z, 
                a.x * b.y - a.y * b.x);
}

float distanceSquared(const vec3& a, const vec3& b)
{
    vec3 d(a - b);
    return dot(d, d);
}

float distance(const vec3& a, const vec3& b)
{
    return sqrtf(distanceSquared(a, b));
}

float length(const vec3& a)
{
    return sqrtf(dot(a, a));
}


vec3& vec3::operator/=(float b)
{
    return *this = *this / b;
}

vec3 vec3::operator/=(const vec3& b)
{
    return *this = *this / b;
}

vec3& vec3::operator*=(float b)
{
    return *this = *this * b;
}

vec3 vec3::operator*=(const vec3& b)
{
    return *this = *this * b;
}

vec3& vec3::operator-=(float b)
{
    return *this = *this - b;
}

vec3& vec3::operator+=(float b)
{
    return *this = *this + b;
}

vec3& vec3::operator+=(const vec3& b)
{
    return *this = *this + b;
}

vec3 vec3::operator/(float b) const
{
    b = 1.0f / b;
    return vec3(x * b, y * b, z * b);
}

vec3 vec3::operator/(const vec3& b) const
{
    return vec3(x / b.x, y / b.y, z / b.z);
}

vec3 vec3::operator*(float b) const
{
    return vec3(x * b, y * b, z * b);
}

vec3 vec3::operator*(const vec3& b) const
{
    return vec3(x * b.x, y * b.y, z * b.z);
}

vec3 vec3::operator-(float b) const
{
    return vec3(x - b, y - b, z - b);
}

vec3 vec3::operator+(float b) const
{
    return vec3(x + b, y + b, z + b);
}

vec3 vec3::operator+(const vec3& b) const
{
    return vec3(x + b.x, y + b.y, z + b.z);
}

vec3 vec3::operator-=(const vec3& b)
{
    return *this = *this - b;
}

vec3 vec3::operator-(const vec3& b) const
{
    return vec3(x - b.x, y - b.y, z - b.z);
}

vec3 vec3::operator-() const
{
    return vec3(-x, -y, -z);
}

const float& vec3::operator[](int n) const
{
    return (&x)[n];
}

float& vec3::operator[](int n)
{
    return (&x)[n];
}

vec3::vec3(float x_, float y_, float z_) :
x(x_), y(y_), z(z_)
{

}

vec3::vec3(float a_) : x(a_), y(a_), z(a_)
{

}

vec3::vec3() : x(0), y(0), z(0)
{

}

bool vec3::operator!=(const vec3& b)
{
    return x != b.x || y != b.y || z != b.z;
}
