#ifndef vector_h__
#define vector_h__

class vec3
{
public:
    float x, y, z;
//     inline const float& x() const { return x; }
//     inline const float& y() const { return y; }
//     inline const float& z() const { return z; }
//     inline float& x() { return x; }
//     inline float& y() { return y; }
//     inline float& z() { return z; }

    vec3();
    explicit vec3(float a_);
    vec3(float x_, float y_, float z_);
    float& operator[](int n);
    const float& operator[](int n) const;

    vec3 operator-() const;
    vec3 operator+ (const vec3& b) const;
    vec3 operator- (const vec3& b) const;
    vec3 operator* (const vec3& b) const;
    vec3 operator/ (const vec3& b) const;
    vec3& operator+= (const vec3& b);
    vec3 operator-= (const vec3& b);
    vec3 operator*= (const vec3& b);
    vec3 operator/= (const vec3& b);
    vec3 operator+ (float b) const;
    vec3 operator- (float b) const;
    vec3 operator* (float b) const;
    vec3 operator/ (float b) const;
    vec3& operator+= (float b);
    vec3& operator-= (float b);
    vec3& operator*= (float b);
    vec3& operator/= (float b);
    bool operator!=(const vec3& b);
};

vec3 normalize(const vec3& a);
float dot(const vec3& a, const vec3& b);
vec3 cross(const vec3& a, const vec3& b);
float distanceSquared(const vec3& a, const vec3& b);
float distance(const vec3& a, const vec3& b);
float length(const vec3& a);

#endif // vector_h__
