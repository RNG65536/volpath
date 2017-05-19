#ifndef light_h__
#define light_h__

class Light
{
    // sphere light
    vec3 center = vec3(2, 2, 1);
    float radius = 0.5f;
    float emission = 50.0f / (radius * radius);

public:
    void sample(vec3 position, float& pdfW, vec3& light_pos, vec3& light_dir, vec3& Li)
    {
        Frame frame(center - position);

        float d2 = dot(position - center, position - center);
        float l2 = d2 - radius * radius;
        float cos_a_max = std::sqrt(l2 / d2);
        float eps1 = rand01(), eps2 = rand01();
        float cos_a = 1 - eps1 + eps1 * cos_a_max;
        float sin_a = std::sqrt(1 - cos_a * cos_a);
        float phi = 2 * M_PI * eps2;

        light_dir = frame.toWorld(vec3(std::cos(phi) * sin_a, std::sin(phi) * sin_a, cos_a));
        light_dir = normalize(light_dir);
        pdfW = 1.0f / (2.0f * M_PI * (1.0f - cos_a_max));
        Li = vec3(emission);
        light_pos = position + light_dir * std::sqrt(l2);
    }
    vec3 Li(const vec3& position, const vec3& direction) const
    {
        return vec3(0.0f); // todo : add intersection
    }
};


#endif // light_h__
