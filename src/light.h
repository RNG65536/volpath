#ifndef light_h__
#define light_h__

class Light
{
    // sphere light
    vec3 center = vec3(2, 2, 1);
//     vec3 center = vec3(-2, -2, -4);
    float radius = 0.5f;
    float emission = 10.0f / (radius * radius);

public:
    void sample(vec3 position, float& pdfW, vec3& light_pos, vec3& light_dir, vec3& Li) const
    {
        Frame frame(center - position);

        float d2 = dot(position - center, position - center);
        float l2 = d2 - radius * radius;
        float cos_a_max = std::sqrt(l2 / d2);

        float rnd0 = rand01(), rnd1 = rand01();
        float cos_a = 1 - rnd0 + rnd0 * cos_a_max;
        float sin_a = std::sqrt(1 - cos_a * cos_a);
        float phi = 2 * M_PI * rnd1;
        light_dir = frame.toWorld(vec3(std::cos(phi) * sin_a, std::sin(phi) * sin_a, cos_a));
        light_dir = normalize(light_dir);

        pdfW = 1.0f / (2.0f * M_PI * (1.0f - cos_a_max));
        Li = vec3(emission);
        light_pos = position + light_dir * std::sqrt(l2);
    }
    vec3 Li(const vec3& position, const vec3& direction, vec3& light_pos) const
    {
        float d2 = dot(position - center, position - center);
        float l2 = d2 - radius * radius;
        if (l2 < 0) // inside sphere
        {
            return vec3(0.0f);
        }
        float cos_a_max = std::sqrt(l2 / d2);
        float cos_a = dot(direction, normalize(center - position));
        if (cos_a < cos_a_max)
        {
            return vec3(0.0f);
        }
        light_pos = position + direction * std::sqrt(l2);
        return vec3(emission); // todo : add intersection
    }
    float pdf(vec3& position, const vec3& direction, vec3& light_pos, vec3& Li) const
    {
        Li = vec3(0.0f);
        float d2 = dot(position - center, position - center);
        float l2 = d2 - radius * radius;
        if (l2 < 0) // inside sphere
        {
            return 0;
        }
        float cos_a_max = std::sqrt(l2 / d2);
        float cos_a = dot(direction, normalize(center - position));
        if (cos_a < cos_a_max)
        {
            return 0;
        }
        light_pos = position + direction * std::sqrt(l2);
        Li = vec3(emission);
        return 1.0f / (2.0f * M_PI * (1.0f - cos_a_max));
    }
};


#endif // light_h__
