#ifndef skydome_h__
#define skydome_h__

#include "hosek/ArHosekSkyModel.h"

vec3 toDirection(float phi, float theta)
{
    return vec3(
        sin(theta) * sin(phi),
        cos(theta),
        sin(theta) * -cos(phi));
}

std::pair<float, float> toAngle(const vec3& direction)
{
    float theta = acos(direction.y);
    float phi = atan(direction.z / direction.x) + M_PI * 0.5f;
    if (direction.x < 0)
    {
        phi += M_PI;
    }
    return std::make_pair(phi, theta);
}
class Skydome
{
    float m_solarElevation;
    float m_solarAzimuth;
    float m_brightness;
    vec3 m_ground_albedo;
    std::unique_ptr<ArHosekSkyModelState> skymodel_state[3];

public:
    Skydome(float turbidity, const vec3& ground_albedo, float solar_elevation, float solar_azimuth, float brightness)
    {
        for (int i = 0; i < 3; i++)
        {
            skymodel_state[i].reset(arhosek_rgb_skymodelstate_alloc_init(turbidity, ground_albedo[i], solar_elevation));
        }
        m_solarElevation = f_max(0.0f, f_min(M_PI_2, solar_elevation));
        m_solarAzimuth = solar_azimuth;
        m_brightness = brightness;
        m_ground_albedo = ground_albedo;
    }
    ~Skydome()
    {
    }

    vec3 sunDirection() const
    {
        return toDirection(m_solarAzimuth, M_PI_2 - m_solarElevation);
    }

    vec3 radiance(const vec3& dir) const
    {
        auto angles = toAngle(dir);
        float phi = angles.first;
        float theta = angles.second;
        vec3 sundir = toDirection(m_solarAzimuth, M_PI_2 - m_solarElevation);
        vec3 viewdir = toDirection(phi, theta);
        float gamma = acos(dot(viewdir, sundir));

        vec3 rad;
        if (theta < M_PI_2)
        {
            for (int i = 0; i < 3; i++)
            {
                rad[i] = (float)arhosek_tristim_skymodel_radiance(skymodel_state[i].get(), theta, gamma, i) * m_brightness;
            }
        }
        else
        {
            for (int i = 0; i < 3; i++)
            {
                rad[i] = m_ground_albedo[i];
            }

//             vec3 viewdir = toDirection(phi, M_PI - theta);
//             float gamma = acos(dot(viewdir, sundir));
//             for (int i = 0; i < 3; i++)
//             {
//                 rad[i] = (float)arhosek_tristim_skymodel_radiance(skymodel_state[i].get(), M_PI - theta, gamma, i);
//             }
//             rad *= m_ground_albedo;
        }

        return rad;
    }
};

#endif // skydome_h__
