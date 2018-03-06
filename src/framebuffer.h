#ifndef framebuffer_h__
#define framebuffer_h__

class FrameBuffer
{
    int m_width, m_height;
    std::vector<vec3> m_buffer;

public:
    FrameBuffer(int width, int height);
    void scaleBrightness(float s);
    std::vector<vec3>& buffer();
    void dump() const;
};

#endif // framebuffer_h__
