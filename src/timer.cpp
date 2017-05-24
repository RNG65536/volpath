#include <Windows.h>
#include "timer.h"

class TimeData
{
public:
    LARGE_INTEGER _freq, _start, _stop;
};

Timer::~Timer()
{
    delete data;
}

void Timer::reset()
{
    QueryPerformanceCounter(&data->_start);
    data->_stop = data->_start;
}

float Timer::elapsed()
{
    QueryPerformanceCounter(&data->_stop);
    return float(data->_stop.QuadPart - data->_start.QuadPart)
        / float(data->_freq.QuadPart);
}

Timer::Timer()
{
    data = new TimeData;
    QueryPerformanceFrequency(&data->_freq);
    QueryPerformanceCounter(&data->_start);
    data->_stop = data->_start;
}
