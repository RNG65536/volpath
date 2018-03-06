#ifndef timer_h__
#define timer_h__

#ifdef _WIN32

class TimeData;

class Timer
{
    TimeData *data;

public:
    Timer();
    ~Timer();

    // Get elapsed time from last reset()
    // or class construction.
    // return the elapsed time
    float elapsed();

    // reset the timer
    void reset();
};

#endif

#endif // timer_h__
