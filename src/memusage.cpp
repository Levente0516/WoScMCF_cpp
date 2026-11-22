#include "memusage.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#if defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <psapi.h>

double currentRssMB()
{
    PROCESS_MEMORY_COUNTERS pmc;

    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        return (double)pmc.WorkingSetSize / (1024.0 * 1024.0);
    }

    return std::nan("");
}

double peakRssMB()
{
    PROCESS_MEMORY_COUNTERS pmc;

    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        return (double)pmc.PeakWorkingSetSize / (1024.0 * 1024.0);
    }

    return std::nan("");
}

#elif defined(__linux__)

static double readStatusKb(const char* key)
{
    FILE* f = std::fopen("/proc/self/status", "r");

    if (!f)
    {
        return std::nan("");
    }

    char line[256];
    double kb = std::nan("");
    const size_t n = std::strlen(key);

    while (std::fgets(line, sizeof(line), f))
    {
        if (std::strncmp(line, key, n) == 0)
        {
            kb = std::atof(line + n);
            break;
        }
    }

    std::fclose(f);

    return kb;
}

double currentRssMB() { return readStatusKb("VmRSS:") / 1024.0; }
double peakRssMB()    { return readStatusKb("VmHWM:") / 1024.0; }

#else

double currentRssMB() { return std::nan(""); }
double peakRssMB()    { return std::nan(""); }

#endif
