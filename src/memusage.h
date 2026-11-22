#pragma once

// Process memory in MB (resident set / working set). NaN if the platform is not supported.
// Kept in its own translation unit because <windows.h> must not be included next to raylib.h.
double currentRssMB();   // right now
double peakRssMB();      // high-water mark since the process started
