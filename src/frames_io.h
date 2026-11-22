#pragma once

// Reading / writing of the per-iteration position files, shared by the solver (main.cpp) and the player (player.cpp).
//
//   <run>/faces.csv              header "v0,v1,v2", zero based vertex indices, one triangle per line
//   <run>/frames/pos_0000.csv    header "x,y,z", one vertex per line, same vertex order as the OBJ
//   <run>/frames/pos_0000.obj    (optional) same positions + faces, e.g. for ParaView

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "raylib.h"

inline std::string framePath(const std::string& runDir, int iter, const char* ext)
{
    char buf[64];

    std::snprintf(buf, sizeof(buf), "/frames/pos_%04d.%s", iter, ext);

    return runDir + buf;
}

inline bool writeFrameCsv(const std::string& path, const std::vector<Vector3>& P)
{
    FILE* f = std::fopen(path.c_str(), "wb");

    if (!f)
    {
        return false;
    }

    std::fputs("x,y,z\n", f);

    for (const Vector3& p : P)
    {
        std::fprintf(f, "%.9g,%.9g,%.9g\n", p.x, p.y, p.z);
    }

    std::fclose(f);

    return true;
}

inline bool writeFrameObj(const std::string& path, const std::vector<Vector3>& P, const std::vector<int>& faces)
{
    FILE* f = std::fopen(path.c_str(), "wb");

    if (!f)
    {
        return false;
    }

    for (const Vector3& p : P)
    {
        std::fprintf(f, "v %.9g %.9g %.9g\n", p.x, p.y, p.z);
    }

    for (size_t k = 0; k + 2 < faces.size(); k += 3)
    {
        std::fprintf(f, "f %d %d %d\n", faces[k] + 1, faces[k + 1] + 1, faces[k + 2] + 1);
    }

    std::fclose(f);

    return true;
}

inline bool writeFacesCsv(const std::string& path, const std::vector<int>& faces)
{
    FILE* f = std::fopen(path.c_str(), "wb");

    if (!f)
    {
        return false;
    }

    std::fputs("v0,v1,v2\n", f);

    for (size_t k = 0; k + 2 < faces.size(); k += 3)
    {
        std::fprintf(f, "%d,%d,%d\n", faces[k], faces[k + 1], faces[k + 2]);
    }

    std::fclose(f);

    return true;
}

// whole file into a string
inline bool readWholeFile(const std::string& path, std::string& out)
{
    FILE* f = std::fopen(path.c_str(), "rb");

    if (!f)
    {
        return false;
    }

    std::fseek(f, 0, SEEK_END);
    long size = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);

    out.resize((size_t)size);

    size_t got = size > 0 ? std::fread(&out[0], 1, (size_t)size, f) : 0;

    std::fclose(f);

    return got == (size_t)size;
}

inline bool readFrameCsv(const std::string& path, std::vector<Vector3>& P)
{
    std::string text;

    if (!readWholeFile(path, text))
    {
        return false;
    }

    P.clear();

    const char* c = text.c_str();
    const char* end = c + text.size();

    while (c < end && *c != '\n') ++c;          // header
    while (c < end && *c == '\n') ++c;

    while (c < end)
    {
        char* next = nullptr;

        float x = std::strtof(c, &next);
        if (next == c) break;
        c = next;
        if (c < end && *c == ',') ++c;

        float y = std::strtof(c, &next);
        c = next;
        if (c < end && *c == ',') ++c;

        float z = std::strtof(c, &next);
        c = next;

        P.push_back({ x, y, z });

        while (c < end && (*c == '\r' || *c == '\n')) ++c;
    }

    return !P.empty();
}

inline bool readFacesCsv(const std::string& path, std::vector<int>& faces)
{
    std::string text;

    if (!readWholeFile(path, text))
    {
        return false;
    }

    faces.clear();

    const char* c = text.c_str();
    const char* end = c + text.size();

    while (c < end && *c != '\n') ++c;
    while (c < end && *c == '\n') ++c;

    while (c < end)
    {
        char* next = nullptr;

        for (int k = 0; k < 3; ++k)
        {
            long v = std::strtol(c, &next, 10);
            if (next == c) return !faces.empty() && faces.size() % 3 == 0;
            faces.push_back((int)v);
            c = next;
            if (c < end && *c == ',') ++c;
        }

        while (c < end && (*c == '\r' || *c == '\n')) ++c;
    }

    return !faces.empty();
}
