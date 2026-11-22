// Plays back the per-iteration positions written by the solver (runs/<run>/frames/pos_XXXX.csv) as an animation.
//
//   player <run folder> [stride]
//
//   Space  play / pause          Left / Right   previous / next frame (while paused)
//   Up / Down  faster / slower   Home or R      back to the first frame
//   L  loop on/off               W  wireframe    S  screenshot (png, next to the executable)
//   mouse: left drag = orbit, wheel = zoom, right/middle drag = pan
//
// stride: only every stride-th frame is loaded (saves memory for long runs).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include "frames_io.h"

namespace fs = std::filesystem;

struct MetricsTable
{
    std::vector<std::string> names;
    std::map<int, std::vector<double>> rows;       //iteration -> values

    bool load(const std::string& path)
    {
        std::ifstream f(path);

        if (!f)
        {
            return false;
        }

        std::string line;
        std::getline(f, line);

        auto split = [](const std::string& s)
        {
            std::vector<std::string> out;
            std::string cell;
            std::stringstream ss(s);

            while (std::getline(ss, cell, ','))
            {
                out.push_back(cell);
            }

            if (!s.empty() && s.back() == ',')
            {
                out.push_back("");
            }

            return out;
        };

        names = split(line);

        while (names.size() && !names.back().empty() && (names.back().back() == '\r'))
        {
            names.back().pop_back();
        }

        while (std::getline(f, line))
        {
            std::vector<std::string> cells = split(line);

            if (cells.empty())
            {
                continue;
            }

            std::vector<double> vals(names.size(), std::nan(""));

            for (size_t k = 0; k < cells.size() && k < names.size(); ++k)
            {
                if (!cells[k].empty())
                {
                    vals[k] = std::atof(cells[k].c_str());
                }
            }

            rows[(int)vals[0]] = vals;
        }

        return true;
    }

    double get(int iter, const std::string& name) const
    {
        auto it = rows.find(iter);

        if (it == rows.end())
        {
            return std::nan("");
        }

        for (size_t k = 0; k < names.size(); ++k)
        {
            if (names[k] == name)
            {
                return it->second[k];
            }
        }

        return std::nan("");
    }
};

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cout << "usage: player <run folder> [stride]\n";

        return 1;
    }

    const std::string dir = argv[1];
    const int stride = argc > 2 ? std::max(1, std::atoi(argv[2])) : 1;

    //topology
    std::vector<int> faces;

    if (!readFacesCsv(dir + "/faces.csv", faces))
    {
        std::cout << "cannot read " << dir << "/faces.csv\n";

        return 1;
    }

    const int T = (int)(faces.size() / 3);

    //frame list: pos_<iter>.csv sorted by iter
    std::vector<int> iters;

    for (const auto& e : fs::directory_iterator(dir + "/frames"))
    {
        const std::string name = e.path().filename().string();

        if (name.rfind("pos_", 0) == 0 && e.path().extension() == ".csv")
        {
            iters.push_back(std::atoi(name.c_str() + 4));
        }
    }

    std::sort(iters.begin(), iters.end());

    std::vector<int> used;

    for (size_t k = 0; k < iters.size(); ++k)
    {
        if (k % stride == 0 || k + 1 == iters.size())
        {
            used.push_back(iters[k]);
        }
    }

    if (used.empty())
    {
        std::cout << "no frames in " << dir << "/frames\n";

        return 1;
    }

    //load everything
    std::vector<std::vector<Vector3>> frames;
    frames.reserve(used.size());

    Vector3 bmin = { 1e30f, 1e30f, 1e30f }, bmax = { -1e30f, -1e30f, -1e30f };

    for (size_t k = 0; k < used.size(); ++k)
    {
        std::vector<Vector3> P;

        if (!readFrameCsv(framePath(dir, used[k], "csv"), P))
        {
            std::cout << "cannot read frame " << used[k] << "\n";

            return 1;
        }

        if (!frames.empty() && P.size() != frames[0].size())
        {
            std::cout << "frame " << used[k] << " has a different vertex count\n";

            return 1;
        }

        for (const Vector3& p : P)
        {
            bmin = Vector3Min(bmin, p);
            bmax = Vector3Max(bmax, p);
        }

        frames.push_back(std::move(P));

        if (k % 20 == 0)
        {
            std::cout << "loading " << k + 1 << " / " << used.size() << "\r" << std::flush;
        }
    }

    std::cout << "loaded " << frames.size() << " frames, " << frames[0].size() << " vertices, " << T << " triangles\n";

    MetricsTable metrics;
    const bool haveMetrics = metrics.load(dir + "/metrics.csv");

    //everything is drawn at unit size around the origin (the raylib near plane would clip a tiny model)
    const Vector3 center = Vector3Scale(Vector3Add(bmin, bmax), 0.5f);
    const float diag = std::max(Vector3Distance(bmin, bmax), 1e-12f);
    const float modelScale = 1.0f / diag;
    const Vector3 modelPos = Vector3Scale(center, -modelScale);

    InitWindow(1200, 800, "WoScMCF player");
    SetTargetFPS(60);

    Camera3D camera = { 0 };
    camera.up = { 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    float yaw = 0.8f, pitch = 0.4f, dist = 2.2f;
    Vector3 target = { 0.0f, 0.0f, 0.0f };

    //unindexed mesh, vertices and colors are re-uploaded every frame
    Mesh mesh = { 0 };
    mesh.triangleCount = T;
    mesh.vertexCount = 3 * T;
    mesh.vertices = (float*)MemAlloc(mesh.vertexCount * 3 * sizeof(float));
    mesh.colors = (unsigned char*)MemAlloc(mesh.vertexCount * 4);

    for (int i = 0; i < mesh.vertexCount; ++i)
    {
        mesh.colors[4 * i + 3] = 255;
    }

    UploadMesh(&mesh, true);

    Model model = LoadModelFromMesh(mesh);

    auto fillMesh = [&](const std::vector<Vector3>& P)
    {
        Mesh& m = model.meshes[0];
        const Vector3 toCam = Vector3Normalize(Vector3Subtract(camera.position, camera.target));

        for (int t = 0; t < T; ++t)
        {
            const Vector3 a = P[faces[3 * t]], b = P[faces[3 * t + 1]], c = P[faces[3 * t + 2]];
            const Vector3 n = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
            const float shade = 0.3f + 0.7f * std::fabs(Vector3DotProduct(n, toCam));    //two sided headlight

            const Vector3 v[3] = { a, b, c };

            for (int k = 0; k < 3; ++k)
            {
                const int i = 3 * t + k;

                m.vertices[3 * i + 0] = v[k].x;
                m.vertices[3 * i + 1] = v[k].y;
                m.vertices[3 * i + 2] = v[k].z;

                m.colors[4 * i + 0] = (unsigned char)(120.0f * shade);
                m.colors[4 * i + 1] = (unsigned char)(170.0f * shade);
                m.colors[4 * i + 2] = (unsigned char)(230.0f * shade);
            }
        }

        UpdateMeshBuffer(m, 0, m.vertices, m.vertexCount * 3 * sizeof(float), 0);
        UpdateMeshBuffer(m, 3, m.colors, m.vertexCount * 4, 0);
    };

    int cur = 0;
    bool playing = true, loop = true, wires = false;
    float fps = 15.0f, acc = 0.0f;
    int shot = 0;

    while (!WindowShouldClose())
    {
        //camera
        const Vector2 d = GetMouseDelta();

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
        {
            yaw -= d.x * 0.005f;
            pitch = std::clamp(pitch + d.y * 0.005f, -1.5f, 1.5f);
        }

        dist = std::clamp(dist * (1.0f - GetMouseWheelMove() * 0.1f), 0.2f, 50.0f);

        camera.position = Vector3Add(target, { dist * cosf(pitch) * sinf(yaw), dist * sinf(pitch), dist * cosf(pitch) * cosf(yaw) });
        camera.target = target;

        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE))
        {
            const Vector3 fwd = Vector3Normalize(Vector3Subtract(camera.target, camera.position));
            const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, camera.up));
            const Vector3 up = Vector3CrossProduct(right, fwd);
            const float s = dist * 0.0015f;

            target = Vector3Add(target, Vector3Scale(right, -d.x * s));
            target = Vector3Add(target, Vector3Scale(up, d.y * s));
        }

        //keys
        if (IsKeyPressed(KEY_SPACE)) playing = !playing;
        if (IsKeyPressed(KEY_L)) loop = !loop;
        if (IsKeyPressed(KEY_W)) wires = !wires;
        if (IsKeyPressed(KEY_UP)) fps = std::min(fps * 2.0f, 480.0f);
        if (IsKeyPressed(KEY_DOWN)) fps = std::max(fps * 0.5f, 0.5f);
        if (IsKeyPressed(KEY_HOME) || IsKeyPressed(KEY_R)) cur = 0;

        if (IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) { playing = false; cur = std::min(cur + 1, (int)frames.size() - 1); }
        if (IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) { playing = false; cur = std::max(cur - 1, 0); }

        if (playing)
        {
            acc += GetFrameTime() * fps;

            while (acc >= 1.0f)
            {
                acc -= 1.0f;

                if (cur + 1 < (int)frames.size())
                {
                    ++cur;
                }
                else if (loop)
                {
                    cur = 0;
                }
                else
                {
                    playing = false;
                }
            }
        }

        fillMesh(frames[cur]);

        BeginDrawing();

        ClearBackground(WHITE);

        BeginMode3D(camera);

        rlDisableBackfaceCulling();

        if (wires)
        {
            DrawModelWires(model, modelPos, modelScale, DARKGRAY);
        }
        else
        {
            DrawModel(model, modelPos, modelScale, WHITE);
        }

        rlEnableBackfaceCulling();

        EndMode3D();

        const int it = used[cur];

        DrawText(TextFormat("iteration %d   (frame %d / %d)   %s   %.1f frames/s   loop %s", it, cur + 1, (int)frames.size(), playing ? "playing" : "paused", fps, loop ? "on" : "off"), 20, 20, 20, DARKGRAY);

        if (haveMetrics)
        {
            DrawText(TextFormat("sphericity = %.4f   sigma_bar = %.4g   min angle = %.2f deg   min area ratio = %.3g",
                                metrics.get(it, "sphericity"), metrics.get(it, "sigma_bar"), metrics.get(it, "min_angle_deg"), metrics.get(it, "min_area_ratio")), 20, 46, 18, DARKGRAY);
        }

        DrawText("Space: play/pause  Left/Right: step  Up/Down: speed  Home/R: restart  L: loop  W: wires  S: screenshot", 20, 770, 18, GRAY);

        EndDrawing();

        if (IsKeyPressed(KEY_S))
        {
            TakeScreenshot(TextFormat("frame_%04d_%d.png", it, shot++));
        }
    }

    UnloadModel(model);
    CloseWindow();

    return 0;
}
