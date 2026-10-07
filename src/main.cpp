#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "raylib.h"
#include "raymath.h"
#include <random>
#include "mantis.h"

struct Vec3D{
    float x, y, z;
};

struct Vertex{
    Vec3D vert;
};

struct Triangle{
    int v0, v1, v2;
};

struct Surface{
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
};

struct OrbitCam {
    float yaw = 0.8f, pitch = 0.5f, dist = 4.5f;
    Vector3 target = { 0, 0, 0 };
};

Surface loadOBJ(const std::string& filename)
{
    Surface surface;

    std::ifstream file(filename);

    if (!file) { throw std::runtime_error("Could not open file"); }

    std::string line;

    while(std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string type;

        ss >> type;

        if (type == "v")
        {
            Vertex v;
            ss >> v.vert.x >> v.vert.y >> v.vert.z;
            surface.vertices.push_back(v);
        }
        else if (type == "f")
        {
            std::string s0, s1, s2;
            ss >> s0 >> s1 >> s2;

            auto getVertexIndex = [](const std::string& s)
            {
                size_t slash = s.find('/');
                return std::stoi(s.substr(0, slash)) - 1;
            };

            Triangle t;

            t.v0 = getVertexIndex(s0);
            t.v1 = getVertexIndex(s1);
            t.v2 = getVertexIndex(s2);

            surface.triangles.push_back(t);
        }
    }

    return surface;
}


std::vector<Vec3D> samplePointsInsideSphere(Vec3D center, float r, int Ns)
{
    std::vector<Vec3D> randomPoints;

    for (int i = 0; i < Ns; i++)
    {
        Vec3D cords;

        float u = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float v = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float phi = 2 * PI * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        float radius;
        float theta = acos(2*v - 1);

        if (Ns > 1)
        {
            radius = r * pow(u, (1.0f/3.0f));
        }
        else if (Ns == 1)
        {
            radius = r;
        }

        cords.x = center.x + radius * sin(theta) * cos(phi);
        cords.y = center.y + radius * sin(theta) * sin(phi);
        cords.z = center.z + radius * cos(theta);

        randomPoints.push_back(cords);
    }

    return randomPoints;
}

Vec3D closestPointOnSurface(const mantis::AccelerationStructure& accel, Vec3D q)
{
    mantis::Result r = accel.calc_closest_point(q.x, q.y, q.z);

    return { r.closest_point[0], r.closest_point[1], r.closest_point[2] };
}

std::tuple<Vec3D, float> MaEtAl(Vec3D n, float R, Vec3D p, const mantis::AccelerationStructure& accel)
{
    const float tol = 1e-4f * R;
    const int maxIter = 1000;                      

    float r = R;
    Vec3D c = { p.x + r*n.x, p.y + r*n.y, p.z + r*n.z };
    Vec3D q = closestPointOnSurface(accel, c);
    Vec3D q_prev = { 1e30f, 1e30f, 1e30f };         

    auto dist = [](Vec3D a, Vec3D b) {
        float dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };

    int it = 0;

    while (dist(q, q_prev) > tol && it++ < maxIter)       
    {
        if (dist(q, p) < tol) 
        {
            break;                      
        }

        float qx = q.x - p.x;
        float qy = q.y - p.y;
        float qz = q.z - p.z;

        float denom = 2.0f * (n.x * qx + n.y * qy + n.z * qz);

        if (denom <= 1e-12f) 
        {
            break;
        }

        r = (qx*qx + qy*qy + qz*qz) / denom;             
        c = { p.x + r*n.x, p.y + r*n.y, p.z + r*n.z };
        q_prev = q;
        q = closestPointOnSurface(accel, c);
    }

    return { c, r };
}

std::tuple<std::vector<Vec3D>, Vec3D, float> MedialAxisPointCloud(Surface surface, int Ns, const mantis::AccelerationStructure& accel)
{
    //Get the bounding box, calculate the half length of its diagonal and take a sphere of a radius with that
    Vec3D bb_min = surface.vertices[0].vert;
    Vec3D bb_max = surface.vertices[0].vert;

    for (const Vertex& v : surface.vertices)
    {
        bb_min.x = std::min(bb_min.x, v.vert.x);
        bb_min.y = std::min(bb_min.y, v.vert.y);
        bb_min.z = std::min(bb_min.z, v.vert.z);

        bb_max.x = std::max(bb_max.x, v.vert.x);
        bb_max.y = std::max(bb_max.y, v.vert.y);
        bb_max.z = std::max(bb_max.z, v.vert.z);
    }

    Vec3D center = {(bb_min.x + bb_max.x) / 2, 
                    (bb_min.y + bb_max.y) / 2,
                    (bb_min.z + bb_max.z) / 2};

    float bb_range_x = bb_max.x -bb_min.x;
    float bb_range_y = bb_max.y -bb_min.y;
    float bb_range_z = bb_max.z -bb_min.z;

    float r = std::sqrt(std::pow(bb_range_x,2) + std::pow(bb_range_y,2) + std::pow(bb_range_z,2)) / 2;

    std::vector<Vec3D> medialBallPoints = samplePointsInsideSphere(center, r, Ns);

    std::vector<Vec3D> medialAxisPoints;

    std::vector<float> r_new;

    for (const Vec3D points : medialBallPoints)
    {
        Vec3D normal;
        Vec3D normalInv;
        Vec3D p = closestPointOnSurface(accel, points);

        float x = points.x - p.x;
        float y = points.y - p.y;
        float z = points.z - p.z;

        float magnitude = sqrt(x*x + y*y + z*z);

        normal.x = x/magnitude;
        normal.y = y/magnitude;
        normal.z = z/magnitude;

        normalInv.x = -1 * normal.x;
        normalInv.y = -1 * normal.y;
        normalInv.z = -1 * normal.z;

        auto [medialAxisPointNormal, normalR] = MaEtAl(normal, r*2, p, accel);
        auto [medialAxisPointNormalInv, normalInvR] = MaEtAl(normalInv, r*2, p, accel);


        float min_r = std::min(normalR, normalInvR);

        if (normalR > normalInvR)
        {
            medialAxisPointNormal = {
                p.x + normalInvR * normal.x,
                p.y + normalInvR * normal.y,
                p.z + normalInvR * normal.z
            };
        }
        else if (normalInvR > normalR)
        {
            medialAxisPointNormalInv = {
                p.x + normalR * normalInv.x,
                p.y + normalR * normalInv.y,
                p.z + normalR * normalInv.z
            };
        }

        r_new.push_back(min_r);
        r_new.push_back(min_r);
        medialAxisPoints.push_back(medialAxisPointNormal);
        medialAxisPoints.push_back(medialAxisPointNormalInv);
    }

    float s = 1.15f;

    std::vector<bool> removed(medialAxisPoints.size(), false);

    for (size_t i = 0; i < medialAxisPoints.size(); ++i)
    {
        for (size_t j = i + 1; j < medialAxisPoints.size(); ++j)
        {
            float dx = medialAxisPoints[i].x - medialAxisPoints[j].x;
            float dy = medialAxisPoints[i].y - medialAxisPoints[j].y;
            float dz = medialAxisPoints[i].z - medialAxisPoints[j].z;

            float distance = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r_new[i] < r_new[j] && distance + s * r_new[i] < s * r_new[j])
            {
                removed[i] = true;
            }     
            else if (r_new[j] < r_new[i] && distance + s * r_new[j] < s * r_new[i])
            {
                removed[j] = true;
            } 
        }
    }

    std::vector<Vec3D> prunedPoints;

    for (size_t i = 0; i < medialAxisPoints.size(); ++i)
    {
        if (!removed[i])
        {
            prunedPoints.push_back(medialAxisPoints[i]);
        }
    }

    return {prunedPoints, center, r};
}

float localFetureSize(const std::vector<Vec3D>& M, Vec3D x)
{
    float minDistance = std::numeric_limits<float>::max();

    for (const Vec3D& m : M)
    {
        float dx = x.x - m.x;
        float dy = x.y - m.y;
        float dz = x.z - m.z;
        
        float d2 = dx * dx + dy * dy + dz * dz;

        minDistance = std::min(minDistance, d2);
    }

    float distance = sqrt(minDistance);

    return (std::max(0.9f * distance, 0.75f));
}


float GreensFunc(Vec3D x, Vec3D z, float sigma, float R)
{
    float r_x = z.x - x.x;
    float r_y = z.y - x.y;
    float r_z = z.z - x.z;
    float r = sqrt(pow(r_x,2) + pow(r_y,2) + pow(r_z,2));

    float value = (1/(4 * PI)) * ((sinh((r-R) * sqrt(sigma)))/(r * sinh(R * sqrt(sigma))));

    return value;
}

float sourceT(Vec3D cpz, const Surface& surface)
{
    Vec3D sourceCenter = surface.vertices[1000].vert;

    float dx = cpz.x - sourceCenter.x;
    float dy = cpz.y - sourceCenter.y;
    float dz = cpz.z - sourceCenter.z;

    float d2 = dx*dx + dy*dy + dz*dz;

    float sigmaSource = 0.25f;

    return std::exp(-d2 / (2.0f * sigmaSource * sigmaSource));
}

float sampling(Vec3D x, Vec3D z)
{
    float r_x = x.x - z.x;
    float r_y = x.y - z.y;
    float r_z = x.z - z.z;
    float r = sqrt(pow(r_x,2) + pow(r_y,2) + pow(r_z,2));

    float value = 1/r;

    return value;
}

float RecursiveEstimate(const Surface& surface, const std::vector<Vec3D>& M, int Nv, Vec3D x, const mantis::AccelerationStructure& accel)
{
    float r = localFetureSize(M,x);

    float u_hat_sphere = 0.0f;
    float sigma = 5;

    float c = sqrt(sigma) * r / std::sinh(sqrt(sigma) * r);

    float u = (static_cast<float>(rand()) + 1.0f) / (static_cast<float>(RAND_MAX) + 1.0f);

    if (u < c)
    {

        Vec3D y = samplePointsInsideSphere(x, r, 1)[0];
        
        u_hat_sphere = RecursiveEstimate(surface, M, Nv, closestPointOnSurface(accel, y), accel);
    }

    std::vector<Vec3D> z = samplePointsInsideSphere(x, r, Nv);

    float u_hat_ball = 0.0f;

    for (int i = 0; i < Nv; i++)
    {
        if (z[i].x == x.x && z[i].y == x.y && z[i].z == x.z)
        {
            continue;
        }

        u_hat_ball += (GreensFunc(x,  z[i], sigma, r) * sourceT(closestPointOnSurface(accel, z[i]), surface))/sampling(x, z[i]);
    }

    u_hat_ball /= Nv;

    return u_hat_sphere + u_hat_ball;
}

float EstimateSolution(int Np, int Nv, Vec3D x, int Ns, const std::vector<Vec3D>& M, const mantis::AccelerationStructure& accel, const Surface& surface)
{
    float u_hat_sum = 0;

    for (int n = 0; n < Np; n++)
    {
        float u_hat = RecursiveEstimate(surface, M, Nv, x, accel);
        
        u_hat_sum += u_hat;
    }
    
    return u_hat_sum/Np;
}


Model buildModel(const std::string& filename)
{
    std::vector<Vector3> pos, nrm; std::vector<Vector2> uv;
    struct Corner { int v, t, n; };
    std::vector<Corner> corners;

    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Could not open file");
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line); std::string type; ss >> type;
        if (type == "v")  { Vector3 p; ss >> p.x >> p.y >> p.z; pos.push_back(p); }
        else if (type == "vt") { Vector2 t; ss >> t.x >> t.y; uv.push_back(t); }
        else if (type == "vn") { Vector3 n; ss >> n.x >> n.y >> n.z; nrm.push_back(n); }
        else if (type == "f") {
            for (int i = 0; i < 3; i++) {
                std::string tok; ss >> tok;
                int v = 0, t = 0, n = 0;
                sscanf(tok.c_str(), "%d/%d/%d", &v, &t, &n); // "v//vn" gives t=0 -> handle below
                corners.push_back({ v - 1, t - 1, n - 1 });
            }
        }
    }

    Mesh mesh = { 0 };
    mesh.vertexCount = (int)corners.size();
    mesh.triangleCount = mesh.vertexCount / 3;
    mesh.vertices  = (float*)MemAlloc(mesh.vertexCount * 3 * sizeof(float));
    mesh.texcoords = (float*)MemAlloc(mesh.vertexCount * 2 * sizeof(float));
    mesh.normals   = (float*)MemAlloc(mesh.vertexCount * 3 * sizeof(float));

    for (int i = 0; i < mesh.vertexCount; i++) {
        const Corner& c = corners[i];
        mesh.vertices[3*i+0] = pos[c.v].x; mesh.vertices[3*i+1] = pos[c.v].y; mesh.vertices[3*i+2] = pos[c.v].z;
        if (c.t >= 0 && c.t < (int)uv.size()) { mesh.texcoords[2*i] = uv[c.t].x; mesh.texcoords[2*i+1] = 1.0f - uv[c.t].y; }
        if (c.n >= 0 && c.n < (int)nrm.size()) { mesh.normals[3*i] = nrm[c.n].x; mesh.normals[3*i+1] = nrm[c.n].y; mesh.normals[3*i+2] = nrm[c.n].z; }
    }

    UploadMesh(&mesh, false);
    return LoadModelFromMesh(mesh);
}

void updateOrbit(Camera3D& cam, OrbitCam& o)
{
    Vector2 d = GetMouseDelta();

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) 
    {
        o.yaw   -= d.x * 0.005f;
        o.pitch += d.y * 0.005f;
        o.pitch = std::clamp(o.pitch, -1.5f, 1.5f);
    }

    o.dist *= 1.0f - GetMouseWheelMove() * 0.1f;
    o.dist = std::clamp(o.dist, 0.2f, 100.0f);

    cam.position = Vector3Add(o.target, {
        o.dist * cosf(o.pitch) * sinf(o.yaw),
        o.dist * sinf(o.pitch),
        o.dist * cosf(o.pitch) * cosf(o.yaw)
    });

    cam.target = o.target;

    if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) 
    {
        Vector3 fwd   = Vector3Normalize(Vector3Subtract(cam.target, cam.position));
        Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, cam.up));
        Vector3 up    = Vector3CrossProduct(right, fwd);
        float s = o.dist * 0.0015f;
        o.target = Vector3Add(o.target, Vector3Scale(right, -d.x * s));
        o.target = Vector3Add(o.target, Vector3Scale(up,     d.y * s));
    }
}

mantis::AccelerationStructure buildAccel(const Surface& S)
{
    std::vector<std::array<float, 3>> P; 
    P.reserve(S.vertices.size());
    std::vector<std::array<uint32_t, 3>> T; 
    T.reserve(S.triangles.size());

    for (const Vertex& v : S.vertices)
    {
        P.push_back({ v.vert.x, v.vert.y, v.vert.z });
    }
    for (const Triangle& t : S.triangles)
    {
        T.push_back({ (uint32_t)t.v0, (uint32_t)t.v1, (uint32_t)t.v2 });
    }

    return mantis::AccelerationStructure(P, T);
}

int main()
{
    const int screenWidth = 1200;
    const int screenHeight = 800;

    int Np = 100;
    int Nv = 8;
    int Ns = 10000;

    Surface surface = loadOBJ("testModels/spot/spot/spot_triangulated.obj");

    InitWindow(screenWidth, screenHeight, "WoS - Medial Axis");

    Camera3D camera = { 0 };

    camera.position = { 3.0f, 2.0f, -3.0f };
    camera.target = { 0.0f, 0.0f, 0.0f };
    camera.up = { 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    Model model = buildModel("testModels/spot/spot/spot_triangulated.obj");

    Texture2D texture = LoadTexture("testModels/spot/spot/spot_texture.png");

    model.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = texture;

    SetTargetFPS(60);

    Vec3D x = surface.vertices[0].vert;

    mantis::AccelerationStructure accel = buildAccel(surface);

    auto [M, center, r] = MedialAxisPointCloud(surface, Ns, accel);

    //float u = EstimateSolution(Np, Nv, x, Ns, M, accel);

    OrbitCam orbit;

    Vec3D bb_min = surface.vertices[0].vert;
    Vec3D bb_max = surface.vertices[0].vert;

    static bool wires = false;
    static bool bb = false;
    static bool mPAC = false;

    std::vector<float> uVals(surface.vertices.size(), 0.0f);
    std::vector<bool>  uDone(surface.vertices.size(), false);
    size_t nextVertex = 0;
    float uMin = 1e30f;
    float uMax = -1e30f;

    bool calculationFinished = false;

    static bool showU = false;

    while (!WindowShouldClose())
    {
        UpdateCamera(&camera, CAMERA_CUSTOM);
        updateOrbit(camera, orbit);
        if (IsKeyPressed(KEY_W))
        {
            wires = !wires; 
        } 
        if (IsKeyPressed(KEY_B))
        {
            bb = !bb; 
        }
        if (IsKeyPressed(KEY_P))
        {
            mPAC = !mPAC; 
        }

        if (IsKeyPressed(KEY_U)) showU = !showU;

        if (showU && !calculationFinished)
        {
            double t0 = GetTime();
            while (nextVertex < surface.vertices.size() && GetTime() - t0 < 0.01)   // ~10 ms of work per frame
            {
                float uv = EstimateSolution(Np, Nv, surface.vertices[nextVertex].vert, Ns, M, accel, surface);
                uVals[nextVertex] = uv;
                uDone[nextVertex] = true;
                uMin = std::min(uMin, uv);
                uMax = std::max(uMax, uv);
                nextVertex++;
            }

            if (nextVertex == surface.vertices.size())
            {
                calculationFinished = true;
            }
        }

        BeginDrawing();

        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

        if (!showU)
        {
            if (wires) 
            {
                DrawModelWires(model, { 0.0f, 0.0f, 0.0f }, 1.0f, LIGHTGRAY);
            }
            else
            {
                DrawModel(model, { 0.0f, 0.0f, 0.0f }, 1.0f, WHITE);
            }
        }
        //Medial point axis cloud
        if (mPAC)
        {
            for (const Vec3D m : M)
            {
                DrawSphere({m.x, m.y, m.z}, 0.01, RED);
            }
        }
        //boundingbox
        if (bb)
        {
            DrawSphere({center.x, center.y, center.z}, r, { 0, 121, 241, 128 });
        }

        if (showU && calculationFinished)
        {
            for (const Triangle& t : surface.triangles)
            {
                float uAvg = (uVals[t.v0] + uVals[t.v1] + uVals[t.v2]) / 3.0f;
                float s = (uMax > uMin) ? (uAvg - uMin) / (uMax - uMin) : 0.5f;
                Color col = ColorFromHSV((1.0f - s) * 240.0f, 0.9f, 1.0f);       // blue = low, red = high



                Vec3D a = surface.vertices[t.v0].vert;
                Vec3D b = surface.vertices[t.v1].vert;
                Vec3D cc = surface.vertices[t.v2].vert;


                DrawTriangle3D({a.x, a.y, a.z}, {b.x, b.y, b.z}, {cc.x, cc.y, cc.z}, col);
            }
        }
                
        EndMode3D();

        EndDrawing();
    }

    UnloadTexture(texture);
    UnloadModel(model);

    CloseWindow();

    return 0;
}