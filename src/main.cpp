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

struct Vertex{
    Vector3 vert;
};

struct Triangle{
    int v0, v1, v2;
};

struct Surface{
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
};

struct OrbitCam {
    float yaw = 0.8f; 
    float pitch = 0.5f; 
    float dist = 4.5f;
    Vector3 target = { 0, 0, 0 };
};

struct Interp {
    int v[3];     
    float w[3];    
};

//generate the seed TODO: seed should be preprocessed
static std::mt19937& rng()
{
    static std::mt19937 gen(123);
    return gen;
}

//generates a float between 2 numbers uniformly
static float uni(float low, float high)
{
    static std::uniform_real_distribution<float> d(low, high);
    return d(rng());
}

//Returns the distance between a point on the surface and its closest medial axis point
float localFetureSize(const std::vector<Vector3>& M, Vector3 x)
{
    float minDistance = std::numeric_limits<float>::max();

    for (const Vector3& m : M)
    {
        float d2 = Vector3Distance(m, x);
        minDistance = std::min(minDistance, d2);
    }

    //TODO: 0.25f should be lambda and be determined by somthing reather than just it beign a burnt in value, because it depends heavly on the size of the modell
    return std::max(0.9f * minDistance, 0.25f);
}


float GreensFunc(Vector3 x, Vector3 z, float sigma, float R)
{
    float r_x = z.x - x.x;
    float r_y = z.y - x.y;
    float r_z = z.z - x.z;
    float r = sqrt(pow(r_x,2) + pow(r_y,2) + pow(r_z,2));

    float value = (1/(4 * PI)) * ((sinh((R-r) * sqrt(sigma)))/(r * sinh(R * sqrt(sigma))));

    return value;
}

float sourceT(Vector3 cpz, const Surface& surface)
{
    Vector3 sourceCenter = surface.vertices[1000].vert;
    
    float dx = cpz.x - sourceCenter.x;
    float dy = cpz.y - sourceCenter.y;
    float dz = cpz.z - sourceCenter.z;
    
    float d2 = dx*dx + dy*dy + dz*dz;
    
    float sigmaSource = 0.25f;

    return std::exp(-d2 / (2.0f * sigmaSource * sigmaSource));
}

float sampling(Vector3 x, Vector3 z)
{
    float r_x = x.x - z.x;
    float r_y = x.y - z.y;
    float r_z = x.z - z.z;
    float r = sqrt(pow(r_x,2) + pow(r_y,2) + pow(r_z,2));

    float value = 1/r;

    return value;
}

//Samples a point on the gives sphere
Vector3 samplePointOnSphere(Vector3 center, float r)
{
    float u = uni(0.0f, 1.0f);
    float v = uni(0.0f, 1.0f);
    float phi = 2 * PI * uni(0.0f, 1.0f);

    float theta = acos(2*v - 1);

    float x = center.x + r * sin(theta) * cos(phi);
    float y = center.y + r * sin(theta) * sin(phi);
    float z = center.z + r * cos(theta);

    return {x, y, z};
}   

//Samples multiple points inside a sphere
std::vector<Vector3> samplePointsInsideSphere(Vector3 center, float r, int Ns)
{
    std::vector<Vector3> randomPoints;
    randomPoints.reserve(Ns);

    for (int i = 0; i < Ns; i++)
    {
        Vector3 cords;

        float u = uni(0.0f, 1.0f);
        float v = uni(0.0f, 1.0f);
        float phi = 2 * PI * uni(0.0f, 1.0f);

        float theta = acos(2*v - 1);

        float radius = r * pow(u, (1.0f/3.0f));

        cords.x = center.x + radius * sin(theta) * cos(phi);
        cords.y = center.y + radius * sin(theta) * sin(phi);
        cords.z = center.z + radius * cos(theta);

        randomPoints.push_back(cords);
    }

    return randomPoints;
}

//Turns Cartesian cordinates to Barycentric cords. while filtering out bad triangles
Vector3 barycentric(const Surface& S, int tri, Vector3 p)
{
    const Triangle& t = S.triangles[tri];

    Vector3 a = S.vertices[t.v0].vert;
    Vector3 b = S.vertices[t.v1].vert;
    Vector3 c = S.vertices[t.v2].vert;

    Vector3 n = Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a));

    if (Vector3DotProduct(n, n) < 1e-20f)
    {
        return { 1.0f, 0.0f, 0.0f };
    }

    Vector3 bc = Vector3Barycenter(p, a, b, c);

    bc.x = std::max(bc.x, 0.0f);
    bc.y = std::max(bc.y, 0.0f);
    bc.z = std::max(bc.z, 0.0f);

    float s = bc.x + bc.y + bc.z;

    return { bc.x / s, bc.y / s, bc.z / s };
}

//Interpolates a tirangles value based on its neighbours
Interp closestInterp(const Surface& S, const mantis::AccelerationStructure& accel, Vector3 q)
{
    mantis::Result r = accel.calc_closest_point(q.x, q.y, q.z);

    Vector3 p = { r.closest_point[0], r.closest_point[1], r.closest_point[2] };

    Interp I;

    if (r.type == mantis::Vertex)          
    {
        I.v[0] = r.primitive_index; 
        I.v[1] = r.primitive_index; 
        I.v[2] = r.primitive_index;

        I.w[0] = 1.0f;              
        I.w[1] = 0.0f;              
        I.w[2] = 0.0f;
    }
    else if (r.type == mantis::Edge)
    {
        auto [a, b] = accel.get_edge(r.primitive_index);

        Vector3 pa = S.vertices[a].vert, pb = S.vertices[b].vert;

        float ex = pb.x - pa.x;
        float ey = pb.y - pa.y; 
        float ez = pb.z - pa.z;

        float len2 = ex*ex + ey*ey + ez*ez;

        float t = len2 > 1e-20f
                ? ((p.x - pa.x)*ex + (p.y - pa.y)*ey + (p.z - pa.z)*ez) / len2
                : 0.0f;

        t = std::clamp(t, 0.0f, 1.0f);

        I.v[0] = (int)a; 
        I.v[1] = (int)b; 
        I.v[2] = (int)a;

        I.w[0] = 1.0f - t; 
        I.w[1] = t;    
        I.w[2] = 0.0f;
    }
    else                                    
    {
        const Triangle& tr = S.triangles[r.primitive_index];

        Vector3 bc = barycentric(S, r.primitive_index, p);

        I.v[0] = tr.v0; 
        I.v[1] = tr.v1; 
        I.v[2] = tr.v2;

        I.w[0] = bc.x;  
        I.w[1] = bc.y;  
        I.w[2] = bc.z;
    }

    return I;
}

//Returns the closest point on a surface from a another point
Vector3 closestPointOnSurface(const mantis::AccelerationStructure& accel, Vector3 q)
{
    mantis::Result r = accel.calc_closest_point(q.x, q.y, q.z);

    return { r.closest_point[0], r.closest_point[1], r.closest_point[2] };
}

//Filter usage
float FilteredEstimate(const Surface& surface, const std::vector<Vector3>& M, int Nv, Vector3 x, const mantis::AccelerationStructure& accel, const std::vector<float>& uPrev)
{
    float r = localFetureSize(M, x);

    float sigma = 5;

    float c = sqrt(sigma) * r / std::sinh(sqrt(sigma) * r);

    float u_hat_sphere = 0.0f;

    float u = uni(0.0f, 1.0f);

    if (u < c)
    {
        Vector3 y = samplePointOnSphere(x, r);

        Interp I = closestInterp(surface, accel, y);

        u_hat_sphere = I.w[0] * uPrev[I.v[0]] + I.w[1] * uPrev[I.v[1]] + I.w[2] * uPrev[I.v[2]];
    }

    std::vector<Vector3> z = samplePointsInsideSphere(x, r, Nv);

    float u_hat_ball = 0.0f;

    for (int i = 0; i < Nv; i++)
    {
        if (z[i].x == x.x && z[i].y == x.y && z[i].z == x.z)
        {
            continue;
        }

        u_hat_ball += (GreensFunc(x, z[i], sigma, r) * sourceT(closestPointOnSurface(accel, z[i]), surface)) / sampling(x, z[i]);
    }

    u_hat_ball /= Nv;

    return u_hat_sphere + u_hat_ball;
}

//Filter usage
float FilteredSolution(int Nw, int Nv, Vector3 x, const std::vector<Vector3>& M, const mantis::AccelerationStructure& accel, const Surface& surface, const std::vector<float>& uPrev)
{
    float sum = 0.0f;

    for (int n = 0; n < Nw; n++)
    {
        sum += FilteredEstimate(surface, M, Nv, x, accel, uPrev);
    }

    return sum / Nw;
}

//Loads in the OBJ file provided
Surface loadOBJ(const std::string& filename)
{
    Surface surface;

    std::ifstream file(filename);

    if (!file) 
    { 
        throw std::runtime_error("Could not open file"); 
    }

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

//Implements the algorithm of Ma et. al.
std::tuple<Vector3, float> MaEtAl(Vector3 n, float R, Vector3 p, const mantis::AccelerationStructure& accel)
{
    const float tol = 1e-4f * R;
    const int maxIter = 1000;                      

    float r = R;

    Vector3 c = Vector3Add(p, (Vector3Scale(n, r))); 

    Vector3 q = closestPointOnSurface(accel, c);

    Vector3 q_prev = { 1e30f, 1e30f, 1e30f };         

    int it = 0;

    while (Vector3Distance(q, q_prev) > tol && it++ < maxIter)       
    {
        if (Vector3Distance(q, p) < tol) 
        {
            break;                      
        }

        Vector3 d = Vector3Subtract(q, p);

        float denom = 2.0f * Vector3DotProduct(n,d);

        if (denom <= 1e-12f) 
        {
            break;
        }

        r = Vector3DotProduct(d,d) / denom; 
        c = Vector3Add(p, (Vector3Scale(n, r)));            
        q_prev = q;
        q = closestPointOnSurface(accel, c);
    }

    return { c, r };
}



//Creates the Medial axis Point Cloud that will be used for how big of a sphere we should choose along the steps of a walk
std::tuple<std::vector<Vector3>, Vector3, float> MedialAxisPointCloud(Model model, int Ns, const mantis::AccelerationStructure& accel)
{
    BoundingBox bb = GetModelBoundingBox(model);

    Vector3 center = Vector3Scale(Vector3Add(bb.max, bb.min), 0.5f);

    Vector3 bb_range = Vector3Subtract(bb.max, bb.min);

    float r = Vector3Length(bb_range) / 2;

    std::vector<Vector3> medialBallPoints = samplePointsInsideSphere(center, r, Ns);

    std::vector<Vector3> medialAxisPoints;
    medialAxisPoints.reserve(Ns);

    std::vector<float> r_new;
    r_new.reserve(Ns);

    for (const Vector3 points : medialBallPoints)
    {

        Vector3 p = closestPointOnSurface(accel, points);

        Vector3 pp = Vector3Subtract(points, p);

        float magnitude = Vector3Length(pp);

        Vector3 normal = Vector3Scale(pp, magnitude);

        Vector3 normalInv = Vector3Scale(normal, -1);

        auto [medialAxisPointNormal, normalR] = MaEtAl(normal, r*2, p, accel);
        auto [medialAxisPointNormalInv, normalInvR] = MaEtAl(normalInv, r*2, p, accel);

        float min_r = std::min(normalR, normalInvR);

        if (normalR > normalInvR)
        {
            medialAxisPointNormal = Vector3Add(p, Vector3Scale(normal, normalInvR));
        }
        else if (normalInvR > normalR)
        {
            medialAxisPointNormalInv = Vector3Add(p, Vector3Scale(normalInv, normalR));
        }

        r_new.push_back(min_r);
        r_new.push_back(min_r);
        medialAxisPoints.push_back(medialAxisPointNormal);
        medialAxisPoints.push_back(medialAxisPointNormalInv);
    }

    //TODO: make s preprocessed , or given as a function arugment
    float s = 1.15f;

    std::vector<bool> removed(medialAxisPoints.size(), false);

    for (size_t i = 0; i < medialAxisPoints.size(); ++i)
    {
        for (size_t j = i + 1; j < medialAxisPoints.size(); ++j)
        {
            float distance = Vector3Distance(medialAxisPoints[i], medialAxisPoints[j]); 

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

    std::vector<Vector3> prunedPoints;
    prunedPoints.reserve(Ns * 2 - removed.size());

    for (size_t i = 0; i < medialAxisPoints.size(); ++i)
    {
        if (!removed[i])
        {
            prunedPoints.push_back(medialAxisPoints[i]);
        }
    }

    return {prunedPoints, center, r};
}

float RecursiveEstimate(const Surface& surface, const std::vector<Vector3>& M, int Nv, Vector3 x, const mantis::AccelerationStructure& accel)
{
    float r = localFetureSize(M,x);

    float u_hat_sphere = 0.0f;
    float sigma = 5;

    float c = sqrt(sigma) * r / std::sinh(sqrt(sigma) * r);

    float u = uni(0.0f, 1.0f);

    if (u < c)
    {
        Vector3 y = samplePointOnSphere(x, r);
        
        u_hat_sphere = RecursiveEstimate(surface, M, Nv, closestPointOnSurface(accel, y), accel);
    }

    std::vector<Vector3> z = samplePointsInsideSphere(x, r, Nv);

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

float EstimateSolution(int Np, int Nv, Vector3 x, int Ns, const std::vector<Vector3>& M, const mantis::AccelerationStructure& accel, const Surface& surface)
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

    const int Np_low = 1;     // kezdőbecslés: régi módszer, kevés minta
    const int Nw     = 32;    // szűrőlépésenként ennyi minta csúcsonként
    const int K      = 5;     // szűrőlépések száma
    int Nv = 16;
    int Ns = 10000;

    //unsigned int seed = 123;

    Surface surface = loadOBJ("testModels/spot/spot/spot_triangulated.obj");

    InitWindow(screenWidth, screenHeight, "WoScMFC");

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

    mantis::AccelerationStructure accel = buildAccel(surface);

    auto [M, center, r] = MedialAxisPointCloud(model, Ns, accel);

    OrbitCam orbit;

    Vector3 bb_min = surface.vertices[0].vert;
    Vector3 bb_max = surface.vertices[0].vert;

    static bool wires = false;
    static bool bb = false;
    static bool mPAC = false;
    static bool showU = false;

    const size_t N = surface.vertices.size();

    std::vector<float> uPrev(N, 0.0f);
    std::vector<float> uNew(N, 0.0f);
    std::vector<float> uVals(N, 0.0f);     

    size_t nextVertex = 0;
    int pass = 0;

    float uMin = 1e30f;
    float uMax = -1e30f;

    bool calculationFinished = false;


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

        if (IsKeyPressed(KEY_U))
        {
            showU = !showU;
        }

        if (showU && !calculationFinished)
        {
            double t0 = GetTime();

            while (nextVertex < N && GetTime() - t0 < 0.01)
            {
                Vector3 xi = surface.vertices[nextVertex].vert;

                if (pass == 0)
                {
                    uNew[nextVertex] = EstimateSolution(Np_low, Nv, xi, Ns, M, accel, surface);
                }
                else
                {
                    uNew[nextVertex] = FilteredSolution(Nw, Nv, xi, M, accel, surface, uPrev);
                }

                nextVertex++;
            }

            if (nextVertex == N)            // egy teljes menet kész
            {
                uPrev.swap(uNew);

                nextVertex = 0;

                pass++;

                if (pass > K)
                {
                    uVals = uPrev;
                    auto [mn, mx] = std::minmax_element(uVals.begin(), uVals.end());
                    uMin = *mn;
                    uMax = *mx;
                    calculationFinished = true;
                }
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
            for (const Vector3 m : M)
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
                Color col = ColorFromHSV((1.0f - s) * 240.0f, 0.9f, 1.0f);       

                Vector3 a = surface.vertices[t.v0].vert;
                Vector3 b = surface.vertices[t.v1].vert;
                Vector3 cc = surface.vertices[t.v2].vert;

                DrawTriangle3D({a.x, a.y, a.z}, {b.x, b.y, b.z}, {cc.x, cc.y, cc.z}, col);
            }

            /*
            Vector3 sc = surface.vertices[1000].vert;
            DrawSphere({sc.x, sc.y, sc.z}, 0.03f, BLACK);  
            */
        }
                
        EndMode3D();

        if (showU && !calculationFinished)
        {
            DrawText(TextFormat("pass %d / %d  (%d / %d vertices)", pass, K, (int)nextVertex, (int)N), 20, 20, 20, DARKGRAY);
        }

        EndDrawing();
    }

    UnloadTexture(texture);
    UnloadModel(model);

    CloseWindow();

    return 0;
}