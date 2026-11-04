#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <memory>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <random>
#include "raylib.h"
#include "raymath.h"
#include "mantis.h"

#pragma region Structs

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

#pragma endregion

#pragma region Random
//generate the seed TODO: seed should be preprocessed
static std::mt19937& rng()
{
    static std::mt19937 gen(123);
    return gen;
}

//uniform number in [0, 1)
static double uni01()
{
    static std::uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng());
}

#pragma endregion

#pragma region Sampling

//uniform random direction on the unit sphere
static Vector3 randomUnitVector()
{
    double z   = 2.0 * uni01() - 1.0;
    double phi = 2.0 * PI * uni01();
    double s   = std::sqrt(std::max(0.0, 1.0 - z * z));

    return { (float)(s * std::cos(phi)), (float)(s * std::sin(phi)), (float)z };
}

//Samples a point on the given sphere
Vector3 samplePointOnSphere(Vector3 center, float r)
{
    return Vector3Add(center, Vector3Scale(randomUnitVector(), r));
}

//Samples multiple points inside a sphere
std::vector<Vector3> samplePointsInsideSphere(Vector3 center, float r, int Ns)
{
    std::vector<Vector3> randomPoints;
    randomPoints.reserve(Ns);

    for (int i = 0; i < Ns; i++)
    {
        float radius = r * (float)std::cbrt(uni01());

        randomPoints.push_back(Vector3Add(center, Vector3Scale(randomUnitVector(), radius)));
    }

    return randomPoints;
}

#pragma endregion 

//Distance-based local feature size from the medial axis point cloud
//TODO: lfsMin should be determined by something rather than just being a burnt in value, because it depends heavily on the size of the model
float localFetureSize(const std::vector<Vector3>& M, Vector3 x, float lfsMin)
{
    float minDistance2 = std::numeric_limits<float>::max();

    for (const Vector3& m : M)
    {
        minDistance2 = std::min(minDistance2, Vector3DistanceSqr(m, x));
    }

    return std::max(0.9f * std::sqrt(minDistance2), lfsMin);
}

#pragma region Triangle methods
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

//Interpolation weights of the closest surface point to q (and optionally the closest point itself)
Interp closestInterp(const Surface& S, const mantis::AccelerationStructure& accel, Vector3 q, Vector3* outPoint = nullptr)
{
    mantis::Result r = accel.calc_closest_point(q.x, q.y, q.z);

    Vector3 p = { r.closest_point[0], r.closest_point[1], r.closest_point[2] };

    if (outPoint) 
    {
        *outPoint = p;
    }

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

//Interp that sits exactly on vertex i
Interp atVertex(int i)
{
    Interp I;

    I.v[0] = I.v[1] = I.v[2] = i;

    I.w[0] = 1.0f;
    I.w[1] = 0.0f;
    I.w[2] = 0.0f;

    return I;
}

float interpScalar(const Interp& I, const std::vector<float>& f)
{
    return I.w[0] * f[I.v[0]] + I.w[1] * f[I.v[1]] + I.w[2] * f[I.v[2]];
}

Vector3 interpVec(const Interp& I, const std::vector<Vector3>& f)
{
    return Vector3Add(Vector3Add(Vector3Scale(f[I.v[0]], I.w[0]),
                                 Vector3Scale(f[I.v[1]], I.w[1])),
                                 Vector3Scale(f[I.v[2]], I.w[2]));
}

#pragma endregion

//Returns the closest point on a surface from a another point
Vector3 closestPointOnSurface(const mantis::AccelerationStructure& accel, Vector3 q)
{
    mantis::Result r = accel.calc_closest_point(q.x, q.y, q.z);

    return { r.closest_point[0], r.closest_point[1], r.closest_point[2] };
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

#pragma region Medial axis

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
std::tuple<std::vector<Vector3>, Vector3, float> MedialAxisPointCloud(const Surface& surface, int Ns, const mantis::AccelerationStructure& accel)
{
    Vector3 bbMin = surface.vertices[0].vert;
    Vector3 bbMax = surface.vertices[0].vert;

    for (const Vertex& v : surface.vertices)
    {
        bbMin = Vector3Min(bbMin, v.vert);
        bbMax = Vector3Max(bbMax, v.vert);
    }

    Vector3 center = Vector3Scale(Vector3Add(bbMax, bbMin), 0.5f);

    Vector3 bb_range = Vector3Subtract(bbMax, bbMin);

    float r = Vector3Length(bb_range) / 2;

    std::vector<Vector3> medialBallPoints = samplePointsInsideSphere(center, r, Ns);

    std::vector<Vector3> medialAxisPoints;
    medialAxisPoints.reserve(Ns * 2);

    std::vector<float> r_new;
    r_new.reserve(Ns * 2);

    for (const Vector3& points : medialBallPoints)
    {
        Vector3 p = closestPointOnSurface(accel, points);

        Vector3 pp = Vector3Subtract(points, p);

        Vector3 normal = Vector3Normalize(pp);

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
    prunedPoints.reserve(medialAxisPoints.size());

    for (size_t i = 0; i < medialAxisPoints.size(); ++i)
    {
        if (!removed[i])
        {
            prunedPoints.push_back(medialAxisPoints[i]);
        }
    }

    return {prunedPoints, center, r};
}

#pragma endregion

#pragma region Screened Poisson ball helpers

//P(no collision inside the ball) = kappa / sinh(kappa), kappa = R * sqrt(sigmaBar)
static double continueProb(double k)
{
    if (k < 1e-4)
    {
        return 1.0 - k * k / 6.0;
    }

    return 2.0 * k * std::exp(-k) / (1.0 - std::exp(-2.0 * k));
}

//Radial density of G_sigma in the ball (t = r / R):  t * sinh(kappa (1 - t))
static double radialDensity(double t, double k)
{
    return t * (std::exp(-k * t) - std::exp(-k * (2.0 - t)));
}

//Maximum of the radial density: tanh(k (1 - t)) = k t
static double radialDensityMax(double k)
{
    double lo = 0.0, hi = 1.0;

    for (int i = 0; i < 40; ++i)
    {
        double mid = 0.5 * (lo + hi);

        if (std::tanh(k * (1.0 - mid)) - k * mid > 0.0)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    return radialDensity(0.5 * (lo + hi), k) * 1.0001;
}

//Samples y ~ G_sigmaBar(x, .) / int G in the ball of radius R around x (rejection sampling on the radius)
//k = R * sqrt(sigmaBar), Mx = radialDensityMax(k)
static Vector3 sampleScreenedBall(Vector3 x, float R, double k, double Mx)
{
    double t;

    do 
    { 
        t = uni01(); 
    } 
    while (uni01() * Mx > radialDensity(t, k));

    return Vector3Add(x, Vector3Scale(randomUnitVector(), (float)(t * R)));
}

#pragma endregion

#pragma region Geometry helpers

//Area belonging to every vertex (1/3 of each adjacent triangle) for the positions P
std::vector<float> vertexAreas(const Surface& S, const std::vector<Vector3>& P)
{
    std::vector<float> A(S.vertices.size(), 0.0f);

    for (const Triangle& t : S.triangles)
    {
        float a = 0.5f * Vector3Length(Vector3CrossProduct(Vector3Subtract(P[t.v1], P[t.v0]), Vector3Subtract(P[t.v2], P[t.v0]))) / 3.0f;

        A[t.v0] += a;
        A[t.v1] += a;
        A[t.v2] += a;
    }

    return A;
}

//sigma_n = (1/delta) * dA_n / dA_0   per vertex  (A_0: areas of the surface the walks live on, A_n: areas of Phi^n)
std::vector<float> computeSigma(const Surface& S, const std::vector<Vector3>& Phi, const std::vector<float>& A0, float delta)
{
    std::vector<float> An = vertexAreas(S, Phi);
    std::vector<float> sigma(An.size());

    for (size_t i = 0; i < sigma.size(); ++i)
    {
        sigma[i] = An[i] / (std::max(A0[i], 1e-12f) * delta);
    }

    return sigma;
}

//Isoperimetric quotient pi^(1/3) (6V)^(2/3) / A : 1 for a perfect sphere, smaller for everything else
double sphericity(const Surface& S, const std::vector<Vector3>& P, double* areaOut = nullptr)
{
    double V = 0.0;
    double A = 0.0;

    for (const Triangle& t : S.triangles)
    {
        Vector3 a = P[t.v0], b = P[t.v1], c = P[t.v2];

        V += Vector3DotProduct(a, Vector3CrossProduct(b, c)) / 6.0;
        A += 0.5 * Vector3Length(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
    }

    if (areaOut) 
    {
        *areaOut = A;
    }

    return std::cbrt((double)PI) * std::pow(6.0 * std::fabs(V), 2.0 / 3.0) / A;
}

//Moves the centroid to the origin and rescales so the total area is targetArea
void normalizeGeometry(const Surface& S, std::vector<Vector3>& P, float targetArea)
{
    std::vector<float> A = vertexAreas(S, P);

    double total = 0.0;
    Vector3 c = { 0.0f, 0.0f, 0.0f };

    for (size_t i = 0; i < P.size(); ++i)
    {
        total += A[i];
        c = Vector3Add(c, Vector3Scale(P[i], A[i]));
    }

    c = Vector3Scale(c, (float)(1.0 / total));

    float s = (float)std::sqrt(targetArea / total);

    for (Vector3& p : P)
    {
        p = Vector3Scale(Vector3Subtract(p, c), s);
    }
}

#pragma endregion

#pragma region cMCF estimators

//Low sample estimate: ONE walk (with null collisions) that starts at vertex 'start', returns Phi^n(Y) at the end point Y
Vector3 WalkEstimate(int start, const Surface& surface, const mantis::AccelerationStructure& accel, const std::vector<float>& lfs, const std::vector<Vector3>& Phi, const std::vector<float>& sigma, float sigmaBar)
{
    const int maxWalkSteps = 10000;

    Vector3 x = surface.vertices[start].vert;
    Interp Ix = atVertex(start);

    for (int step = 0; step < maxWalkSteps; ++step)
    {
        float R = interpScalar(Ix, lfs); //feature size at x
        double k = std::sqrt((double)sigmaBar) * R;

        if (uni01() < continueProb(k))
        {
            //no collision in the ball: uniform point on the sphere, projected back to the surface
            Vector3 y = samplePointOnSphere(x, R);
            Vector3 p;

            Ix = closestInterp(surface, accel, y, &p);
            x  = p;
        }

        else
        {
            //collision in the ball: y ~ G, evaluated at its projection onto the surface
            Vector3 y = sampleScreenedBall(x, R, k, radialDensityMax(k));
            Vector3 p;

            Interp Iy = closestInterp(surface, accel, y, &p);

            if (uni01() < interpScalar(Iy, sigma) / sigmaBar)
            {
                return interpVec(Iy, Phi);                 //real collision: Phi^n(Y)
            }

            x  = p;                                        //null collision: continue from y
            Ix = Iy;
        }
    }

    return interpVec(Ix, Phi);                             //safety net, should (almost) never happen
}

Vector3 EstimateSolution(int Np, int start, const Surface& surface, const mantis::AccelerationStructure& accel, const std::vector<float>& lfs, const std::vector<Vector3>& Phi, const std::vector<float>& sigma, float sigmaBar)
{
    Vector3 sum = { 0.0f, 0.0f, 0.0f };

    for (int n = 0; n < Np; n++)
    {
        sum = Vector3Add(sum, WalkEstimate(start, surface, accel, lfs, Phi, sigma, sigmaBar));
    }

    return Vector3Scale(sum, 1.0f / (float)Np);
}

//Filtered estimate at vertex i: ONE step, the rest comes from the interpolated previous estimate uPrev
Vector3 FilteredEstimate(int i, int Nv, const Surface& surface, const mantis::AccelerationStructure& accel, const std::vector<float>& lfs, const std::vector<Vector3>& Phi, const std::vector<float>& sigma, float sigmaBar, const std::vector<Vector3>& uPrev)
{
    Vector3 x = surface.vertices[i].vert;
    float R = lfs[i];

    double k = std::sqrt((double)sigmaBar) * R;
    double Mx = radialDensityMax(k);
    double c = continueProb(k);

    //sphere term:  c * u_prev(z),  z uniform on the sphere
    Interp Iz = closestInterp(surface, accel, samplePointOnSphere(x, R));

    Vector3 sphereTerm = Vector3Scale(interpVec(Iz, uPrev), (float)c);

    //ball term:  (1 - c) * mean over y ~ G of  [ p Phi^n(y) + (1 - p) u_prev(y) ],   p = sigma_n(y) / sigmaBar
    Vector3 ball = { 0.0f, 0.0f, 0.0f };

    for (int j = 0; j < Nv; ++j)
    {
        Interp Iy = closestInterp(surface, accel, sampleScreenedBall(x, R, k, Mx));

        float p = std::clamp(interpScalar(Iy, sigma) / sigmaBar, 0.0f, 1.0f);

        Vector3 val = Vector3Add(Vector3Scale(interpVec(Iy, Phi),   p), Vector3Scale(interpVec(Iy, uPrev), 1.0f - p));

        ball = Vector3Add(ball, val);
    }

    ball = Vector3Scale(ball, (float)((1.0 - c) / Nv));

    return Vector3Add(sphereTerm, ball);
}

Vector3 FilteredSolution(int Nw, int Nv, int i, const Surface& surface, const mantis::AccelerationStructure& accel, const std::vector<float>& lfs, const std::vector<Vector3>& Phi, const std::vector<float>& sigma, float sigmaBar, const std::vector<Vector3>& uPrev)
{
    Vector3 sum = { 0.0f, 0.0f, 0.0f };

    for (int n = 0; n < Nw; n++)
    {
        sum = Vector3Add(sum, FilteredEstimate(i, Nv, surface, accel, lfs, Phi, sigma, sigmaBar, uPrev));
    }

    return Vector3Scale(sum, 1.0f / (float)Nw);
}

#pragma endregion

Model buildModel(const std::string& filename)
{
    std::vector<Vector3> pos, nrm; std::vector<Vector2> uv;
    struct Corner { int v, t, n; };
    std::vector<Corner> corners;

    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Could not open file");
    std::string line;
    while (std::getline(file, line)) 
    {
        std::stringstream ss(line); std::string type; ss >> type;
        if (type == "v")  
        { 
            Vector3 p; 
            ss >> p.x >> p.y >> p.z; 
            pos.push_back(p); 
        }
        else if (type == "vt") 
        { 
            Vector2 t; 
            ss >> t.x >> t.y; 
            uv.push_back(t); 
        }
        else if (type == "vn") 
        { 
            Vector3 n; 
            ss >> n.x >> n.y >> n.z; 
            nrm.push_back(n); 
        }
        else if (type == "f") 
        {
            for (int i = 0; i < 3; i++) 
            {
                std::string tok; ss >> tok;
                int v = 0, t = 0, n = 0;
                sscanf(tok.c_str(), "%d/%d/%d", &v, &t, &n);
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

    for (int i = 0; i < mesh.vertexCount; i++) 
    {
        const Corner& c = corners[i];
        mesh.vertices[3*i+0] = pos[c.v].x; mesh.vertices[3*i+1] = pos[c.v].y; mesh.vertices[3*i+2] = pos[c.v].z;
        if (c.t >= 0 && c.t < (int)uv.size()) 
        { 
            mesh.texcoords[2*i] = uv[c.t].x;
            mesh.texcoords[2*i+1] = 1.0f - uv[c.t].y; 
        }

        if (c.n >= 0 && c.n < (int)nrm.size()) 
        { 
            mesh.normals[3*i] = nrm[c.n].x; 
            mesh.normals[3*i+1] = nrm[c.n].y; 
            mesh.normals[3*i+2] = nrm[c.n].z; 
        }
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
    o.dist = std::clamp(o.dist, 0.2f, 300.0f);

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

std::unique_ptr<mantis::AccelerationStructure> buildAccel(const Surface& S)
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

    return std::make_unique<mantis::AccelerationStructure>(P, T);
}

int main()
{
    const int screenWidth = 1200;
    const int screenHeight = 800;

    const int Ns = 20000; //medial axis points
    const float lfsMin = 10.0f; //localFetureSize lower bound TODO: should be based on the model

    //flow
    const float delta = 0.01f; // time step of one iteration: kappa = R * sqrt(sigmaBar) = R * sqrt(max(dA_n/dA_0) / delta)
    const int numIterations = 40; // how many times the geometry is moved

    //estimator
    const int Np_low = 20; // pass 0: low sample estimate
    const int Nw = 30; // filter passes per vertex
    const int K = 3; // number of filter passes
    const int Nv = 5; // ball samples per filtered estimate

    Surface surface = loadOBJ("testModels/armadillo/armadillo.obj");
    //Surface surface = loadOBJ("testModels/spot/spot/spot_triangulated.obj");

    InitWindow(screenWidth, screenHeight, "WoScMFC");

    Camera3D camera = { 0 };

    camera.position = { 3.0f, 2.0f, -3.0f };
    camera.target = { 0.0f, 0.0f, 0.0f };
    camera.up = { 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    Model model = buildModel("testModels/armadillo/armadillo.obj");

    //Texture2D texture = LoadTexture("testModels/spot/spot/spot_texture.png");

    //model.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = texture;

    SetTargetFPS(60);

    OrbitCam orbit;

    static bool wires = false;
    static bool bb = false;
    static bool mPAC = false;
    static bool flowOn = false;

    const size_t N = surface.vertices.size();

    std::vector<Vector3> P_orig(N);
    
    for (size_t i = 0; i < N; ++i)
    {
        P_orig[i] = surface.vertices[i].vert;
    }

    const std::vector<float> A_orig = vertexAreas(surface, P_orig);

    const float targetArea = std::accumulate(A_orig.begin(), A_orig.end(), 0.0f);

    std::unique_ptr<mantis::AccelerationStructure> accel;
    std::vector<Vector3> M;
    Vector3 bsCenter = { 0.0f, 0.0f, 0.0f };
    float bsRadius = 1.0f;
    std::vector<float> lfs(N);             
    float meanLfs = 0.0f;
    std::vector<float> A0;                  
    std::vector<Vector3> Phi;               

    //(re)builds everything that depends on the geometry: surface, acceleration structure, medial axis, feature size, areas
    auto rebuildBase = [&](const std::vector<Vector3>& positions)
    {
        for (size_t i = 0; i < N; ++i)
        {
            surface.vertices[i].vert = positions[i];
        }

        accel = buildAccel(surface);

        auto res = MedialAxisPointCloud(surface, Ns, *accel);

        M = std::get<0>(res);
        bsCenter = std::get<1>(res);
        bsRadius = std::get<2>(res);

        for (size_t i = 0; i < N; ++i)
        {
            lfs[i] = localFetureSize(M, positions[i], lfsMin);
        }

        meanLfs = std::accumulate(lfs.begin(), lfs.end(), 0.0f) / (float)N;

        A0  = vertexAreas(surface, positions);
        Phi = positions;
    };

    rebuildBase(P_orig);

    std::vector<Vector3> uPrev(N); 
    std::vector<Vector3> uNew(N);
    std::vector<Vector3> pendingPositions;
    std::vector<float> sigma;

    float sigmaBar = 0.0f;
    bool sigmaReady = false;
    bool needRebuild = false;
    int iteration = 0;
    int pass = 0;
    size_t nextVertex  = 0;

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
        if (IsKeyPressed(KEY_SPACE))
        {
            flowOn = !flowOn;         
        }
        if (IsKeyPressed(KEY_R))
        {
            flowOn = false;
            rebuildBase(P_orig);
            iteration = 0;
            pass = 0;
            nextVertex = 0;
            sigmaReady = false;
            needRebuild = false;
        }

        if (flowOn && iteration < numIterations)
        {
            //start of an iteration: sigma_n and sigmaBar
            if (!sigmaReady)
            {
                sigma = computeSigma(surface, Phi, A0, delta);
                sigmaBar = *std::max_element(sigma.begin(), sigma.end());
                sigmaReady = true;
            }

            double t0 = GetTime();

            while (nextVertex < N && GetTime() - t0 < 0.05)
            {
                if (pass == 0)
                {
                    uNew[nextVertex] = EstimateSolution(Np_low, (int)nextVertex, surface, *accel, lfs, Phi, sigma, sigmaBar);
                }
                else
                {
                    uNew[nextVertex] = FilteredSolution(Nw, Nv, (int)nextVertex, surface, *accel, lfs, Phi, sigma, sigmaBar, uPrev);
                }

                nextVertex++;
            }

            if (nextVertex == N)        
            {
                uPrev.swap(uNew);

                nextVertex = 0;

                pass++;

                if (pass > K)
                {
                    //iteration finished: uPrev = Phi^{n+1}, move the geometry
                    std::vector<Vector3> newPos = uPrev;

                    normalizeGeometry(surface, newPos, targetArea);                    

                    Phi = newPos;

                    iteration++;
                    pass = 0;
                    sigmaReady = false;
                }
            }
        }

        BeginDrawing();

        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

        bool showPhi = flowOn || iteration > 0;

        if (!showPhi)
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
        else
        {
            //the current shape Phi, simple headlight shading
            Vector3 toCam = Vector3Normalize(Vector3Subtract(camera.position, camera.target));

            for (const Triangle& t : surface.triangles)
            {
                Vector3 a = Phi[t.v0];
                Vector3 b = Phi[t.v1];
                Vector3 c = Phi[t.v2];

                if (wires)
                {
                    DrawLine3D(a, b, DARKGRAY);
                    DrawLine3D(b, c, DARKGRAY);
                    DrawLine3D(c, a, DARKGRAY);
                }
                else
                {
                    Vector3 n = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));

                    float shade = 0.3f + 0.7f * std::max(0.0f, Vector3DotProduct(n, toCam));

                    Color col = { (unsigned char)(120.0f * shade), (unsigned char)(170.0f * shade), (unsigned char)(230.0f * shade), 255 };

                    DrawTriangle3D(a, b, c, col);
                }
            }
        }

        //Medial point axis cloud (of the surface the walks currently live on)
        if (mPAC)
        {
            for (const Vector3& m : M)
            {
                DrawSphere({m.x, m.y, m.z}, 0.01, RED);
            }
        }

        //boundingbox
        if (bb)
        {
            DrawSphere({bsCenter.x, bsCenter.y, bsCenter.z}, bsRadius, { 0, 121, 241, 128 });
        }
                
        EndMode3D();

        if (showPhi)
        {
            double area = 0.0;
            double sph = sphericity(surface, Phi, &area);

            DrawText(TextFormat("iteration %d / %d   pass %d / %d   (%d / %d vertices)", iteration, numIterations, pass, K, (int)nextVertex, (int)N), 20, 20, 20, DARKGRAY);

            if (sigmaBar > 0.0f)
            {
                double kMean = std::sqrt((double)sigmaBar) * meanLfs;

                DrawText(TextFormat("sigmaBar = %.1f   kappa(mean R) = %.2f", sigmaBar, kMean), 20, 46, 20, DARKGRAY);
            }

            DrawText(TextFormat("area = %.3f   sphericity = %.4f  (sphere = 1)", area, sph), 20, 72, 20, DARKGRAY);
        }

        DrawText("Space: start/pause   R: reset   W: wires   P: medial axis   B: bounding sphere", 20, screenHeight - 30, 18, GRAY);

        EndDrawing();
    }

    //UnloadTexture(texture);
    UnloadModel(model);

    CloseWindow();

    return 0;
}