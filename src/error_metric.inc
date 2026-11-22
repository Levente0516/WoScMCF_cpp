#pragma region Error metric (deterministic reference step)

// ---------------------------------------------------------------------------
// Reference for ONE backward Euler step of the flow, computed deterministically on S0:
//
//        ( diag(A_n) + delta * L0 ) u = diag(A_n) * Phi^n
//
// L0   : cotangent stiffness matrix of the ORIGINAL surface S0 (positive semi-definite)
// A_n  : vertex areas of the current shape Phi^n   (so sigma_n * A_0 = A_n / delta, as in computeSigma)
// u    : reference for Phi^{n+1} (before normalizeGeometry)
//
// The Monte Carlo result (uPrev after the last pass) is compared with u.
// Needs: Surface, Triangle, Vector3 (raylib/raymath), <vector>, <cmath>, <algorithm>.
// ---------------------------------------------------------------------------

struct CotanOp
{
    std::vector<int> i, j;        // one entry per (triangle, corner): edge (i, j) opposite to the corner
    std::vector<double> w;        // 0.5 * cot(angle at that corner)
    std::vector<double> diag;     // diagonal of L0
};

// L0 x = sum_e w_e (x_i - x_j) (e_i - e_j)   -- build once from the positions of S0
CotanOp buildCotanOp(const Surface& S, const std::vector<Vector3>& P)
{
    CotanOp L;
    L.diag.assign(S.vertices.size(), 0.0);

    for (const Triangle& t : S.triangles)
    {
        const int idx[3] = { t.v0, t.v1, t.v2 };

        for (int k = 0; k < 3; ++k)
        {
            const int o = idx[k], a = idx[(k + 1) % 3], b = idx[(k + 2) % 3];

            const double ux = (double)P[a].x - P[o].x, uy = (double)P[a].y - P[o].y, uz = (double)P[a].z - P[o].z;
            const double vx = (double)P[b].x - P[o].x, vy = (double)P[b].y - P[o].y, vz = (double)P[b].z - P[o].z;

            const double cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx;
            const double cross = std::sqrt(cx * cx + cy * cy + cz * cz);

            if (cross < 1e-18)
            {
                continue;                                   // degenerate triangle
            }

            const double w = 0.5 * (ux * vx + uy * vy + uz * vz) / cross;

            L.i.push_back(a);
            L.j.push_back(b);
            L.w.push_back(w);

            L.diag[a] += w;
            L.diag[b] += w;
        }
    }

    return L;
}

// y = (diag(M) + delta * L0) x
static void applyOperator(const CotanOp& L, const std::vector<double>& M, double delta, const std::vector<double>& x, std::vector<double>& y)
{
    for (size_t k = 0; k < x.size(); ++k)
    {
        y[k] = M[k] * x[k];
    }

    for (size_t e = 0; e < L.w.size(); ++e)
    {
        const double d = delta * L.w[e] * (x[L.i[e]] - x[L.j[e]]);

        y[L.i[e]] += d;
        y[L.j[e]] -= d;
    }
}

// Jacobi preconditioned conjugate gradient for one coordinate
static int solveCG(const CotanOp& L, const std::vector<double>& M, double delta, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-10, int maxIter = 5000)
{
    const size_t N = b.size();

    std::vector<double> r(N), z(N), p(N), Ap(N), invD(N);

    for (size_t k = 0; k < N; ++k)
    {
        invD[k] = 1.0 / (M[k] + delta * L.diag[k]);
    }

    applyOperator(L, M, delta, x, Ap);

    double bnorm = 0.0;

    for (size_t k = 0; k < N; ++k)
    {
        r[k] = b[k] - Ap[k];
        bnorm += b[k] * b[k];
    }

    bnorm = std::sqrt(bnorm) + 1e-300;

    double rz = 0.0;

    for (size_t k = 0; k < N; ++k)
    {
        z[k] = invD[k] * r[k];
        p[k] = z[k];
        rz += r[k] * z[k];
    }

    int it = 0;

    for (; it < maxIter; ++it)
    {
        double rnorm = 0.0;

        for (size_t k = 0; k < N; ++k)
        {
            rnorm += r[k] * r[k];
        }

        if (std::sqrt(rnorm) < tol * bnorm)
        {
            break;
        }

        applyOperator(L, M, delta, p, Ap);

        double pAp = 0.0;

        for (size_t k = 0; k < N; ++k)
        {
            pAp += p[k] * Ap[k];
        }

        const double alpha = rz / pAp;

        for (size_t k = 0; k < N; ++k)
        {
            x[k] += alpha * p[k];
            r[k] -= alpha * Ap[k];
        }

        double rzNew = 0.0;

        for (size_t k = 0; k < N; ++k)
        {
            z[k] = invD[k] * r[k];
            rzNew += r[k] * z[k];
        }

        const double beta = rzNew / rz;
        rz = rzNew;

        for (size_t k = 0; k < N; ++k)
        {
            p[k] = z[k] + beta * p[k];
        }
    }

    return it;
}

// Reference for Phi^{n+1}: An = vertexAreas(surface, Phi), Phi = current shape
std::vector<Vector3> referenceStep(const CotanOp& L0, const std::vector<float>& An, const std::vector<Vector3>& Phi, double delta)
{
    const size_t N = Phi.size();

    std::vector<double> M(An.begin(), An.end());
    std::vector<Vector3> u(N);

    for (int c = 0; c < 3; ++c)
    {
        std::vector<double> b(N), x(N);

        for (size_t k = 0; k < N; ++k)
        {
            const double phi = (c == 0) ? Phi[k].x : (c == 1) ? Phi[k].y : Phi[k].z;

            b[k] = M[k] * phi;
            x[k] = phi;                                      // initial guess: no movement
        }

        solveCG(L0, M, delta, b, x);

        for (size_t k = 0; k < N; ++k)
        {
            if (c == 0) u[k].x = (float)x[k];
            else if (c == 1) u[k].y = (float)x[k];
            else u[k].z = (float)x[k];
        }
    }

    return u;
}

struct StepError
{
    double relDisp;     // RMS|u - uRef| / RMS|uRef - Phi|   (1 = as wrong as not moving at all, 0 = exact)
    double relBBox;     // RMS|u - uRef| / bounding box diagonal of Phi
    double maxBBox;     // max|u - uRef| / bounding box diagonal of Phi
};

StepError stepError(const std::vector<Vector3>& u, const std::vector<Vector3>& uRef, const std::vector<Vector3>& Phi)
{
    const size_t N = Phi.size();

    double e2 = 0.0, d2 = 0.0, eMax = 0.0;
    Vector3 bmin = Phi[0], bmax = Phi[0];

    for (size_t i = 0; i < N; ++i)
    {
        const double e = Vector3Distance(u[i], uRef[i]);
        const double d = Vector3Distance(uRef[i], Phi[i]);

        e2 += e * e;
        d2 += d * d;
        eMax = std::max(eMax, e);

        bmin = Vector3Min(bmin, Phi[i]);
        bmax = Vector3Max(bmax, Phi[i]);
    }

    const double rmsE = std::sqrt(e2 / N);
    const double rmsD = std::sqrt(d2 / N);
    const double diag = Vector3Distance(bmin, bmax);

    return { rmsE / std::max(rmsD, 1e-30), rmsE / diag, eMax / diag };
}

#pragma endregion
