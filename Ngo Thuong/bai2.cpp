#include <bits/stdc++.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

MatrixXd P;
MatrixXd C;
MatrixXd B;
MatrixXd S;
MatrixXd D;
MatrixXd K;

VectorXd d;
VectorXd x0;
VectorXd x1;
VectorXd x2;

double alpha;
double sigma;
double v;
double r;
double L;
double alpha_n;
double m;
double M;
double M1;
double M2;
double gamma_ = 1.5;
double theta = 0.99;

// F(x) = K*x
MatrixXd Fx(MatrixXd x)
{
    return K * x;
}

// C(x) = 1/2 * ( (x-d)T * P * (x-d) - r*r)
double Cx(MatrixXd x)
{
    MatrixXd x_minus_d = x - d;
    MatrixXd Cx = 0.5 * (((x_minus_d.transpose()) * P) * (x - d));
    return Cx.sum() - 0.5 * r * r;
}

// Tích vô hướng của x và y
double scalar(MatrixXd x, MatrixXd y)
{
    MatrixXd x_prod_y = x.array() * y.array();
    return x_prod_y.sum();
}

// Tạo ma trận xác định dương ngẫu nhiên cấp n
MatrixXd generatePositiveDefiniteMatrix(int n)
{
    // generate a random n x n matrix
    MatrixXd randomMatrix = MatrixXd::Random(n, n);
    MatrixXd transposeMatrix = randomMatrix.transpose();

    // construct a symmetric matrix
    const double symCoef = 0.5;
    MatrixXd symMatrix = symCoef * (transposeMatrix * randomMatrix);

    MatrixXd positiveDefiniteMatrix = symMatrix + n * MatrixXd::Identity(n, n);
    return positiveDefiniteMatrix;
}

double genRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char const *argv[])
{
    freopen("out2.csv", "w", stdout);
    cout << "step,error" << endl;

    srand(int(time(0)));
    int n = 50;

    // Preamble
    // Ma trận B ngẫu nhiên cấp n
    B = MatrixXd::Random(n, n);

    /**
     Ma trận C ngẫu nhiên cấp n
     S = C - C.T nên S là ma trận đối xứng lệch
     */
    C = MatrixXd::Random(n, n);
    S = C - C.transpose();

    /**
     t là vector n chiều ngẫu nhiên
     D là ma trận đường chéo có đường chéo chính là vector t
     */
    VectorXd t = VectorXd::Random(n);
    D = t.array().matrix().asDiagonal();

    // K = B * B.T + S + D
    K = B * B.transpose() + S + D;

    // d là vector ngẫu nhiên n chiều
    d = VectorXd::Random(n);

    // P xác định dương
    P = generatePositiveDefiniteMatrix(n);

    alpha = genRand(1e-4, 1.0 / 3 - 1e-4);
    sigma = genRand(1e-4, 1 - 3 * alpha - 1e-4);
    v = genRand(1e-4, (1 - 3 * alpha - sigma) / (1 - alpha - 2 * alpha * alpha + sigma) - 1e-4);
    r = double(rand());
    L = K.norm();
    M1 = P.norm();
    theta = 0.99;
    gamma_ = 1.5;

    VectorXcd eigenValuesOfP = P.eigenvalues();
    double gtr_max = eigenValuesOfP.real().maxCoeff();
    double gtr_min = eigenValuesOfP.real().minCoeff();
    M2 = (L * sqrt(gtr_max) * (r + sqrt(gtr_min) * d.norm())) / (pow(gtr_min, 1.5) * r);
    M = M1 * M2;

    x0 = VectorXd::Random(n);
    int count = 0;
    double m = x0.norm();
    while (count < 50000)
    {
        cout << count << "," << m << endl;

        double beta_n = genRand(0 + 1e-4, theta / (M + sqrt(M * M + L * L)) - 1e-4);
        VectorXd y = x0 - beta_n * Fx(x0);
        double con_y = Cx(x0) + scalar(P * (x0 - d), -beta_n * Fx(x0));
        if (con_y > 0)
        {
            double lamda = con_y / scalar(P * (x0 - d), P * (x0 - d));
            y = y - lamda * P * (x0 - d);
        }

        VectorXd dx = (x0 - y) - beta_n * (Fx(x0) - Fx(y));
        double lamda_n = (scalar(x0 - y, dx) - M * beta_n * ((x0 - y).norm() * (x0 - y).norm())) / (dx.norm() * dx.norm());
        VectorXd x1 = x0 - gamma_ * lamda_n * dx;
        x0 = x1;
        m = x0.norm();

        count++;
    }
    return 0;
}
