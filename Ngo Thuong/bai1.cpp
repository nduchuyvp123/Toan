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
VectorXd x1;
VectorXd x2;

double alpha;
double sigma;
double v;
double r;
double tau;
double alpha_n;
double m;

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
     MatrixXd symMatrix = symCoef * (randomMatrix + transposeMatrix);

     MatrixXd positiveDefiniteMatrix = symMatrix + n * MatrixXd::Identity(n, n);
     return positiveDefiniteMatrix;
}

MatrixXd generatePositiveDefiniteMatrix_2(int n)
{
     MatrixXd P = MatrixXd::Random(n, n);
     MatrixXd P_t = P.transpose();
     MatrixXd PplusP_t = P + P_t;
     P = PplusP_t;
     VectorXcd eigenValuesOfPplusP_t = PplusP_t.eigenvalues();
     double gtr = eigenValuesOfPplusP_t.real().minCoeff() - 1;
     MatrixXd PminusId = P - gtr * MatrixXd::Identity(n, n);
     return PminusId;
}

double genRand(double fMin, double fMax)
{
     double f = (double)rand() / RAND_MAX;
     return fMin + f * (fMax - fMin);
}
int main()
{
     freopen("out1.csv", "w", stdout);
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
     tau = genRand(1e-4, v / K.norm() - 1e-4);

     // x1, x2 ban đầu ngẫu nhiên
     x1 = VectorXd::Random(n);
     x2 = VectorXd::Random(n);

     // Indicamble
     m = x2.norm();
     int count = 0;

     while (count < 50000)
     {
          cout << count << "," << m << endl;

          alpha_n = genRand(1e-4, alpha - 1e-4);

          // w = x_n + alpha * (x_n - x_n-1)
          VectorXd w = x2 + alpha_n * (x2 - x1);

          // y = w - tau * F(w) - ......
          VectorXd y = w - tau * Fx(w);

          // nễu con_y > 0 thì y -= lamda * (P * (w - d))
          double con_y = Cx(w) + scalar(P * (w - d), (-tau) * Fx(w));
          if (con_y > 0)
          {
               double lamda = con_y / scalar(P * (w - d), P * (w - d));
               y = y - lamda * (P * (w - d));
          }

          // x_n+1 = w - tau * Fx(y)
          VectorXd x3 = w - tau * Fx(y);

          // nếu con_x > 0 thì x -= lamda * (P * (w - d))
          double con_x = Cx(w) + scalar(P * (w - d), (-tau) * Fx(y));
          if (con_x > 0)
          {
               double lamda = con_x / scalar(P * (w - d), P * (w - d));
               x3 = x3 - lamda * (P * (w - d));
          }

          x1 = x2;
          x2 = x3;
          m = x2.norm();
          count++;
     }

     return 0;
}
