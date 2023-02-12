// https://www.zybuluo.com/rfhongyi/note/597242 (Simultaneous over-relaxation(SOR))
// https://www.particleincell.com/2015/dielectric-potential-solver/

#include <bits/stdc++.h>
using namespace std;

const double pi = acos(-1.0);
const double eps0 = 8.85e-12;

struct Field2D
{
	int n, m;
	double dt;
	double h; // block size (m)
	vector<vector<double>> phi, sigma, epsr, rho;
	// potential, conductivity, relativie permittivity, charge density

	void Init(double l_n, double l_m, double _h, double _dt)
	{
		n = int(l_n / _h), m = int(l_m / _h), h = _h, dt = _dt;
		phi.clear();
		sigma.clear();
		epsr.clear();
		rho.clear();
		vector<double> v(m);
		for (int i = 0; i < m; i++)
			v[i] = 0;
		for (int i = 0; i < n; i++)
		{
			phi.push_back(v);
			sigma.push_back(v);
			epsr.push_back(v);
			rho.push_back(v);
		}
		{ // setting
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					double x = i * h, y = j * h;
					if ((x - 10) * (x - 10) + (y - 5) * (y - 5) > 2 * 2)
					{
						sigma[i][j] = 0;
						epsr[i][j] = 1;
					}
					else
					{
						sigma[i][j] = 0;
						epsr[i][j] = 3;
					}
				}
			}
			for (int j = 0; j < m; j++)
			{
				phi[0][j] = -30;
				phi[n - 1][j] = 30;
			}
		}
	}
	int SOR()
	{
		int cnt = 0;
		double _max;
		do
		{
			_max = 0;
			for (int i = 1; i < n - 1; i++)
			{
				phi[i][0] = phi[i][1];
				phi[i][m - 1] = phi[i][m - 2];
			}
			for (int i = 1; i < n - 1; i++)
			{
				for (int j = 1; j < m - 1; j++)
				{
					double phi_new = (epsr[i + 1][j] * phi[i + 1][j] + epsr[i][j] * phi[i - 1][j] +
									  epsr[i][j + 1] * phi[i][j + 1] + epsr[i][j] * phi[i][j - 1] + rho[i][j] * h * h / eps0) /
									 (epsr[i + 1][j] + epsr[i][j + 1] + 2 * epsr[i][j]);
					double dphi = phi_new - phi[i][j];
					_max = max(_max, fabs(dphi));
					phi[i][j] += dphi * 2 / (1 + pi / (n * m));
				}
			}
			cnt++;
		} while (_max > 0.001);
		return cnt;
	}
	void Print()
	{
		FILE *f = fopen("../f.txt", "w");
		for (int i = 0; i < n - 1; i++)
		{
			for (int j = 0; j < m - 1; j++)
			{
				double Ex = (phi[i + 1][j] - phi[i][j]) / h;
				double Ey = (phi[i][j + 1] - phi[i][j]) / h;
				fprintf(f, "%lf ", sqrt(Ex * Ex + Ey * Ey));
			}
			fprintf(f, "\n");
		}
	}
} F;

int main()
{
	F.Init(20, 10, 0.1, 1);
	cout << F.SOR() << endl;
	F.Print();
	return 0;
}
