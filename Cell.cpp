// https://www.zybuluo.com/rfhongyi/note/597242 (Simultaneous over-relaxation(SOR))
// https://www.particleincell.com/2015/dielectric-potential-solver/

#include <bits/stdc++.h>
using namespace std;

const double pi = acos(-1.0);
const double eps0 = 8.85e-12;

struct Field2D
{
	int n, m;
	double t, dt;
	double h; // block size (m)
	vector<vector<double>> phi, sigma, epsr, rho;
	// potential, conductivity, relativie permittivity, charge density
	void Init(double l_n, double l_m, double _h, double _dt)
	{
		{ // initialize
			n = (int)ceil(l_n / _h), m = (int)ceil(l_m / _h), h = _h, dt = _dt;
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
		}
		{ // setting
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					double x = i * h, y = j * h - 1e-6;
					if (i > 0 && i < n - 1)
						phi[i][j] = 0.007 / n * i;
					double d = 200e-9; // thickness of the cell membrane
					double R = 4e-6;   // radius of the cell
					double delta = 0;  // distance between the top of the substrate and the bottom of the cell
					if (x > 1.5e-6 && x <= 13.5e-6 && y > 2e-6)
					{
						double r = sqrt((x - 7.5e-6) * (x - 7.5e-6) + (y - 2e-6 - delta) * (y - 2e-6 - delta));
						if ((x > 7.5e-6 - R && x <= 7.5e-6 + R && y <= 2e-6 + d + delta && y > 2e-6 + delta) || (r > R - d && r <= R && y > 2e-6 + delta))
						{
							sigma[i][j] = 0;
							epsr[i][j] = 11;
						}
						else
						{
							sigma[i][j] = 1.5;
							epsr[i][j] = 72;
						}
					}
					else if (x > 0.5e-6 && x <= 14.5e-6 && y > 1e-6)
					{
						sigma[i][j] = 0;
						epsr[i][j] = 4;
					}
					else
					{
						sigma[i][j] = 0;
						epsr[i][j] = 1;
					}
				}
			}
			for (int j = 0; j < m; j++)
			{
				phi[0][j] = 0;
				phi[n - 1][j] = 0.007;
			}
		}
	}
	// type 0: n,m,h,dt,t
	// type 1: phi, type 2: sigma, type 3: epsr, type 4: rho
	// type 5: Electric Field Strength
	void Load(string s, int type)
	{
		FILE *f = fopen(s.c_str(), "r");
		if (type == 0)
		{
			fscanf(f, "%d%d%lf%lf%lf", &n, &m, &h, &dt, &t);
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
		}
		else if (type == 1)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					fscanf(f, "%lf", &phi[i][j]);
		}
		else if (type == 2)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					fscanf(f, "%lf", &sigma[i][j]);
		}
		else if (type == 3)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					fscanf(f, "%lf", &epsr[i][j]);
		}
		else if (type == 4)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					fscanf(f, "%lf", &rho[i][j]);
		}
		fclose(f);
	}
	void LoadAll(string s)
	{
		if (s.back() != '/')
			s.push_back('/');
		for (int i = 0; i <= 4; i++)
		{
			s.push_back('0' + i);
			s += ".txt";
			Load(s, i);
			s.pop_back(), s.pop_back(), s.pop_back(), s.pop_back(), s.pop_back();
		}
	}
	void Save(string s, int type)
	{
		FILE *f = fopen(s.c_str(), "w");
		if (type == 0)
			fprintf(f, "%d %d %E %E %E", n, m, h, dt, t);
		else if (type == 1)
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					fprintf(f, "%E ", phi[i][j]);
				fprintf(f, "\n");
			}
		}
		else if (type == 2)
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					fprintf(f, "%E ", sigma[i][j]);
				fprintf(f, "\n");
			}
		}
		else if (type == 3)
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					fprintf(f, "%E ", epsr[i][j]);
				fprintf(f, "\n");
			}
		}
		else if (type == 4)
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					fprintf(f, "%E ", rho[i][j]);
				fprintf(f, "\n");
			}
		}
		else if (type == 5)
		{
			for (int i = 1; i < n - 1; i++)
			{
				for (int j = 1; j < m - 1; j++)
				{
					double Ex = -(phi[i + 1][j] - phi[i][j]) / h;
					double Ey = -(phi[i][j + 1] - phi[i][j]) / h;
					fprintf(f, "%E ", sqrt(Ex * Ex + Ey * Ey));
				}
				fprintf(f, "\n");
			}
		}
		fclose(f);
	}
	void SaveAll(string s) // Create the folder first!
	{
		if (s.back() != '/')
			s.push_back('/');
		for (int i = 0; i <= 5; i++)
		{
			s.push_back('0' + i);
			s += ".txt";
			Save(s, i);
			s.pop_back(), s.pop_back(), s.pop_back(), s.pop_back(), s.pop_back();
		}
	}
	int SOR(double thres)
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
		} while (_max > thres);
		return cnt;
	}
	void Evolve(string s, double T)
	{
		SaveAll(s);
		int cnt = 0;
		FILE *f = fopen((s + "/6.csv").c_str(), "w+");
		for (; t < T; t += dt)
		{
			if (cnt % 1 == 0)
				SaveAll(s);
			cnt++;
			SOR(1e-7);
			double E = 0; // maximum electric field strength on the cell membrane
			for (int i = 1; i < n - 1; i++)
			{
				for (int j = 1; j < m - 1; j++)
				{
					double x = i * h, y = j * h - 1e-6;
					double d = 200e-9; // thickness of the cell membrane
					double R = 4e-6;   // radius of the cell
					double delta = 0;  // distance between the top of the substrate and the bottom of the cell
					if (x > 1.5e-6 && x <= 13.5e-6 && y > 2e-6)
					{
						double r = sqrt((x - 7.5e-6) * (x - 7.5e-6) + (y - 2e-6 - delta) * (y - 2e-6 - delta));
						if ((x > 7.5e-6 - R && x <= 7.5e-6 + R && y <= 2e-6 + d + delta && y > 2e-6 + delta) || (r > R - d && r <= R && y > 2e-6 + delta))
						{
							double Ex = -(phi[i + 1][j] - phi[i][j]) / h;
							double Ey = -(phi[i][j + 1] - phi[i][j]) / h;
							E = max(E, sqrt(Ex * Ex + Ey * Ey));
						}
					}
				}
			}
			fprintf(f, "%E,%E\n", t, E);
			printf("%d %lf\n", cnt, E);
			for (int i = 1; i < n - 2; i++)
			{
				for (int j = 1; j < m - 2; j++)
				{
					double Ex = -(phi[i + 1][j] - phi[i][j]) / h;
					double Ey = -(phi[i][j + 1] - phi[i][j]) / h;
					rho[i + 1][j] += min(sigma[i + 1][j], sigma[i][j]) * Ex / h * dt;
					rho[i][j] -= min(sigma[i + 1][j], sigma[i][j]) * Ex / h * dt;
					rho[i][j + 1] += min(sigma[i][j + 1], sigma[i][j]) * Ey / h * dt;
					rho[i][j] -= min(sigma[i][j + 1], sigma[i][j]) * Ey / h * dt;
				}
			}
		}
		fclose(f);
	}
} F;

int main()
{
	//F.Init(15e-6, 9e-6, 50e-9, 1e-11);
	F.LoadAll("../EFS/1");
	F.dt = 1e-11;
	F.Evolve("../EFS/2", 1e-7);
	return 0;
}
