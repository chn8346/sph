#pragma once
#include "lib_data.h"
#include <fstream>

int bpos(double x, double y, double z)
{
	int a = x;
	int b = y;
	int c = z;
	return(a * block_num * block_num + b * block_num + c);
}

int bpos(int x, int y, int z)
{
	return(x * block_num * block_num + y * block_num + z);
}

int neighbors(int *l, int bpos_)
{
	l = new int[27];
	int x = bpos_ / (block_num * block_num);
	int y = bpos_ % (block_num * block_num)/ block_num;
	int z = bpos_ % block_num;
	int t = 0;
	for (int a = -1; a < 2; a++)
	{
		for (int b = -1; b < 2; b++)
		{
			for (int c = -1; c < 2; c++)
			{
				if (x + a >= 0 && x + a < block_num)
				{
					if (y + b >= 0 && y + b < block_num)
					{
						if (z + c >= 0 && z + c < block_num)
						{
							l[t++] = bpos(x + a, y + b, z + c);
						}
					}
				}
			}
		}
	}
	return t;
}

namespace sph {
	bool is_init = false;
	block* blocks = nullptr;
	bit* atom = nullptr;
	double lu0 = 0;

	void changepos(double x[3], double y[3])
	{
		x[0] = y[0];
		x[1] = y[1];
		x[2] = y[2];
	}
	void changepos(double* x, double y1, double y2, double y3)
	{
		x[0] = y1;
		x[1] = y2;
		x[2] = y3;
	}

	void init()
	{
		// init block
		int blocksize = pow(block_num ,3);
		blocks = new block[blocksize];
		while (--blocksize)
		{
			//blocks[blocksize].atom = new short[at_num];
			blocks[blocksize].num = 0;
		}

		// init atom
		atom = new bit[atom_num];
		int t = at_num;
		while (--t)
		{
			atom[t].peng = false;
			atom[t].type = water;
		}
		t = at_num;
		while (t++ < atom_num)
		{
			atom[t].type = sand;
		}

		// init sand
		for (int i = at_num; i < atom_num; i++)
		{
			atom[i].type = sand;
			changepos(atom[i].v, 0, 0, 0);
			changepos(atom[i].a_, g_v);
		}

		// init lu0
		lu0 = at_num / 2 / block_size / block_num / height;

		is_init = true;
	}

	void makeRain()
	{
		if (!is_init)
		{
			init();
		}
		double rain_height = 80;
		for (int i = 0; i < at_num; i++)
		{
			srand(i);
			double limit = block_num*block_size - rain_height;
			int size = size_;
			atom[i].pos[0] = rand() % (size/2);
			atom[i].pos[0] = atom[i].pos[0] + (size / 2);
			atom[i].pos[1] = rand() % size;
			atom[i].pos[2] = rand();
			atom[i].pos[2] = atom[i].pos[2] / RAND_MAX * limit + rain_height;
			changepos(atom[i].v, 0, 0, -10);
			changepos(atom[i].a_, g_v);
			/*int bsn = bpos(atom[i].pos[0], atom[i].pos[1], atom[i].pos[0]);			
			blocks[bsn].atom[blocks[bsn].num] = i;
			blocks[bsn].num++;*/
		}
	}

	void makewall(int ytan, int xtan, int h)
	{
		for (int i = at_num; i < atom_num; i++)
		{
			double r = rand();
			atom[i].pos[0] = 50 + r / RAND_MAX * 50;
			r = rand();
			atom[i].pos[1] = r / RAND_MAX * 100;
			r = rand();
			atom[i].pos[2] = r / RAND_MAX * 60;
		}
	}

	void makeUp()
	{
		if (!is_init)
		{
			init();
		}

		// 造浪 暂时用正弦代替
		int chipNum = at_num / block_num;
		int step = x_width / x_point;
		int point_num = chipNum / x_point;
		int sn = 0;
		for (int i = 0; i < block_num; i++)
		{
			for (int j = 0; j < x_point; j++)
			{
				double dv = cos(pi / x_point * j);
				double dv2 = sin(pi / x_point * j);
				double part = dv2 / point_num;
				for (int k = 0; k < point_num; k++)
				{
					changepos(atom[sn].a_, a_vt);
					double r = rand();
					changepos(atom[sn].v, init_v, r/RAND_MAX*init_v/2, 0);
					changepos(atom[sn].pos,x_start+j*step ,i, part*k*height);
					/*
					int bsn = bpos(atom[sn].pos[0], atom[sn].pos[1], atom[sn].pos[2]);
					blocks[bsn].atom[blocks[bsn].num] = sn;
					blocks[bsn].num++;
					atom[sn].blockpos = bsn;
					sn++;
					*/
					sn++;
				}
			}
		}
		std::fstream fp;
		fp.open("data.csv", std::ios::app);

		for (int o = 0; o < at_num;o++)
		{
			for (int i = 0; i < 3; i++)
			{
				fp << atom[o].pos[i] << ",";
			}
			for (int i = 0; i < 3; i++)
			{
				fp << atom[o].v[i] << ",";
			}
			for (int i = 0; i < 3; i++)
			{
				fp << atom[o].a_[i] << ",";
			}
			fp << "\n";
		}
		fp.close();
	}

	void update_lu()
	{
		for (int i = 0; i < at_num; i++)
		{
			std::cout << "update lu " << i << " ";
			double h = r_atom;
			double new_lu = 0;
			for (int j = 0; j < at_num; j++)
			{
				
				double t[3];
				for (int k = 0; k < 3; k++)
				{
					double t1 = atom[j].pos[k];
					t[k] = atom[i].pos[k] - t1;
				}

				if (h * h <= realvalue2(3, t))
					continue;

				new_lu = new_lu + powl(h*h - realvalue2(3, t), 1.0/3);
			}
			if (atom[i].pos[1] < r_atom)
			{
				for (int j = 0; j < at_num; j++)
				{
					double t[3];
					for (int k = 0; k < 3; k++)
					{
						double t1 = atom[j].pos[k] + block_num * block_size;
						t[k] = atom[i].pos[k] - t1;
					}
					if (h * h <= realvalue2(3, t))
						continue;

					new_lu = new_lu + powl(h * h - realvalue2(3, t), 1.0 / 3);
				}
			}
			if (atom[i].pos[1] + r_atom > block_num* block_size)
			{
				for (int j = 0; j < at_num; j++)
				{
					double t[3];
					for (int k = 0; k < 3; k++)
					{
						double t1 = atom[j].pos[k] - block_num * block_size;
						t[k] = atom[i].pos[k] - t1;
					}
					if (h * h <= realvalue2(3, t))
						continue;

					new_lu = new_lu + powl(h * h - realvalue2(3, t), 1.0 / 3);
				}
			}
			atom[i].lu = new_lu*315*(atom[i].mess)/64/pi/powl(h, 6);
			std::cout << atom[i].lu << std::endl;
		}
	}

	void update_press()
	{
		double h = r_atom;
		for (int i = 0; i < at_num; i++)
		{
			atom[i].p = K * (atom[i].lu - lu0);
		}

		for (int i = 0; i < at_num; i++)
		{
			std::cout << "update press " << i << " ";
			double new_a[3] = { 0,0,0 };
			for (int j = 0; j < at_num; j++)
			{
				if (i == j)
					continue;
				double t[3];
				for (int k = 0; k < 3; k++)
				{
					t[k] = atom[i].pos[k] - atom[j].pos[k];
				}

				if (h * h < realvalue2(3, t))
					continue;

				double real = realvalue(3, t);
				if (real == 0)
					continue;
				double ac = (atom[i].p + atom[j].p) * powl((h - real), 2) / atom[i].lu / atom[j].lu / 6 / real;
				for (int k = 0; k < 3; k++)
				{
					new_a[k] = new_a[k] + ac * t[k];
				}
			}
			if (atom[i].pos[1] < r_atom)
			{
				for (int j = 0; j < at_num; j++)
				{
					if (i == j)
						continue;
					double t[3];
					for (int k = 0; k < 3; k++)
					{
						t[k] = atom[i].pos[k] - atom[j].pos[k];
					}
					t[1] = t[1] + block_num * block_size;

					if (h * h < realvalue2(3, t))
						continue;

					double real = realvalue(3, t);
					if (real == 0)
						continue;
					double ac = (atom[i].p + atom[j].p) * powl((h - real), 2) / atom[i].lu / atom[j].lu / 6 / real;
					for (int k = 0; k < 3; k++)
					{
						new_a[k] = new_a[k] + ac * t[k];
					}
				}
			}
			if (atom[i].pos[1] + r_atom > block_num* block_size)
			{
				for (int j = 0; j < at_num; j++)
				{
					if (i == j)
						continue;
					double t[3];
					for (int k = 0; k < 3; k++)
					{
						t[k] = atom[i].pos[k] - atom[j].pos[k];
					}
					t[1] = t[1] - block_num * block_size;

					if (h * h < realvalue2(3, t))
						continue;

					double real = realvalue(3, t);
					if (real == 0)
						continue;
					double ac = (atom[i].p + atom[j].p) * powl((h - real), 2) / atom[i].lu / atom[j].lu / 6 / real;
					for (int k = 0; k < 3; k++)
					{
						new_a[k] = new_a[k] + ac * t[k];
					}
				}
			}
			double ac = atom[i].mess / pi / powl(h, 6);
			for (int k = 0; k < 3; k++)
			{
				atom[i].next_a[k] = -new_a[k] * ac;
			}
			std::cout << atom[i].next_a[0] << atom[i].next_a[1] << atom[i].next_a[2] << std::endl;
		}
	}

	void update_stick()
	{
		double h = r_atom;
		for (int i = 0; i < at_num; i++)
		{
			std::cout << "update stick" << i ;
			double new_sta[3] = { 0,0,0 };
			for (int j = 0; j < at_num; j++)
			{
				if (i == j)
					continue;
				double t[3];
				for (int k = 0; k < 3; k++)
				{
					t[k] = atom[i].pos[k] - atom[j].pos[k];
				}

				if (h * h < realvalue2(3, t))
					continue;
			
				double ac = (h - realvalue(3, t)) / atom[i].lu / atom[j].lu;

				for (int k = 0; k < 3; k++)
				{
					new_sta[k] = new_sta[k] + ac * (atom[i].v[k] - atom[j].v[k]);
				}
			}

			if (atom[i].pos[1] < r_atom)
			{
				for (int j = 0; j < at_num; j++)
				{
					if (i == j)
						continue;
					double t[3];
					for (int k = 0; k < 3; k++)
					{
						t[k] = atom[i].pos[k] - atom[j].pos[k];
					}
					t[1] = t[1] + block_num * block_size;
					if (h * h < realvalue2(3, t))
						continue;

					double ac = (h - realvalue(3, t)) / atom[i].lu / atom[j].lu;

					for (int k = 0; k < 3; k++)
					{
						new_sta[k] = new_sta[k] + ac * (atom[i].v[k] - atom[j].v[k]);
					}
				}
			}
			if (atom[i].pos[1] + r_atom > block_num* block_size)
			{
				for (int j = 0; j < at_num; j++)
				{
					if (i == j)
						continue;
					double t[3];
					for (int k = 0; k < 3; k++)
					{
						t[k] = atom[i].pos[k] - atom[j].pos[k];
					}
					t[1] = t[1] - block_num * block_size;
					if (h * h < realvalue2(3, t))
						continue;

					double ac = (h - realvalue(3, t)) / atom[i].lu / atom[j].lu;

					for (int k = 0; k < 3; k++)
					{
						new_sta[k] = new_sta[k] + ac * (atom[i].v[k] - atom[j].v[k]);
					}
				}
			}

			double ac = (atom[i].mess) * miu_water * 5 / pi / powl(h, 6);
			for (int k = 0; k < 3l; k++)
			{
				atom[i].next_a[k] = atom[i].next_a[k] + ac * new_sta[k];
			}
			std::cout << atom[i].next_a[0] << atom[i].next_a[1] << atom[i].next_a[2] << std::endl;

		}
	}

	void update()
	{
		update_lu();
		update_press();
		update_stick();

		for (int p = 0; p < at_num; p++)
		{
			std::cout << p << std::endl;
			// pos update
			for (int i = 0; i < 3; i++)
			{
				atom[p].pos[i] = atom[p].pos[i] + atom[p].v[i]*time;
			}

			// v update
			for (int i = 0; i < 3; i++)
			{
				atom[p].v[i] = atom[p].v[i] + atom[p].a_[i]*time;
			}

			// a update
			for (int i = 0; i < 3; i++)
			{
				atom[p].a_[i] = atom[p].next_a[i] + a_vt[i];
			}

			// pos adjust
			if (atom[p].pos[0] > block_num*block_size)
			{
				atom[p].pos[0] = block_num * block_size;
				atom[p].a_[0] = -atom[p].a_[0];
				atom[p].v[0] = -atom[p].v[0];
			}

			if (atom[p].pos[0] < 0)
			{
				atom[p].pos[0] = 0;
				atom[p].a_[2] = 0;
				atom[p].v[2] = 10;
			}

			if (atom[p].pos[1] <= 0)
				atom[p].pos[1] = block_num * block_size;

			if (atom[p].pos[1] >= block_num* block_size)
			{
				atom[p].pos[1] = 0;
			}

			if (atom[p].pos[2] < 0)
			{
				atom[p].pos[2] = 0;
				atom[p].a_[2] = -atom[p].a_[2]*0.1;
				atom[p].v[2] = -atom[p].v[2]*0.1;
			}

			if (atom[p].pos[2] > block_num* block_size)
			{
				atom[p].pos[2] = block_num * block_size;
				atom[p].a_[2] = -atom[p].a_[2];
				atom[p].v[2] = -atom[p].v[2];
			}

			double h = 5;
			for (int i = 0; i < at_num; i++)
			{
				if (!atom[i].peng)
				{
					for (int j = at_num; i < atom_num; i++)
					{
						double t[3];
						for (int k = 0; k < 3; k++)
						{
							t[k] = atom[i].pos[k] - atom[j].pos[k];
						}


						if (realvalue2(3, t) < h * h)
						{
							//atom[i].peng = true;
							changepos(atom[i].v, 0, 0, 0);
							changepos(atom[i].a_, a_vt);
							atom[j].pos[0] = atom[j].pos[0] + t[0];
							atom[j].pos[1] = atom[j].pos[1] + t[1];
							atom[j].pos[2] = atom[j].pos[2] + t[2];
						}
					}
				}
			}
		}

		// file write
		std::fstream fp;
		fp.open("data.csv", std::ios::app);
		for (int o = 0; o < at_num; o++)
		{
			for (int i = 0; i < 3; i++)
			{
				fp << atom[o].pos[i] << ",";
			}
			for (int i = 0; i < 3; i++)
			{
				fp << atom[o].v[i] << ",";
			}
			for (int i = 0; i < 3; i++)
			{
				fp << atom[o].a_[i] << ",";
			}
			fp << "\n";
		}
		fp.close();

		fp.open("data2.csv", std::ios::app);
		for (int o = at_num; o < atom_num; o++)
		{
			for (int i = 0; i < 3; i++)
				fp << atom[o].pos[i] << ",";
			fp << "\n";
		}
		fp.close();
	}
}
