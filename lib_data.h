#pragma once
#include<math.h>
#include<iostream>
// 粒子数量定义
# define at_num	2000
# define sand_num 4000
# define atom_num 6000

// 重力加速度、圆周率定义
#define time 0.5;
#define grivaty  9.8
double g_v[3] = { 0, 0, -9.8 };
double a_vt[3] = { -1, 0, -9.8};//倾斜的重力加速度
#define pi 3.14159265358979

// 粒子作用范围
#define r_atom 2
#define r_sand 2

// 区块量
#define block_num 100
#define block_size 1
#define size_ 100

// 浪参数
//   z↑  ↗y
//    | /
//    |————————→ x
#define x_start 0// 开始位置的y坐标
#define height 80 // 浪高
#define x_width 40 // 浪的x轴宽度
#define x_point 10 // 浪的取点数
#define x_coast 50 //海岸线位置
#define init_v  20 // 初始速度
#define K  0.001  // 流体常量
#define miu_water 0.00009 //水的粘度系数

// 用于快速搜索使用的单个块
struct block
{
	short num = 0;
	short* atom = nullptr;// 使用时再开拓空间
};

// 用于表示单个粒子的结构

	// 用于确定粒子类型
	#define sand true
	#define water false

struct bit
{
	bool type;	// 类型
	bool peng = false;
	int blockpos; // 回查块的位置

	double p;	// 压力
	double lu;	// 密度
	double mess = 0.0006; // 质量

	double pos[3] = { 0,0,0 }; // 位置
	double v[3] = { 0,0,0 }; // 速度矢量
	double a_[3] = { 0,0,0 };	//加速度方向

	double next_a[3] = { 0 };//计算所得加速度
};

double realvalue2(int n, double* x)
{
	double x2 = 0;
	for (int i = 0; i < n; i++)
	{
		x2 = x2 + x[i] * x[i];
		if (x2 < 0)
			return -1;
	}
	return x2;
}

double realvalue(int n, double* x)
{
	double x2 = realvalue2(n, x);
	return sqrt(x2);
}