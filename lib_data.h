#pragma once
#include<math.h>
#include<iostream>
// ������������
# define at_num	2000
# define sand_num 4000
# define atom_num 6000

// �������ٶȡ�Բ���ʶ���
#define time 0.5;
#define grivaty  9.8
double g_v[3] = { 0, 0, -9.8 };
double a_vt[3] = { -1, 0, -9.8};//��б���������ٶ�
#define pi 3.14159265358979

// �������÷�Χ
#define r_atom 2
#define r_sand 2

// ������
#define block_num 100
#define block_size 1
#define size_ 100

// �˲���
//   z��  �Jy
//    | /
//    |������������������ x
#define x_start 0// ��ʼλ�õ�y����
#define height 80 // �˸�
#define x_width 40 // �˵�x����
#define x_point 10 // �˵�ȡ����
#define x_coast 50 //������λ��
#define init_v  20 // ��ʼ�ٶ�
#define K  0.001  // ���峣��
#define miu_water 0.00009 //ˮ��ճ��ϵ��

// ���ڿ�������ʹ�õĵ�����
struct block
{
	short num = 0;
	short* atom = nullptr;// ʹ��ʱ�ٿ��ؿռ�
};

// ���ڱ�ʾ�������ӵĽṹ

	// ����ȷ����������
	#define sand true
	#define water false

struct bit
{
	bool type;	// ����
	bool peng = false;
	int blockpos; // �ز���λ��

	double p;	// ѹ��
	double lu;	// �ܶ�
	double mess = 0.0006; // ����

	double pos[3] = { 0,0,0 }; // λ��
	double v[3] = { 0,0,0 }; // �ٶ�ʸ��
	double a_[3] = { 0,0,0 };	//���ٶȷ���

	double next_a[3] = { 0 };//�������ü��ٶ�
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