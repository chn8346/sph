#pragma once
#include "Vector3DF.h"
#include "MFC_OSG.h"

#define PARTICLE_NUM 4000


class Particle
{
public:
	Particle();
	~Particle();
	Vector3DF m_pos;   //λ��
	double color[3];      //��ɫ

	float m_pressure;  //ѹ��
	float m_density;   //�ܶ�

	Vector3DF m_Vel;   //�ٶ�
	Vector3DF m_VelEval;  //Ԥ�����ٶ�
	Vector3DF m_Acc;   //���ٶ�

};


class CSPH
{
public:
	CSPH(void);
	~CSPH(void);
	void Init(int particleNum); //��ʼ�� ָ�����ӵĸ�������ĸ���, �Լ���������λ�õĳ�ʼ��
	void InitPos();				//��ʼ�����ӵ�λ��
	void FindBound();			//�������ı߽磬�������ӵ����ֵ����Сֵ���
	void FindNeighbour();		//����ÿһ�����ӵ���������

	void ComputeDensity();			//�õ����ӵ��ܶ�
	void ComputeFroce();			//�������ӵ�����
	void DetectCollision(float timeStep);			//��ײ���
	void GlassCollision(Vector3DF &prePos, Vector3DF &colPos, Vector3DF &vel);			//��ײ���
	
	void UpdatePos(float timeStep); //�������ӵ�λ��
	void MoveGlass();				//�ƶ�����
	void ComputeTouchFroce();		//����Ӵ���

	void saveData(); //��������
	

	Particle *m_particle; //����
	int m_particleNum;

	float m_K;			//�ܶ�ת��ѹ����ϵ��
	float m_mass;       //��������
	float m_Kpoly6;     //�ܶȺ˺���ϵ��
	float m_Kspiky;     //ѹ���˺���ϵ��
	float m_Kviscosity; //ճ�Ȳ��� 0.2f

	Vector3DF m_gravity;  //����

	//
	float m_smoothLen;//�⻬�˺����뾶h
	float m_timeStep; //ʱ�䲽�� ʹ����һ���̶����� 0.001
	int m_curStep;	  //��ǰ����


	//����
	Vector3DF m_boundMin; //�±߽�
	Vector3DF m_boundMax; //�ϱ߽�
	int m_gridNumX,m_gridNumY,m_gridNumZ; //����Ĵ�С��x*y*z; ��������λ�ü���
	int m_gridNum;
	int *Grid; //���������д�ŵ����ӵ�������ÿһ�����������һ�ٸ� new [x][y][z][100]
	int *GridNum;   //���������Ӧ�����д�����ӵĸ��� new [x][y][z]

	//
	int m_particleIndexInGrid[PARTICLE_NUM]; //�����������е�λ��
	int m_particleNeighbour[PARTICLE_NUM][100]; //�����ھӵ�id;
	int m_particleNeighbourNum[PARTICLE_NUM];   //�����ھӵ�����;
	float m_particleNeighbourDistance[PARTICLE_NUM][100];  //�����ھӵľ���

	//����
	float m_GlassRange[6]; //x(-+) y(-+) z(-+)
	Vector3DF m_GlassNormal[6];  //�������Ӧ�ķ�����
	int m_iMoveGlassDir;
	double ballPos[6]; //С���λ��
	double ballPosPer[3];

	//
	FILE *file;
	int countNum3;
};

