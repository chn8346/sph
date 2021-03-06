#include "StdAfx.h"
#include "SPH.h"
#include <math.h>
#define PI 3.14159265358979323846264338327950L

Particle::Particle()
{
	m_Vel.Set(0,0,0);
	m_VelEval.Set(0,0,0);
	m_Acc.Set(0,0,0);
}

Particle::~Particle()
{

}

CSPH::CSPH(void)
{
	m_smoothLen = 0.01;
	m_timeStep = 0.001;
	m_curStep = 0;
	m_mass = 0.00020543f;
	m_K = 1.5f;
	m_Kpoly6 = 315 * m_mass/((64*PI)*pow(m_smoothLen,9));
	m_Kspiky = 45.0f/(PI*pow(m_smoothLen,6));
	m_Kviscosity = 0.2f;

	m_GlassRange[0] = -0.05;
	m_GlassRange[1] = 0.05;
	m_GlassRange[2] = 0;
	m_GlassRange[3] = 1;
	m_GlassRange[4] = -0.05;
	m_GlassRange[5] = 0.05;

	m_GlassNormal[0].Set(1,0,0);
	m_GlassNormal[1].Set(-1,0,0);
	m_GlassNormal[2].Set(0,1,0);
	m_GlassNormal[3].Set(0,-1,0);
	m_GlassNormal[4].Set(0,0,1);
	m_GlassNormal[5].Set(0,0,-1);

	m_gravity.Set(0,-9.8,0);
	m_iMoveGlassDir = 0;

	ballPos[0] = 0;
	ballPos[1] = 0;
	ballPos[2] = 0;
	ballPos[3] = 0;
	ballPos[4] = 0;
	ballPos[5] = 0;

	fopen_s(&file,"force.txt","w");
	countNum3 = 0;
}


CSPH::~CSPH(void)
{
	if (m_particle)
	{
		delete[] m_particle;
	}

	//�ֿ���
	fclose(file);
}

//��ʼ��
void CSPH::Init(int particleNum)
{
	m_particleNum = particleNum;
	m_particle = new Particle[m_particleNum];
	InitPos();
}

//��ʼ�����ӵ�λ��
void CSPH::InitPos()
{
	float dis = m_smoothLen/**0.6*/; //����֮��ľ���ȡ�⻬�뾶��0.6��ʹ�øտ�ʼ�������õ���
	int cx = int(powf(m_particleNum/4,1.0/3)); 
	int cy = 4*cx;
	int cz = cx;
	for (int i = 0; i < cy; i++)
	{
		for (int j = 0; j < cx; j++)
		{
			for (int k = 0; k < cz; k++)
			{					
				int index = i*cx*cz + j*cz + k;
				if(index>=m_particleNum) 
					return;
				float y = i * dis/* * 0.6*/ + 0.1;
				float x = (j - cx*0.5) * dis;
				float z = (k - cz*0.5) * dis;
				m_particle[index].m_pos.Set(x,y,z);  //��ʼ�����ӵ�λ��		
			}
		}
	}
}

//���߽磬������m_boundMin,m_boundMax��
void CSPH::FindBound()
{
	float minX,minY,minZ;
	float maxX,maxY,maxZ;
	minX = minY = minZ = 10000;
	maxX = maxY = maxZ = -10000;

	//�����������ӵ�λ�ã������С��x��y��z
	for (int i = 0; i < m_particleNum; ++i)
	{
		if (m_particle[i].m_pos.x < minX)
		{
			minX = m_particle[i].m_pos.x;
		}
		if (m_particle[i].m_pos.x > maxX)
		{
			maxX = m_particle[i].m_pos.x;
		}

		if (m_particle[i].m_pos.y < minY)
		{
			minY = m_particle[i].m_pos.y;
		}
		if (m_particle[i].m_pos.y > maxY)
		{
			maxY = m_particle[i].m_pos.y;
		}

		if (m_particle[i].m_pos.z < minZ)
		{
			minZ = m_particle[i].m_pos.z;
		}
		if (m_particle[i].m_pos.z > maxZ)
		{
			maxZ = m_particle[i].m_pos.z;
		}
	}
	m_boundMin.Set(minX,minY,minZ);
	m_boundMax.Set(maxX,maxY,maxZ);
}

void CSPH::FindNeighbour()
{
	//��������ĸ���
	FindBound();
	m_gridNumX = (m_boundMax.x - m_boundMin.x + m_smoothLen)/m_smoothLen;
	m_gridNumY = (m_boundMax.y - m_boundMin.y + m_smoothLen)/m_smoothLen;
	m_gridNumZ = (m_boundMax.z - m_boundMin.z + m_smoothLen)/m_smoothLen;
	m_gridNum = m_gridNumX*m_gridNumY*m_gridNumZ;
	const int gridSaveNum = 100; //ÿ�������д�ŵ�������
	//Ϊ��������ڴ棬�����б���particle������

	Grid = new int[m_gridNum*gridSaveNum];
	GridNum = new int[m_gridNum];
	memset(Grid,0,sizeof(int)*m_gridNum*gridSaveNum);
	memset(GridNum,0,sizeof(int)*m_gridNum);


	for (int i = 0; i < m_particleNum; ++i)
	{
		Vector3DF gridIndexTemp = m_particle[i].m_pos - m_boundMin; //��ǰ������߽�Ĳ�ֵ
		int cx = int(gridIndexTemp.x/m_smoothLen);
		int cy = int(gridIndexTemp.y/m_smoothLen);
		int cz = int(gridIndexTemp.z/m_smoothLen);
		int gridIndex = cx*m_gridNumY*m_gridNumZ + cy*m_gridNumZ + cz;//��ǰ���Ӵ���������
		Grid[gridIndex*gridSaveNum + GridNum[gridIndex]] = i;             //�������¼��������
		m_particleIndexInGrid[i] = gridIndex;
		GridNum[gridIndex]++;                                         //��¼�������Ѿ������˶�������
		if (GridNum[gridIndex] > 100)
		{

			return;
		}
	}
	//�����ھ�
	memset(m_particleNeighbour,0,sizeof(int)*PARTICLE_NUM*200);
	memset(m_particleNeighbourNum,0,sizeof(int)*PARTICLE_NUM);

	for (int i = 0; i < m_particleNum; ++i)
	{
		Vector3DF gridIndexTemp = m_particle[i].m_pos - m_boundMin; //��ǰ������߽�Ĳ�ֵ
		int cx = int(gridIndexTemp.x/m_smoothLen);
		int cy = int(gridIndexTemp.y/m_smoothLen);
		int cz = int(gridIndexTemp.z/m_smoothLen);
		//�����������Ҹ�27������Ԫ
		for (int ix = cx - 1; ix <= cx + 1; ix++)
		{
			for (int iy = cy - 1; iy <= cy + 1; iy++)
			{
				for (int iz = cz - 1; iz <= cz + 1; iz++)
				{
					int gridIndex = ix*m_gridNumY*m_gridNumZ + iy*m_gridNumZ + iz;//��ǰ���Ӵ���������
					if(gridIndex < 0 || gridIndex > m_gridNum)
						continue;
					for (int j = 0; j < GridNum[gridIndex]; ++j)
					{
						int particleId = Grid[gridIndex*gridSaveNum + j]; //���ӵ�����id
						if (particleId == i) //����Ǹ����ӱ���
							continue;
						float disToNeighbour = m_particle[i].m_pos.DistanceTo(m_particle[particleId].m_pos); //���������������ӵľ���
						if (disToNeighbour <= m_smoothLen)
						{
							m_particleNeighbour[i][m_particleNeighbourNum[i]] = particleId;
							m_particleNeighbourDistance[i][m_particleNeighbourNum[i]] = disToNeighbour;
							m_particleNeighbourNum[i]++;
						}
					}
				}
			}
		}
	}
	if (Grid)
	{ 
		delete[] Grid;
	}
	if(GridNum)
	{
		delete[] GridNum;
	}
	
}

//����ÿ�����ӵ��ܶȡ�ѹ��
void CSPH::ComputeDensity()
{
	//�����ܶȹ�ʽ 315/(64*PAI*h^9)Sigma(h^2-r^2)^3	
	//�����315/(64*PAI*h^9)����Ϊ������
	float h6 = pow(m_smoothLen,6);
	
	for (int i = 0; i < m_particleNum; ++i)
	{
		m_particle[i].m_density = 0;
	}

	//
	for (int i = 0; i < m_particleNum; ++i)
	{
		float sum = h6; //�Ȱ��Լ����ܶȼ���
		for (int j = 0; j < m_particleNeighbourNum[i]; ++j)
		{
			float h2_r2 = m_smoothLen*m_smoothLen - m_particleNeighbourDistance[i][j]*m_particleNeighbourDistance[i][j];
			sum += pow(h2_r2,3);
		}

		m_particle[i].m_density = sum * m_Kpoly6;
		m_particle[i].m_pressure = m_K*(m_particle[i].m_density - 1000.0f); //ѹ��=K*(��-��0��
	}
}

void CSPH::ComputeFroce()
{
	for (int i = 0; i < m_particleNum; i++)
	{
		m_particle[i].m_Acc.Set(0,0,0);	
	}
	float h = m_smoothLen;
	Vector3DF vdiff;
	Vector3DF diff;
	
	Vector3DF force1(0,0,0); //ѹ��
	Vector3DF force2(0,0,0); //ճ����

	for (int i = 0; i < m_particleNum; ++i)
	{
		Vector3DF force(0,0,0);
		for (int j = 0; j < m_particleNeighbourNum[i]; ++j)
		{
			float r = m_particleNeighbourDistance[i][j];
			float r_h = h-r; 
			int otherParticleId = m_particleNeighbour[i][j];
			diff = m_particle[i].m_pos - m_particle[otherParticleId].m_pos; //λ�õ�΢��
			vdiff = m_particle[otherParticleId].m_Vel - m_particle[i].m_Vel; //�ٶȵ�΢��
			//ѹ��
			float param = (m_particle[i].m_pressure + m_particle[otherParticleId].m_pressure)*0.5*m_Kspiky*r_h*r_h/ r;
			force1 =diff * (m_mass * param / (m_particle[i].m_density * m_particle[otherParticleId].m_density)); 
			//ճ����
			force2 = vdiff * (m_mass * r_h  * m_Kviscosity * m_Kspiky / (m_particle[i].m_density * m_particle[otherParticleId].m_density));

			force += force1;
			force += force2;
		}	
		m_particle[i].m_Acc = force;
	}

}

void CSPH::DetectCollision(float timeStep)
{
	Vector3DF prePos;
	Vector3DF colAcc;
	for (int i = 0; i < m_particleNum; ++i)
	{
		colAcc.Set(0,0,0);
		prePos = m_particle[i].m_pos + m_particle[i].m_VelEval * timeStep; //���㵱ǰʱ�����ӵ�λ�ã����û����ײӦ���ߵ�����
		GlassCollision(prePos,colAcc,m_particle[i].m_Vel);
		m_particle[i].m_Acc += colAcc;
	}
}


void CSPH::GlassCollision(Vector3DF &prePos, Vector3DF &colAcc, Vector3DF &vel)
{
	float radius = 0.004;
	float EPSILON = 0.00001f;
	float stiff = 30000;
	float damp = 128;
	float r = 0.010; //С��İ뾶

	for (int i = 0; i < 6; i++)
	{
		int j = i / 2;
		float diff;
		if(i%2) 
			diff = 2 * radius - (m_GlassRange[i]- prePos[j]);//1,3,5�������������ʾ�ұ߽硣	���Խ���������������Ǹ�ֵ��diff�ʹ���2*radius
		else 
			diff = 2 * radius - (prePos[j] - m_GlassRange[i]);	//0,2,4��ʾ��߽硣				���Խ���������������Ǹ�ֵ��diff�ʹ���2*radius
		if (diff > EPSILON)//ֻҪ����EPSILON��Ҫִ�У�˵��ǰ���㹻�����߽�ͻ���롣��ǰ���غϱ߽�ʱ��Ӧ���Ѿ����ڵ���2*radius�ˡ�
		{
			float v0 = m_GlassNormal[i]*vel;
			float reverse = stiff*diff - damp*v0;//���ƶ��ٶ�ת��������������	,m_norΪ������ײ��ķ��������ѵ�ǰ���е��ٶ�ת�����������ϡ���ǰ�������ײ�����ֵӦ����<0��
			colAcc += m_GlassNormal[i] * reverse;//һ���Ǿ�������һ�����˶��������أ�������Ӱ������������ײ��ķ������˶���
		}
	}
}


void CSPH::UpdatePos(float timeStep)
{
	for (int i = 0; i < m_particleNum; i++)// ����ģ��    
	{
		Vector3DF finalAcc = m_particle[i].m_Acc + m_gravity;
		Vector3DF velEval = m_particle[i].m_VelEval + finalAcc*timeStep;
		m_particle[i].m_pos += velEval*timeStep;
		m_particle[i].m_Vel = (m_particle[i].m_VelEval + velEval)*0.5;
		m_particle[i].m_VelEval = velEval;
	}
	//MoveGlass();
}

//�ƶ�����
void CSPH::MoveGlass()
{
	if (m_iMoveGlassDir == 0)//��������
	{
		if (m_GlassRange[5] < 0.2)
			m_GlassRange[5] += 0.0005;
		else
			m_iMoveGlassDir = 1;
	}
	else if (m_iMoveGlassDir == 1)
	{
		if (m_GlassRange[5] > 0.01)
			m_GlassRange[5] -= 0.0005;
		else m_iMoveGlassDir = 0;
	}
}

void CSPH::ComputeTouchFroce()
{

};



void CSPH::saveData()
{
	FILE *file;
	fopen_s(&file,"data.txt","a+");
	if (file == NULL)
	{
		return;
	}
	fprintf_s(file,"step:%d\n",m_curStep);
	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 3; j++)
		{
			fprintf_s(file,"%.3d ",m_particle[i].m_Acc[j]);

		}
		fprintf_s(file, "\n");
	}
	fclose(file);
}

	