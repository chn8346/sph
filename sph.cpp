// sph.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 

#include <iostream>
#include "sph_cal.h"

using namespace sph;
using namespace std;
int main()
{
	std::fstream fp;
	fp.open("data.csv", std::ios::out);
	fp.close();
	fp.open("data2.csv", std::ios::out);
	fp.close();


	init();
	makeUp();
	//makeRain();
	makewall(1,2,3);
	for (int i = 0; i < 10; i++)
	{
		update();
	}
	system("pause");
}