#ifndef TPOINT_H
#define TPOINT_H

#include "synonymous_types.h"

struct TPoint {
	uint F_idThis; //������� ����� � ��������� Fpoints
	uint F_idCoords;  //������� ������ ���������� ����� ��������� ����� � ��������� Fcoords 
	uint F_idEvaluations; //������� ������ ���������� ����� �������� ������� ����� � ���������
	static uint F_dimension; // ����������� ������
	static uint F_constraints; // ����������� �����������
	static uint F_xxxx; //????
	std::vector<std::list<uint>> inc_coords; //������ ������� �� ���������� ����� � ������� ���������� ���������
	std::vector<std::list<uint>> dec_coords; //������ ������� �� ���������� ����� � ������� ���������� ���������       
};

#endif