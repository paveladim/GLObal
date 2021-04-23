#ifndef TPOINT_H
#define TPOINT_H

#include "synonymous_types.h"

struct TPoint {
	uint F_idThis; //ѕозици€ точки в хранилище Fpoints
	uint F_idCoords;  //ѕозици€ начала размещени€ блока координат точки в хранилище Fcoords 
	uint F_idEvaluations; //ѕозици€ начала размещени€ блока значений функций точки в хранилище
	static uint F_dimension; // размерность задачи
	static uint F_constraints; // размерность ограничений
	static uint F_xxxx; //????
	std::vector<std::list<uint>> inc_coords; //¬ектор списков на порождЄнные точки в сторону увеличени€ координат
	std::vector<std::list<uint>> dec_coords; //¬ектор списков на порождЄнные точки в сторону уменьшени€ координат       
};

#endif