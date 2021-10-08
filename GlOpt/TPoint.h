#ifndef TPOINT_H
#define TPOINT_H

#include "synonymous_types.h"

struct TPoint {
	uint F_idThis; //ѕозици€ точки в хранилище Fpoints
	static uint F_dimension; // размерность задачи
	static uint F_constraints; // размерность ограничений
	// static uint F_xxxx; //????
	std::vector<std::list<uint>> inc_coords; //¬ектор списков на порождЄнные точки в сторону увеличени€ координат
	std::vector<std::list<uint>> dec_coords; //¬ектор списков на порождЄнные точки в сторону уменьшени€ координат    

	uint get_id_coord() const { return F_idThis * F_dimension; }
	uint get_id_evaluations() const { return F_idThis * (F_constraints + 1); }

	enum class direction { DECREASE, INCREASE };
	enum class TPointTypeError { OK, ERROR };

	TPoint() : F_idThis(0) {
		if (F_dimension == 0) TPointTypeError::ERROR;
		else TPointTypeError::OK;
		dec_coords.resize(F_dimension);
		inc_coords.resize(F_dimension);
	}

	uint does_point_exist(const uint& des_val, 
						  const direction& dir,
						  const uint& divide_axis, 
						  const std::deque<uint>& coords) {
		if (dir == direction::INCREASE) {
			auto beg = inc_coords[divide_axis].begin();
			auto end = inc_coords[divide_axis].end();
			for (auto it = beg; it != end; ++it)
				if (coords[(*it) * F_dimension + divide_axis] == des_val) return *it;
		}
		else {
			auto beg = dec_coords[divide_axis].begin();
			auto end = dec_coords[divide_axis].end();
			for (auto it = beg; it != end; ++it)
				if (coords[(*it) * F_dimension + divide_axis] == des_val) return *it;
		}

		return 0;
	}

	void connect_points(const uint& id_child, 
		                const uint& divide_axis, 
		                const TPoint::direction& dir) {
		if (dir == direction::DECREASE)
			dec_coords[divide_axis].push_back(id_child);
		else
			inc_coords[divide_axis].push_back(id_child);
	}
};

#endif