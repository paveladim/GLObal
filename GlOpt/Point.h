#ifndef POINT_H
#define POINT_H

#include "synonymous_types.h"

class Point {
	// идентификатор точки
	uint _id;

	//¬ектор списков на порождЄнные точки в сторону увеличени€ координат
	std::vector<std::list<uint>> _inc_coord;
	//¬ектор списков на порождЄнные точки в сторону уменьшени€ координат
	std::vector<std::list<uint>> _dec_coord;

	// параметры задачи
	static uint _dimension;
	static uint _constraints;
public:
	// направлени€ изменени€ координат
	enum class Direction { DECREASE, INCREASE };
public:
	Point();
	uint get_id() const;
	uint get_id_coord() const;
	uint get_id_evaluations() const;
	void set_id(const uint& id);
	uint does_point_exist(const uint&,
						  const Direction&,
						  const uint&,
						  const std::deque<uint>&);
	void connect_points(const uint&,
						const uint&,
						const Direction&);
};

#endif // POINT_H

