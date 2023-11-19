#pragma once
struct Point
{
	public:
	double X, Y;

	Point():Point(0, 0)
	{
	}

	Point(double x, double y)
	{
		X = x;
		Y = y;
	}
};